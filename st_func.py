import streamlit as st

import pandas as pd
import os
import time
from Bio import SeqIO
import requests
from bs4 import BeautifulSoup

from cross_allergen import *
from parser import *

def search_allergen_org_with_real_uniprot_st(source_keyword="peanut", fasta_dir="./fasta_files"):
    base_url = "https://www.allergen.org"
    search_url = f"{base_url}/search.php"
    params = {
        "allergenname": "",
        "allergensource": source_keyword,
        "TaxSource": "",
        "TaxOrder": "",
        "foodallerg": "all",
        "bioname": "",
    }

    with st.spinner(f"🔍 正在搜索 {source_keyword} 相关过敏原..."):
        resp = requests.get(search_url, params=params)
        if resp.status_code != 200:
            st.error(f"❌ 请求失败: {resp.status_code}")
            return []

    soup = BeautifulSoup(resp.text, "html.parser")
    table = soup.find("table", {"border": "1"})
    if not table:
        st.warning("⚠️ 未找到结果表格")
        return []

    rows = table.find_all("tr")[2:]
    total_allergens = len(rows)
    st.info(f"共找到 {total_allergens} 个过敏原条目")

    os.makedirs(fasta_dir, exist_ok=True)
    results = []

    progress_bar = st.progress(0)  # 进度条
    status_text = st.empty()       # 动态文本

    missing_ids = []

    for idx, row in enumerate(rows, 1):
        cols = row.find_all("td")
        try:
            allergen_name = cols[1].text.strip()
            biological_name = cols[2].text.strip()
            iuis_id = cols[1].find("a")["href"].split("=")[-1].strip()
            uniprot_id = get_uniprot_id_from_detail_page(iuis_id)
        except:
            uniprot_id = None

        sequence = None
        if uniprot_id:
            fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            fasta_resp = requests.get(fasta_url)
            if fasta_resp.ok:
                sequence = fasta_resp.text
                fasta_path = os.path.join(fasta_dir, f"{uniprot_id}.fasta")
                with open(fasta_path, "w", encoding="utf-8") as f:
                    f.write(sequence)
            else:
                missing_ids.append(uniprot_id)
        else:
            missing_ids.append(allergen_name)

        results.append({
            "allergen_name": allergen_name,
            "biological_name": biological_name,
            "iuis_id": iuis_id,
            "uniprot_id": uniprot_id,
            "sequence": sequence
        })

        progress_bar.progress(idx / total_allergens)
        status_text.text(f"[{idx}/{total_allergens}] 已获取: {allergen_name}")
    if missing_ids:
        st.warning(f"⚠️ 未找到 UniProt ID 的过敏原: {', '.join(missing_ids)}")

    return results

def predict_cross_allergen_streamlit(species_name, query_fasta_dir="fasta_files", db_fasta="./naive/uniprotkb_allergen_2025_07_22.fasta", k=5, cpu=None, local_enabled=True):

    # 1. 读取缓存的 CSV
    species_csv = f"species_cache/{species_name}_allergens.csv"
    if not os.path.exists(species_csv):
        st.error(f"❌ 未找到缓存文件: {species_csv}，请先运行物种过敏原搜索。")
        return

    allergen_df = pd.read_csv(species_csv)
    if "uniprot_id" not in allergen_df.columns:
        st.error("❌ CSV 文件缺少 'uniprot_id' 列，无法继续。")
        return

    allergen_ids = allergen_df["uniprot_id"].tolist()
    st.info(f"🔍 从缓存读取到 {len(allergen_ids)} 个 {species_name} 过敏原 ID")

    # 2. 加载数据库序列
    db_records = list(SeqIO.parse(db_fasta, "fasta"))

    # 3. 找到对应的 fasta 文件
    query_files = []
    for uid in allergen_ids:
        fasta_path = os.path.join(query_fasta_dir, f"{uid}.fasta")
        if os.path.exists(fasta_path):
            query_files.append(fasta_path)

    if not query_files:
        st.error("❌ 没有找到可用的 query FASTA 文件")
        return

    st.info(f"📂 找到 {len(query_files)} 条待测序列")

    # 4. 运行检测
    all_results = []
    total_start = time.time()
    log_box = st.empty()
    log_text = ""

    for idx, qfile in enumerate(query_files, 1):
        query_seq = SeqIO.read(qfile, "fasta")
        log_text += f"\n[{idx}/{len(query_files)}] 处理序列: {query_seq.id} ({qfile})"
        log_box.code(log_text)

        start = time.time()
        results = check_cross_parallel(query_seq, db_records, win=80, id_thresh=0.35, n_jobs=cpu, local_enabled=local_enabled)
        end = time.time()

        if results:
            results = dedup_by_os(results)
            results.sort(key=lambda x: x["score"], reverse=True)
            for r in results:
                r['query_id'] = query_seq.id
            all_results.extend(results)
            log_text += f"\n  ✔️ 找到 {len(results)} 个潜在交叉过敏原，耗时 {round(end-start, 2)} s"
        else:
            log_text += "\n  ⚠️ 未发现潜在交叉过敏原"

        log_box.code(log_text)

    total_end = time.time()

    if not all_results:
        st.warning("❌ 全部序列均未发现潜在交叉过敏原")
        return

    # 5. 保存 txt
    os.makedirs("result_cache", exist_ok=True)
    result_txt = f"result_cache/{species_name}.txt"
    with open(result_txt, 'w', encoding='utf-8') as f:
        current_query = None
        for i, r in enumerate(all_results, 1):
            if r['query_id'] != current_query:
                current_query = r['query_id']
                f.write(f"\n\n==== Query: {current_query} ====\n")
            f.write(f"{i}. {r['id']}\n"
                    f"   📛 名称: {r['description']}\n"
                    f"   🧬 Identity: {r['identity']:.3f}\n"
                    f"   🔗 Has 6-mer match: {r['has_6mer']}\n"
                    f"   📊 综合得分: {r['score']:.3f}\n"
                    f"   🔍 匹配模式: {r['match_mode']}\n")

    # 6. 转 DataFrame
    content = {key: [] for key in all_results[0].keys()}
    for item in all_results:
        for key in content.keys():
            if key == 'id':
                item[key] = item[key].split('|')[1]
            content[key].append(item[key])
    df = pd.DataFrame(content)

    st.success(f"🎉 所有查询序列处理完成，总耗时: {round(total_end - total_start, 2)} s")

    return df
