#!/usr/bin/env python3
import os
from pickle import FALSE
import re
import sys
import time
import json
import edlib
import argparse
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

with open('data/epitope_map.json', encoding='utf-8') as f:
    epitope_map = json.load(f)

with open('data/epitope_info.json', encoding='utf-8') as f:
    epitope_info = json.load(f)

def sliding_window_identity_fast(seq1, seq2, win=80):
    best = 0.0
    for i in range(len(seq1) - win + 1):
        sub = seq1[i:i+win]
        result = edlib.align(sub, seq2, mode='HW', task='distance')
        if result["editDistance"] == -1:
            continue
        identity = 1.0 - result["editDistance"] / win
        if identity > best:
            best = identity
        if best == 1.0:
            return best
    return best


def sliding_window_identity(seq1, seq2, win=80):
    best = 0.0
    for i in range(len(seq1) - win + 1):
        sub = seq1[i:i+win]
        score = pairwise2.align.globalxx(sub, seq2, one_alignment_only=True, score_only=True)
        identity = score / win
        if identity > best:
            best = identity
    return best


def has_exact8mer(seq1, seq2):
    sixmers = {seq1[i:i+8] for i in range(len(seq1) - 8 + 1)}
    for mer in sixmers:
        if mer in seq2:
            return True
    return False

def best_epitope_identity_edlib(query_str, epitope_list):
    best_identity = 0.0
    matched_epitopes = []

    for epitope in epitope_list:
        if len(epitope) < 6:
            continue
        result = edlib.align(epitope, query_str, mode='HW', task='distance')
        if result["editDistance"] == -1:
            continue
        identity = 1.0 - result["editDistance"] / len(epitope)
        if identity > best_identity:
            best_identity = identity
            matched_epitopes = [epitope]
        elif identity == best_identity:
            matched_epitopes.append(epitope)

        if identity == 1.0:
            return 1.0, [epitope]

    return best_identity, matched_epitopes


def count_substring_hits(subseq, target_str):
    """统计subseq在target_str中出现了多少次（允许重叠）"""
    count = start = 0
    while True:
        start = target_str.find(subseq, start)
        if start == -1:
            break
        count += 1
        start += 1  # 允许重叠匹配
    return count

def check_one_target(args, exact6_bonus=0.1, epitope_factor=0.8):
    query_str, rec, win, id_thresh, local_enabled = args
    target_str = str(rec.seq).upper()
    if target_str == query_str:
        return None
    
    uniprot_id = rec.id.split("|")[1]

    best_identity = 0.0
    best_start = -1
    for i in range(len(query_str) - win + 1):
        sub = query_str[i:i+win]
        result = edlib.align(sub, target_str, mode='HW', task='distance')
        if result["editDistance"] == -1:
            continue
        identity = 1.0 - result["editDistance"] / win
        if identity > best_identity:
            best_identity = identity
            best_start = i
        if best_identity == 1.0:
            break

    if best_identity < id_thresh:
        return None

    best_window_seq = query_str[best_start:best_start+win]

    if local_enabled:
        # epitope列表
        epitope_list = epitope_map.get(uniprot_id, [])
        # 判断该window是否包含epitope
        window_contains_epitope = any(epitope in best_window_seq for epitope in epitope_list)
        # 调整identity
        adjusted_identity = best_identity if window_contains_epitope else best_identity * epitope_factor
    else:
        window_contains_epitope = False
        adjusted_identity = best_identity

    # 统计该window在目标序列hit次数
    window_hit_count = count_substring_hits(best_window_seq, target_str)

    # 统计所有8mer在目标序列的精准匹配次数
    eightmers = {query_str[i:i+8] for i in range(len(query_str) - 8 + 1)}
    eightmer_hit_count = 0
    for mer in eightmers:
        eightmer_hit_count += count_substring_hits(mer, target_str)

    # 是否有6mer匹配
    exact6 = has_exact8mer(query_str, target_str)

    score = adjusted_identity + (exact6_bonus if exact6 else 0.0)

    match_mode = "global_epitope_adjusted" if local_enabled else "global_no_epitope"

    return {
        "id": rec.id,
        "description": rec.description,
        "identity": best_identity,
        "identity_adjusted": adjusted_identity,
        "has_6mer": exact6,
        "score": score,
        "match_mode": match_mode,
        "window_contains_epitope": window_contains_epitope,
        "best_window_seq": best_window_seq,
        "window_hit_count": window_hit_count,
        "eightmer_hit_count": eightmer_hit_count
    }



# multi-process optimization
def check_cross_parallel(query_seq, db_records, win=80, id_thresh=0.35, n_jobs=None, local_enabled=FALSE):
    query_str = str(query_seq.seq).upper()
    if n_jobs is None:
        n_jobs = cpu_count()

    args_list = [(query_str, rec, win, id_thresh, local_enabled) for rec in db_records]

    print(f"💻 Using {n_jobs} processes, total candidates: {len(db_records)}")
    results = []
    with Pool(processes=n_jobs) as pool:
        for r in tqdm(pool.imap_unordered(check_one_target, args_list), total=len(args_list), desc="🔍 比对中"):
            if r:
                results.append(r)
    return results


def extract_os(description):
    """从description中提取OS=物种信息"""
    match = re.search(r'OS=([^=]+?)\s*(?:OX=|GN=|PE=|SV=|$)', description)
    return match.group(1).strip() if match else "Unknown"


def dedup_by_os(results):
    os_best = {}
    for r in results:
        os_name = extract_os(r["description"])
        if os_name not in os_best or r["score"] > os_best[os_name]["score"]:
            os_best[os_name] = r
    return list(os_best.values())

def predict_cross_allergen(allergen_name, query_fasta_dir="fasta_files", db_fasta="./naive/uniprotkb_allergen_2025_07_22.fasta", k=5, cpu=None, local_enabled=True):
    db_records = list(SeqIO.parse(db_fasta, "fasta"))

    query_files = [os.path.join(query_fasta_dir, f) for f in os.listdir(query_fasta_dir)
                   if f.lower().endswith(('.fasta', '.fa'))]    # filter, list generator

    print(f"🔎 发现 {len(query_files)} 条待测序列文件")

    all_results = []
    total_start = time.time()

    for idx, qfile in enumerate(query_files, 1):
        query_seq = SeqIO.read(qfile, "fasta")
        print(f"\n[{idx}/{len(query_files)}] 处理序列: {query_seq.id} ({qfile})")

        start = time.time()
        results = check_cross_parallel(query_seq, db_records, win=80, id_thresh=0.35, n_jobs=cpu, local_enabled=local_enabled)
        end = time.time()

        if results:
            results = dedup_by_os(results)
            results.sort(key=lambda x: x["score"], reverse=True)    # key: input=each item in the list, output: value(int, float, str)
            for r in results:
                r['query_id'] = query_seq.id  # 标记是哪个query的结果
            all_results.extend(results)
            print(f"  ✔️ 找到 {len(results)} 个潜在交叉过敏原，耗时 {round(end-start, 2)} s")
        else:
            print("  ⚠️ 未发现潜在交叉过敏原")

        # break

    total_end = time.time()

    if not all_results:
        print("\n❌ 全部序列均未发现潜在交叉过敏原")
        return

    # 汇总结果写文件，按 query_id + score 排序
    all_results.sort(key=lambda x: (x['query_id'], -x['score']))

    with open(f'result_cache/{allergen_name}.txt', 'w', encoding='utf-8') as f:
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

    print(f"\n🎉 所有查询序列处理完成，结果保存在: output_all_queries.txt")
    print(f"总耗时: {round(total_end - total_start, 2)} s")

    content = {}
    keys = all_results[0].keys()
    for key in keys:
        content[key] = []
        
    for item in all_results:
        for key in keys:
            if key == 'id':
                item[key] = item[key].split('|')[1]
            content[key].append(item[key])

    df = pd.DataFrame(content)

    return df


def main():
    parser = argparse.ArgumentParser(description="预测多个潜在交叉过敏原（Top-K）")
    parser.add_argument("--query_fasta_dir", help="待检测蛋白序列所在文件夹，每个文件一个FASTA", default="./fasta_files")
    parser.add_argument("--db_fasta", help="本地数据库（FASTA）", default="./naive/uniprotkb_allergen_2025_07_22.fasta")
    parser.add_argument("-k", type=int, default=5, help="输出最可能交叉反应的 Top-K 个蛋白")
    parser.add_argument("--cpu", type=int, default=None, help="使用的 CPU 核数")
    parser.add_argument("--local_enabled", type=bool, default=False, help="是否启用局部表位匹配")
    args = parser.parse_args()

    db_records = list(SeqIO.parse(args.db_fasta, "fasta"))

    query_files = [os.path.join(args.query_fasta_dir, f) for f in os.listdir(args.query_fasta_dir)
                   if f.lower().endswith(('.fasta', '.fa'))]    # filter, list generator

    print(f"🔎 发现 {len(query_files)} 条待测序列文件")

    all_results = []
    total_start = time.time()

    for idx, qfile in enumerate(query_files, 1):
        query_seq = SeqIO.read(qfile, "fasta")
        print(f"\n[{idx}/{len(query_files)}] 处理序列: {query_seq.id} ({qfile})")

        start = time.time()
        results = check_cross_parallel(query_seq, db_records,  win=80, id_thresh=0.35, n_jobs=args.cpu, local_enabled=args.local_enabled)
        end = time.time()

        if results:
            results = dedup_by_os(results)
            results.sort(key=lambda x: x["score"], reverse=True)    # key: input=each item in the list, output: value(int, float, str)
            for r in results:
                r['query_id'] = query_seq.id  # 标记是哪个query的结果
            import pdb; pdb.set_trace()
            all_results.extend(results)
            print(f"  ✔️ 找到 {len(results)} 个潜在交叉过敏原，耗时 {round(end-start, 2)} s")
        else:
            print("  ⚠️ 未发现潜在交叉过敏原")

        break


    total_end = time.time()

    # if not all_results:
    #     print("\n❌ 全部序列均未发现潜在交叉过敏原")
    #     return

    # 汇总结果写文件，按 query_id + score 排序
    all_results.sort(key=lambda x: (x['query_id'], -x['score']))

    with open('output_all_queries.txt', 'w', encoding='utf-8') as f:
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

    print(f"\n🎉 所有查询序列处理完成，结果保存在: output_all_queries.txt")
    print(f"总耗时: {round(total_end - total_start, 2)} 秒")

    all_results = [{'id': 'tr|D3K177|D3K177_ARAHY', 'description': 'tr|D3K177|D3K177_ARAHY Profilin OS=Arachis hypogaea OX=3818 PE=1 SV=1', 'identity': 0.9625, 'has_6mer': True, 'score': 1.0625, 'match_mode': 'global', 'query_id': 'sp|Q9SQI9|PROF_ARAHY'}]

    content = {}
    keys = all_results[0].keys()
    for key in keys:
        content[key] = []
    for item in all_results:
        for key in keys:
            content[key].append(item[key])

    df = pd.DataFrame(content)



if __name__ == "__main__":
    main()

