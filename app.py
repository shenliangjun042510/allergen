import streamlit as st
import pandas as pd
import json
import time
import re
import os
from collections import defaultdict
import random
from io import StringIO

from parser import search_allergen_org_with_real_uniprot
from cross_allergen import search_cross_allergen


# ========== 模拟预测函数 ==========
def simulate_prediction(input_species: str, delay_sec=3) -> pd.DataFrame:
    """模拟一个预测函数，延迟几秒后返回伪造数据"""
    time.sleep(delay_sec)  # 模拟长时间计算

    mock_species = [
        "Pistacia vera", "Corylus avellana", "Juglans regia",
        "Glycine max", "Oryza sativa", "Triticum aestivum",
        "Solanum lycopersicum", "Malus domestica", "Quercus suber"
    ]
    random.shuffle(mock_species)

    n = random.randint(5, 9)
    df = pd.DataFrame({
        "species": mock_species[:n],
        "max_score": [round(random.uniform(0.35, 0.95), 3) for _ in range(n)],
        "avg_score": [round(random.uniform(0.3, 0.9), 3) for _ in range(n)],
        "hits": [random.randint(1, 5) for _ in range(n)]
    })

    # 添加属名和是否为本属
    input_genus = input_species.strip().split()[0]
    df["genus"] = df["species"].apply(lambda x: x.split()[0])
    df["related_to_input"] = df["genus"].str.lower() == input_genus.lower()
    return df


@st.cache_data
def load_results(file_path="output.txt"):
    """解析 output.txt 中的预测结果"""
    results = []
    with open(file_path, encoding='utf-8') as f:
        lines = f.readlines()

    current_desc = None
    current_score = None

    for line in lines:
        if line.startswith("   📛 名称:"):
            current_desc = line.strip().replace("📛 名称:", "").strip()
        elif line.startswith("   📊 综合得分:"):
            score_str = line.strip().replace("📊 综合得分:", "").strip()
            try:
                current_score = float(score_str)
            except ValueError:
                current_score = None

            if current_desc and current_score is not None:
                results.append({
                    "description": current_desc,
                    "score": current_score
                })
                current_desc = None
                current_score = None
    return results


def extract_os(description):
    match = re.search(r'OS=([^=]+?)\s*(?:OX=|GN=|PE=|SV=|$)', description)
    return match.group(1).strip() if match else "Unknown"

def get_genus(species_name):
    return species_name.strip().split()[0]

def parse_allergen_text(path: str) -> pd.DataFrame:
    with open(path, 'r') as f:
        text = f.read()
        
    entries = re.split(r"\n\d+\.\s+", text.strip())
    records = []

    for entry in entries:
        lines = entry.strip().split("\n")
        if not lines:
            continue
        
        record = {
            "id": None,
            "name": None,
            "identity": None,
            "has_6mer_match": None,
            "score": None,
            "mode": None
        }

        # 第一行为 ID
        record["id"] = lines[0].strip().split("|")[1]

        for line in lines[1:]:
            line = line.strip()
            if line.startswith("📛"):
                record["name"] = extract_os(line.split("📛 名称:")[1].strip())
            elif line.startswith("🧬"):
                match = re.search(r"Identity:\s*([\d.]+)", line)
                if match:
                    record["identity"] = float(match.group(1))
            elif line.startswith("🔗"):
                match = re.search(r"Has 6-mer match:\s*(True|False)", line)
                if match:
                    record["has_6mer_match"] = match.group(1) == "True"
            elif line.startswith("📊"):
                match = re.search(r"综合得分:\s*([\d.]+)", line)
                if match:
                    record["score"] = float(match.group(1))
            elif line.startswith("🔍"):
                match = re.search(r"匹配模式:\s*(\w+)", line)
                if match:
                    record["mode"] = match.group(1)
        if None in record.values():
            continue
        records.append(record)

    return pd.DataFrame(records)


def summarize_species_dedup_genus(results):
    species_scores = defaultdict(list)

    for r in results:
        os_name = extract_os(r["description"])
        species_scores[os_name].append(r["score"])

    summary = []
    for species, scores in species_scores.items():
        summary.append({
            "species": species,
            "genus": get_genus(species),
            "max_score": max(scores),
            "avg_score": sum(scores) / len(scores),
            "hits": len(scores)
        })

    best_per_genus = {}
    for item in summary:
        genus = item["genus"]
        if genus not in best_per_genus or item["max_score"] > best_per_genus[genus]["max_score"]:
            best_per_genus[genus] = item

    return sorted(best_per_genus.values(), key=lambda x: x["max_score"], reverse=True)

# ========== Streamlit 页面 ==========
st.set_page_config(page_title="交叉过敏原预测器", layout="wide")
st.title("🧬 蛋白质交叉过敏原预测")

# html, css, javascripts
st.markdown(
    """
    <div style="text-align: center;">
        <img src="https://www.cusabio.com/manage/upload/202109/antibody-cross-reaction.png" width="600"/>
        <p>cross allergen</p>
    </div>
    """,
    unsafe_allow_html=True
)

input_species = st.text_input("🔍 输入物种名称（例如 peanut, shrimp, soy）：", value="peanut")
top_n = st.slider("🔢 显示 Top-N 结果", 5, 20, 10)

predict_button = st.button("🚀 开始预测")

if predict_button:
    if os.path.exists(f"result_cache/{input_species}.txt"):
        df = parse_allergen_text(f"result_cache/{input_species}.txt")
    else:
        results = search_allergen_org_with_real_uniprot(input_species)
        # context handler
        with st.spinner("🧠 正在预测可能的交叉过敏原物种，请稍候..."):
            df = search_cross_allergen(allergen_name="input_species")
            df = df.sort_values(by="score", ascending=False).reset_index(drop=True)

            # df = simulate_prediction(input_species, delay_sec=3)
            # df = df.sort_values(by="max_score", ascending=False).reset_index(drop=True)

    st.success("✅ 预测完成！以下是结果：")

    st.subheader("📋 Top-N 预测结果")
    top_df = df.head(top_n)
    # st.dataframe(top_df.style.applymap(
    #     lambda v: "background-color: #ffd8d8" if isinstance(v, bool) and v else "",
    #     subset=["related_to_input"]
    # ))
    st.dataframe(top_df)

    # 下载按钮
    csv_data = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 下载全部预测结果（CSV）",
        data=csv_data,
        file_name=f"{input_species}_prediction.csv",
        mime="text/csv"
    )
