import streamlit as st
import pandas as pd
import json
import time
import re
import os
from collections import defaultdict
import random

import requests
from bs4 import BeautifulSoup

from parser import get_uniprot_id_from_detail_page

from io import StringIO

from st_func import *

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
            "has_8mer_match": None,
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
                    record["has_8mer_match"] = match.group(1) == "True"
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
print(f"input species: {input_species}")
top_n = st.slider("🔢 显示 Top-N 结果", 5, 20, 10)

predict_button = st.button("🚀 开始预测")

# ========== 主流程 ==========
if predict_button:
    # 1. 搜索并保存该物种的过敏原
    os.makedirs("species_cache", exist_ok=True)
    allergen_csv_path = f"species_cache/{input_species}_allergens.csv"

    # 1. 优先从 cache 读取物种过敏原
    if os.path.exists(allergen_csv_path):
        allergen_df = pd.read_csv(allergen_csv_path)
    else:
        st.info(f"🔍 正在搜索 {input_species} 的过敏原...")
        allergens = search_allergen_org_with_real_uniprot_st(input_species)
        allergen_df = pd.DataFrame(allergens)
        allergen_df.to_csv(allergen_csv_path, index=False, encoding="utf-8")

    # 展示过敏原列表
    st.subheader(f"📋 {input_species} 的过敏原列表")
    st.dataframe(allergen_df,height=200)

    # 下载按钮
    csv_data = allergen_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 下载该物种过敏原列表（CSV）",
        data=csv_data,
        file_name=f"{input_species}_allergens.csv",
        mime="text/csv"
    )

    # 2. 交叉过敏原预测
    if os.path.exists(f"result_cache/{input_species}.txt"):
        df = parse_allergen_text(f"result_cache/{input_species}.txt")
    else:
        with st.spinner("🧠 正在预测可能的交叉过敏原物种，请稍候..."):
            df = predict_cross_allergen_streamlit(species_name=input_species)
            df = df.sort_values(by="score", ascending=False).reset_index(drop=True)

    # 展示过敏原列表
    st.subheader(f"📋 {input_species} 的交叉过敏原预测结果")
    st.dataframe(df,height=200)

    # 下载按钮
    csv_data = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label="📥 下载结果 CSV",
        data=csv_data,
        file_name=f"{input_species}_cross_allergen_results.csv",
        mime="text/csv"
    )

