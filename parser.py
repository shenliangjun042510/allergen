import os
import requests
from bs4 import BeautifulSoup

def get_uniprot_id_from_detail_page(iuis_id: str) -> str | None:
    url = f"https://www.allergen.org/viewallergen.php?aid={iuis_id}"
    resp = requests.get(url)
    if resp.status_code != 200:
        print(f"❌ 详情页请求失败: {iuis_id}")
        return None

    soup = BeautifulSoup(resp.text, "html.parser")
    table = soup.find("table", {"id": "isotable"})
    if not table:
        print(f"⚠️ 没有找到 id='isotable' 的表格")
        return None

    for a in table.find_all("a", href=True):
        href = a["href"]
        if "uniprot.org" in href:
            return href.split("/")[-1].strip()
    return None


def search_allergen_org_with_real_uniprot(source_keyword="peanut", fasta_dir="./fasta_files"):
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

    resp = requests.get(search_url, params=params)
    # resp.text: html 网页源码

    if resp.status_code != 200:
        print("❌ 请求失败:", resp.status_code)
        return []

    soup = BeautifulSoup(resp.text, "html.parser")
    table = soup.find("table", {"border": "1"})
    if not table:
        print("❌ 未找到结果表格")
        return []

    rows = table.find_all("tr")[2:]

    total_allergens = len(rows)
    print(f"共找到 {total_allergens} 个过敏原条目")

    # 创建保存fasta的文件夹
    os.makedirs(fasta_dir, exist_ok=True)

    results = []

    for idx, row in enumerate(rows, 1):
        cols = row.find_all("td")

        allergen_name = cols[1].text.strip()
        biological_name = cols[2].text.strip()
        iuis_id = cols[1].find("a")["href"].split("=")[-1].strip()

        uniprot_id = get_uniprot_id_from_detail_page(iuis_id)

        sequence = None
        if uniprot_id:
            fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
            fasta_resp = requests.get(fasta_url)
            if fasta_resp.ok:
                sequence = fasta_resp.text

                # 打印当前进度
                print(f"[{idx}/{total_allergens}] 已获取 UniProt 序列: {uniprot_id} ({allergen_name})")

                # 保存到文件
                fasta_path = os.path.join(fasta_dir, f"{uniprot_id}.fasta")
                with open(fasta_path, "w", encoding="utf-8") as f:
                    f.write(sequence)
            else:
                print(f"⚠️ UniProt请求失败: {uniprot_id}")
        else:
            print(f"⚠️ 未找到 UniProt ID: {allergen_name}")

        results.append({
            "allergen_name": allergen_name,
            "biological_name": biological_name,
            "iuis_id": iuis_id,
            "uniprot_id": uniprot_id,
            "sequence": sequence
        })

    print(f"\n所有可用序列已保存到目录: {os.path.abspath(fasta_dir)}")
    return results


if __name__ == "__main__":
    results = search_allergen_org_with_real_uniprot("peanut")
