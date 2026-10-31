#!/usr/bin/env python3
"""
FluJump — 从 NCBI Protein 下载流感 HA 蛋白序列及元数据
下载 H1/H3/H5/H7 亚型，原样保存不做清洗。
"""

import csv
import os
import re
import time
import calendar
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

# ── 配置 ──────────────────────────────────────────────────
Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY") or None
BATCH_SIZE = 500
SLEEP_INTERVAL = 0.4  # 请求间隔（秒）
OUTPUT_DIR = Path("data/raw")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "H1": "Influenza A virus hemagglutinin H1",
    "H3": "Influenza A virus hemagglutinin H3",
    "H5": "Influenza A virus hemagglutinin H5",
    "H7": "Influenza A virus hemagglutinin H7",
}

# ── 共享字段列表 ──────────────────────────────────────────
FIELD_NAMES = [
    "accession", "strain_name", "subtype", "serotype",
    "host_species", "host_category", "is_human", "organism",
    "collection_date_raw", "collection_year", "collection_month", "collection_day",
    "geo_location", "country",
    "isolation_source", "lab_host", "passage_history", "note",
    "ha_sequence", "seq_length", "segment", "definition",
]


# ══════════════════════════════════════════════════════════
#  helper: 日期解析
# ══════════════════════════════════════════════════════════

# 月份缩写 → 数字
MONTH_MAP = {m.lower(): i for i, m in enumerate(calendar.month_abbr) if i > 0}


def parse_collection_date(raw: str | None):
    """返回 (year, month, day)，每个都是 int 或 None。"""
    if not raw:
        return None, None, None
    raw = raw.strip()

    # 跳过含有 "/" 的年份范围 如 "2013/2014"
    if re.match(r"^\d{4}/\d{4}$", raw):
        return None, None, None

    patterns = [
        # 05-Mar-2013
        (r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$",
         lambda g: (int(g[2]), MONTH_MAP.get(g[1].lower()), int(g[0]))),
        # Mar-2013
        (r"^([A-Za-z]{3})-(\d{4})$",
         lambda g: (int(g[1]), MONTH_MAP.get(g[0].lower()), None)),
        # 2013-03-05
        (r"^(\d{4})-(\d{1,2})-(\d{1,2})$",
         lambda g: (int(g[0]), int(g[1]), int(g[2]))),
        # 2013/03/05
        (r"^(\d{4})/(\d{1,2})/(\d{1,2})$",
         lambda g: (int(g[0]), int(g[1]), int(g[2]))),
        # 2013-03
        (r"^(\d{4})-(\d{1,2})$",
         lambda g: (int(g[0]), int(g[1]), None)),
        # 2013
        (r"^(\d{4})$",
         lambda g: (int(g[0]), None, None)),
    ]

    for pat, extractor in patterns:
        m = re.match(pat, raw)
        if m:
            return extractor(m.groups())
    # 无法解析
    return None, None, None


# ══════════════════════════════════════════════════════════
#  helper: host_category 推断
# ══════════════════════════════════════════════════════════

# 注意：顺序重要 — 更具体的先匹配
# 09-19：匹配改为 \b 词边界（允许复数 s），防止裸子串误伤
# （cat→furcata/flycatcher、pig→pigeon、cow→Moscow 等，见 clean_data.py）
HOST_RULES = [
    (["homo sapiens", "human", "homo", "patient"], "human"),
    (["chicken", "duck", "goose", "quail", "turkey", "mallard",
      "pheasant", "pigeon", "poultry", "wild bird", "teal",
      "gull", "shorebird", "swan", "fowl", "avian",
      "waterbird", "waterfowl", "bluebird", "flycatcher", "guineafowl"], "avian"),
    (["swine", "pig", "porcine", "hog"], "swine"),
    (["equine", "horse"], "equine"),
    (["canine", "dog"], "canine"),
    (["bovine", "cattle", "cow", "dairy"], "bovine"),
    (["feline", "cat"], "feline"),
    (["ferret", "mink"], "mustelid"),
    (["seal", "whale", "marine"], "marine_mammal"),
    (["mouse", "bat", "rodent"], "other_mammal"),
    (["environment", "water", "air", "sewage"], "environmental"),
]

_HOST_RULES_C = [([re.compile(r'\b' + re.escape(kw) + r's?\b') for kw in kws], cat)
                 for kws, cat in HOST_RULES]


def infer_host_category(host_species: str | None) -> str:
    """从 host_species 字符串推断 host_category。"""
    if not host_species:
        return "unknown"
    hs = host_species.lower()
    for patterns, category in _HOST_RULES_C:
        for pat in patterns:
            if pat.search(hs):
                return category
    return "unknown"


# ══════════════════════════════════════════════════════════
#  helper: 从 qualifiers 安全获取
# ══════════════════════════════════════════════════════════

def _get(qualifiers: dict, key: str) -> str:
    """返回 qualifier 的第一个值，不存在返回空字符串。"""
    vals = qualifiers.get(key, [])
    return vals[0] if vals else ""


# ══════════════════════════════════════════════════════════
#  helper: 提取 country
# ══════════════════════════════════════════════════════════

def extract_country(geo_location: str) -> str:
    """从 geo_location 提取 country（取冒号前部分）。"""
    if not geo_location:
        return ""
    return geo_location.split(":")[0].strip()


# ══════════════════════════════════════════════════════════
#  helper: 提取 strain_name
# ══════════════════════════════════════════════════════════

def extract_strain_name(qualifiers: dict, definition: str) -> str:
    """/strain > /isolate > definition 首段。"""
    strain = _get(qualifiers, "strain")
    if strain:
        return strain
    isolate = _get(qualifiers, "isolate")
    if isolate:
        return isolate
    # definition 首段（到第一个逗号或分号）
    if definition:
        m = re.match(r"^([^,;]+)", definition)
        if m:
            return m.group(1).strip()
    return ""


# ══════════════════════════════════════════════════════════
#  解析单条 GenBank 记录
# ══════════════════════════════════════════════════════════

def parse_record(record, subtype: str) -> dict:
    """从 Bio.SeqRecord 提取所有字段。"""
    feats = record.features
    quals = feats[0].qualifiers if feats else {}

    # 序列
    ha_seq = str(record.seq).upper()

    # qualifier 字段
    host_species = _get(quals, "host") or _get(quals, "organism")
    organism = _get(quals, "organism")
    collection_date_raw = _get(quals, "collection_date")
    geo_location = _get(quals, "country") or _get(quals, "geo_loc_name")
    isolation_source = _get(quals, "isolation_source")
    lab_host = _get(quals, "lab_host")
    passage_history = _get(quals, "note") or _get(quals, "passage_history")
    note = _get(quals, "note")
    segment = _get(quals, "segment")
    serotype = _get(quals, "serotype")

    # 日期解析
    year, month, day = parse_collection_date(collection_date_raw)

    # host_category
    host_category = infer_host_category(host_species)

    # country
    country = extract_country(geo_location)

    # strain_name
    strain_name = extract_strain_name(quals, record.description)

    return {
        "accession": record.annotations.get("accessions", [""])[0] if record.annotations.get("accessions") else record.id,
        "strain_name": strain_name,
        "subtype": subtype,
        "serotype": serotype,
        "host_species": host_species,
        "host_category": host_category,
        "is_human": 1 if host_category == "human" else 0,
        "organism": organism,
        "collection_date_raw": collection_date_raw,
        "collection_year": year if year else "",
        "collection_month": month if month else "",
        "collection_day": day if day else "",
        "geo_location": geo_location,
        "country": country,
        "isolation_source": isolation_source,
        "lab_host": lab_host,
        "passage_history": passage_history,
        "note": note,
        "ha_sequence": ha_seq,
        "seq_length": len(ha_seq),
        "segment": segment,
        "definition": record.description,
    }


# ══════════════════════════════════════════════════════════
#  下载 & 解析
# ══════════════════════════════════════════════════════════

def download_subtype(subtype: str, query: str):
    """下载一个亚型的所有记录，写入 CSV。"""
    output_path = OUTPUT_DIR / f"{subtype}_raw.csv"
    print(f"\n{'='*60}")
    print(f"  开始下载 {subtype}: {query}")
    print(f"{'='*60}")

    # Step 1 — ESearch（获取 WebEnv / QueryKey）
    print("  [1/3] ESearch 查询中...")
    handle = Entrez.esearch(
        db="protein",
        term=query,
        retmax=0,
        usehistory="y",
    )
    result = Entrez.read(handle)
    handle.close()
    count = int(result["Count"])
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    print(f"  → 共 {count} 条记录")
    if count == 0:
        print(f"  ⚠ 无结果，跳过。")
        return

    # Step 2 — 分批 EFetch + 解析
    print(f"  [2/3] 开始分批下载（batch_size={BATCH_SIZE}）...")

    all_rows = []
    for start in range(0, count, BATCH_SIZE):
        end = min(start + BATCH_SIZE, count)
        retries = 3
        data = None
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(
                    db="protein",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=BATCH_SIZE,
                    webenv=webenv,
                    query_key=query_key,
                )
                data = handle.read()
                handle.close()
                break
            except Exception as e:
                print(f"    ⚠ 第 {attempt+1} 次尝试失败: {e}")
                time.sleep(2)
        if data is None:
            print(f"    ✗ 批次 {start}-{end} 下载失败，跳过。")
            continue

        # 解析 GenBank 文本
        try:
            records = list(SeqIO.parse(StringIO(data), "genbank"))
        except Exception as e:
            print(f"    ✗ 解析失败: {e}")
            continue

        for rec in records:
            row = parse_record(rec, subtype)
            all_rows.append(row)

        print(f"    ✓ {start+1}-{end}/{count}  ({len(all_rows)} 已解析)")

        if end < count:
            time.sleep(SLEEP_INTERVAL)

    print(f"  → 共解析 {len(all_rows)} 条")

    # Step 3 — 写 CSV
    print(f"  [3/3] 写入 {output_path} ...")
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELD_NAMES)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"  ✓ 完成: {output_path}")


# ══════════════════════════════════════════════════════════
#  main
# ══════════════════════════════════════════════════════════

def main():
    if not Entrez.email:
        raise SystemExit("请先通过 NCBI_EMAIL 环境变量设置联系邮箱")
    total_start = time.time()
    for subtype, query in QUERIES.items():
        subtype_start = time.time()
        download_subtype(subtype, query)
        elapsed = time.time() - subtype_start
        print(f"  ⏱  {subtype} 耗时: {elapsed:.1f} 秒")

    total_elapsed = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"  全部完成！总耗时: {total_elapsed:.1f} 秒 ({total_elapsed/60:.1f} 分钟)")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
