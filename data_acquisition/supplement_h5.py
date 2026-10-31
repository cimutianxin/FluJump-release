#!/usr/bin/env python3
"""
H5 补充下载：牛源 (Bos taurus) + 人源 (Homo sapiens) HA 序列
与现有 H5_raw.csv 合并去重后，覆盖原文件。
"""

import csv
import os
import re
import time
import calendar
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

Entrez.email = os.environ.get("NCBI_EMAIL", "")
Entrez.api_key = os.environ.get("NCBI_API_KEY") or None
BATCH_SIZE = 500
SLEEP = 0.4

QUERIES = {
    "cattle": 'H5N1 hemagglutinin "Bos taurus"',
    "human": 'H5N1 hemagglutinin "Homo sapiens"',
}

# ── 字段（与 download_flu.py 一致）──
FIELDS = [
    "accession", "strain_name", "subtype", "serotype",
    "host_species", "host_category", "is_human", "organism",
    "collection_date_raw", "collection_year", "collection_month", "collection_day",
    "geo_location", "country",
    "isolation_source", "lab_host", "passage_history", "note",
    "ha_sequence", "seq_length", "segment", "definition",
]

# ── 复用的 helper ──
MONTH_MAP = {m.lower(): i for i, m in enumerate(calendar.month_abbr) if i > 0}

HOST_RULES = [
    (["homo sapiens", "human", "homo", "patient"], "human"),
    (["chicken", "duck", "goose", "quail", "turkey", "mallard",
      "pheasant", "pigeon", "poultry", "wild bird", "teal",
      "gull", "shorebird", "swan", "fowl",
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

GENUS_RULES = [
    (["anas", "gallus", "meleagris", "coturnix", "larus", "calidris",
      "tyto", "cairina"], "avian"),
    (["sus"], "swine"),
    (["homo"], "human"),
    (["canis"], "canine"),
    (["equus"], "equine"),
    (["bos"], "bovine"),
    (["felis"], "feline"),
    (["mustela", "neovison"], "mustelid"),
    (["phoca", "tursiops", "balaenoptera"], "marine_mammal"),
    (["mus", "rattus", "myotis", "pteropus"], "other_mammal"),
]

_HOST_RULES_C = [([re.compile(r'\b' + re.escape(kw) + r's?\b') for kw in kws], cat)
                 for kws, cat in HOST_RULES]


def _get(q, key):
    vals = q.get(key, [])
    return vals[0] if vals else ""


def infer_host(host_species, strain_name):
    if not host_species:
        host_species = ""
    hs = host_species.lower()
    if "synthetic" in hs:
        return "unknown"
    for pats, cat in _HOST_RULES_C:
        for p in pats:
            if p.search(hs):
                return cat
    for genera, cat in GENUS_RULES:
        for g in genera:
            if re.search(r'\b' + re.escape(g) + r'\b', hs):
                return cat
    if strain_name:
        sn = strain_name.lower()
        patterns = [
            (r"/swine/|/pig/", "swine"),
            (r"/duck/|/chicken/|/goose/|/quail/|/turkey/|/mallard/|/pheasant/"
             r"|/pigeon/|/teal/|/gull/|/avian/", "avian"),
            (r"/cattle/|/bovine/|/dairy/|/cow/", "bovine"),
            (r"/human/|/homo /|/patient/", "human"),
            (r"/equine/|/horse/", "equine"),
            (r"/canine/|/dog/", "canine"),
            (r"/feline/|/cats?/", "feline"),
        ]
        for pat, cat in patterns:
            if re.search(pat, sn):
                return cat
    return "unknown"


def parse_date(raw):
    if not raw:
        return "", "", ""
    raw = raw.strip()
    if re.match(r"^\d{4}/\d{4}$", raw):
        return "", "", ""
    tests = [
        (r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})$",
         lambda g: (g[2], str(MONTH_MAP.get(g[1].lower(), "")), g[0])),
        (r"^([A-Za-z]{3})-(\d{4})$",
         lambda g: (g[1], str(MONTH_MAP.get(g[0].lower(), "")), "")),
        (r"^(\d{4})-(\d{1,2})-(\d{1,2})$",
         lambda g: (g[0], g[1], g[2])),
        (r"^(\d{4})/(\d{1,2})/(\d{1,2})$",
         lambda g: (g[0], g[1], g[2])),
        (r"^(\d{4})-(\d{1,2})$",
         lambda g: (g[0], g[1], "")),
        (r"^(\d{4})$",
         lambda g: (g[0], "", "")),
    ]
    for pat, fn in tests:
        m = re.match(pat, raw)
        if m:
            return fn(m.groups())
    return "", "", ""


def extract_strain(quals, defn):
    s = _get(quals, "strain")
    if s: return s
    iso = _get(quals, "isolate")
    if iso: return iso
    if defn:
        m = re.match(r"^([^,;]+)", defn)
        if m: return m.group(1).strip()
    return ""


def parse_record(rec, subtype):
    feats = rec.features
    quals = feats[0].qualifiers if feats else {}
    seq = str(rec.seq).upper()
    hs = _get(quals, "host") or _get(quals, "organism")
    org = _get(quals, "organism")
    date_raw = _get(quals, "collection_date")
    geo = _get(quals, "country") or _get(quals, "geo_loc_name")
    strain = extract_strain(quals, rec.description)
    hc = infer_host(hs, strain)
    y, m, d = parse_date(date_raw)
    country = geo.split(":")[0].strip() if geo else ""
    acc = rec.annotations.get("accessions", [rec.id])
    if isinstance(acc, list): acc = acc[0] if acc else rec.id

    return {
        "accession": acc,
        "strain_name": strain,
        "subtype": subtype,
        "serotype": _get(quals, "serotype"),
        "host_species": hs,
        "host_category": hc,
        "is_human": 1 if hc == "human" else 0,
        "organism": org,
        "collection_date_raw": date_raw,
        "collection_year": y,
        "collection_month": m,
        "collection_day": d,
        "geo_location": geo,
        "country": country,
        "isolation_source": _get(quals, "isolation_source"),
        "lab_host": _get(quals, "lab_host"),
        "passage_history": _get(quals, "note") or _get(quals, "passage_history"),
        "note": _get(quals, "note"),
        "ha_sequence": seq,
        "seq_length": len(seq),
        "segment": _get(quals, "segment"),
        "definition": rec.description,
    }


def download_query(label, query, subtype):
    print(f"\n  [{label}] {query}")
    handle = Entrez.esearch(db="protein", term=query, retmax=0, usehistory="y")
    r = Entrez.read(handle); handle.close()
    count = int(r["Count"])
    webenv, qkey = r["WebEnv"], r["QueryKey"]
    print(f"  → {count} 条")

    rows = []
    for start in range(0, count, BATCH_SIZE):
        end = min(start + BATCH_SIZE, count)
        data = None
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="protein", rettype="gb", retmode="text",
                                  retstart=start, retmax=BATCH_SIZE,
                                  webenv=webenv, query_key=qkey)
                data = h.read(); h.close()
                break
            except Exception as e:
                print(f"    ⚠ retry {attempt+1}: {e}")
                time.sleep(2)
        if not data: continue
        for rec in SeqIO.parse(StringIO(data), "genbank"):
            rows.append(parse_record(rec, subtype))
        print(f"    {start+1}-{end}/{count} ({len(rows)} parsed)")
        if end < count: time.sleep(SLEEP)
    return rows


def main():
    if not Entrez.email:
        raise SystemExit("请先通过 NCBI_EMAIL 环境变量设置联系邮箱")
    print("=" * 60)
    print("  H5 补充下载: 牛源 + 人源")
    print("=" * 60)

    # 加载现有 H5 数据
    h5_path = Path("data/raw/H5_raw.csv")
    existing = {}
    if h5_path.exists():
        with open(h5_path, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                existing[row["accession"]] = row
    print(f"\n  现有 H5: {len(existing)} 条")

    # 补充下载
    new_total = 0
    for label, query in QUERIES.items():
        rows = download_query(label, query, "H5")
        new = 0
        for row in rows:
            acc = row["accession"]
            if acc not in existing:
                existing[acc] = row
                new += 1
            # 如果已存在但 host 信息缺失，用新数据补充
            elif existing[acc].get("host_category") in ("unknown", "") and \
                 row.get("host_category") not in ("unknown", ""):
                existing[acc] = row
                new += 1
        new_total += new
        print(f"  → 新增 {new} 条（去重后）")

    print(f"\n  最终 H5: {len(existing)} 条（+{new_total}）")

    # 写回
    backup_path = Path("data/raw/H5_raw_backup.csv")
    if h5_path.exists():
        h5_path.rename(backup_path)
        print(f"  备份: {backup_path}")

    with open(h5_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        # 保持原有顺序 + 新条目
        all_accs = list(existing.keys())
        writer.writerows(existing[acc] for acc in all_accs)

    print(f"  写入: {h5_path}")
    print(f"\n{'='*60}")
    print(f"  完成！H5 数据已更新")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
