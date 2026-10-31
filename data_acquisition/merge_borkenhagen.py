#!/usr/bin/env python3
"""
解析 Borkenhagen 2024 补充数据，合并 Data S1 (FASTA) + Data S2 (标签)。
三轮匹配: 精确归一化 → token 重叠 → 相似度兜底。
"""

import pandas as pd, re
from difflib import SequenceMatcher
from pathlib import Path

DATA_DIR = Path("data/dataset_borkenhagen_raw")
LABEL_CSV = DATA_DIR / "irv70044-sup-0003-supplementary_file2.csv"
FASTA = DATA_DIR / "DataS1_sequences.fasta"
OUT = Path("data/dataset_borkenhagen/borkenhagen_merged.csv")
OUT.parent.mkdir(parents=True, exist_ok=True)

def norm(s):
    return re.sub(r'[^a-z0-9]', '', s.lower())

def tokens(s):
    return set(re.findall(r'[a-z0-9]{4,}', norm(s)))

def parse_fasta(path):
    seqs = {}
    cid, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if cid: seqs[cid] = ''.join(buf)
                cid = line[1:].split('|')[0].strip(); buf = []
            elif cid:
                buf.append(line)
    if cid: seqs[cid] = ''.join(buf)
    return seqs

# ── Load ──
labels = pd.read_csv(LABEL_CSV, encoding='latin-1')
seqs = parse_fasta(FASTA)
print(f"Labels: {len(labels)}, Sequences: {len(seqs)}")

# ── Match ──
match = {}
used = set()

# R1: exact normalized
for lid in labels['ID']:
    nl = norm(lid)
    for fk in seqs:
        if norm(fk) == nl and fk not in used:
            match[lid] = fk; used.add(fk); break

# R2: token overlap (>=3 shared tokens, similarity > 0.7)
for lid in labels['ID']:
    if lid in match: continue
    lt = tokens(lid)
    best_k, best_n = None, 0
    for fk in seqs:
        if fk in used: continue
        n = len(lt & tokens(fk))
        if n >= 3 and n > best_n and SequenceMatcher(None, norm(lid), norm(fk)).ratio() > 0.7:
            best_n, best_k = n, fk
    if best_k:
        match[lid] = best_k; used.add(best_k)

# R3: similarity fallback
for lid in labels['ID']:
    if lid in match: continue
    nl = norm(lid)
    best_k, best_s = None, 0
    for fk in seqs:
        if fk in used: continue
        s = SequenceMatcher(None, nl, norm(fk)).ratio()
        if s > best_s: best_s, best_k = s, fk
    if best_s > 0.6 and best_k:
        match[lid] = best_k; used.add(best_k)

# ── Build ──
merged = labels.copy()
merged['fasta_key'] = merged['ID'].map(match)
merged['ha_sequence'] = merged['fasta_key'].map(seqs)
merged = merged.dropna(subset=['ha_sequence'])
merged['seq_len'] = merged['ha_sequence'].str.len()

# ── Report ──
pos = merged['Label'].sum()
print(f"Merged: {len(merged)} 条  |  a2,6={pos}  |  a2,3={len(merged)-pos}")
print(f"Seq len: {merged['seq_len'].value_counts().to_dict()}")
unmatched = [k for k in labels['ID'] if k not in match]
if unmatched:
    print(f"Unmatched ({len(unmatched)}): {unmatched}")
merged.to_csv(OUT, index=False)
print(f"✓ {OUT}")
