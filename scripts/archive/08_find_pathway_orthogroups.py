#!/usr/bin/env python3
"""
08_find_pathway_orthogroups.py

Locate the carotenoid and MEP pathway orthogroups, and count how many taro
genes fall into each.

The orthogroup table identifies sequences by protein accession, because that
is what the FASTA headers carry. Gene symbols therefore have to be resolved to
protein accessions first, which is done by querying NCBI rather than by
hard-coding accession numbers, so that nothing here depends on my memory.

Arabidopsis is used as the anchor species because its carotenoid pathway is
the best characterised and its gene symbols are stable.

Usage:
    python scripts/08_find_pathway_orthogroups.py
"""
import json, os, subprocess, sys
from collections import defaultdict

REFS = os.environ.get("TARO_REFS")
WORK = os.environ.get("TARO_WORK")
if not REFS or not WORK:
    sys.exit("source envs/activate.sh first")

# Arabidopsis gene symbols for the pathway. Resolved to protein accessions
# below; nothing is assumed about the accession numbers themselves.
PATHWAY = {
    "MEP backbone":   ["DXS", "DXR", "ISPD", "ISPE", "ISPF", "ISPG", "ISPH"],
    "carotenoid core":["PSY", "PDS3", "ZDS", "CRTISO", "Z-ISO"],
    "cyclisation":    ["LUT2", "LCYB", "LYC"],
    "xanthophyll":    ["LUT1", "LUT5", "CHY1", "CHY2", "ABA1", "NPQ1", "ABA4"],
    "degradation":    ["CCD1", "CCD4", "NCED3"],
    "sink / stability":["OR", "ORL", "PGL35"],
}

def find(pattern, root):
    for dirpath, _, files in os.walk(root):
        for f in files:
            if f == pattern:
                return os.path.join(dirpath, f)
    return None

OG = find("Orthogroups.tsv", os.path.join(WORK, "05_orthology"))
if not OG:
    sys.exit("Orthogroups.tsv not found; has OrthoFinder finished writing?")
print(f"orthogroups: {OG}\n")

# ---- resolve symbols to Arabidopsis protein accessions via NCBI -----------
def proteins_for(symbol):
    try:
        out = subprocess.run(
            ["datasets", "summary", "gene", "symbol", symbol,
             "--taxon", "Arabidopsis thaliana"],
            capture_output=True, text=True, timeout=60).stdout
        d = json.loads(out)
    except Exception:
        return []
    accs = []
    for r in d.get("reports", []):
        g = r.get("gene", {})
        if g.get("symbol", "").upper() != symbol.upper():
            continue
        for t in g.get("transcripts", []) or []:
            p = (t.get("protein") or {}).get("accession_version")
            if p:
                accs.append(p)
    return sorted(set(accs))

# ---- read orthogroups ----------------------------------------------------
header = None
acc2og = {}
og_rows = {}
with open(OG) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        og = parts[0]
        og_rows[og] = dict(zip(header[1:], parts[1:]))
        for cell in parts[1:]:
            for item in cell.split(", "):
                item = item.strip()
                if item:
                    acc2og[item.split("|")[-1]] = og

species = header[1:]
print(f"species in table: {', '.join(species)}\n")

# ---- walk the pathway ----------------------------------------------------
rows = []
for group, symbols in PATHWAY.items():
    print(f"\n=== {group} ===")
    for sym in symbols:
        accs = proteins_for(sym)
        if not accs:
            print(f"  {sym:<10} symbol not resolved at NCBI")
            continue
        ogs = {acc2og[a] for a in accs if a in acc2og}
        if not ogs:
            print(f"  {sym:<10} {len(accs)} protein(s), none in any orthogroup")
            continue
        for og in sorted(ogs):
            r = og_rows[og]
            counts = {s: (len([x for x in r[s].split(", ") if x.strip()])
                          if r.get(s) else 0) for s in species}
            taro = counts.get("colesc", 0)
            arab = counts.get("aratha", 0)
            flag = "" if arab == 1 else f"  <-- CHECK: {arab} Arabidopsis genes"
            print(f"  {sym:<10} {og}  taro={taro:<3} arabidopsis={arab}{flag}")
            rows.append((group, sym, og, taro, arab, counts))

# ---- the headline --------------------------------------------------------
print("\n" + "="*66)
psy = [r for r in rows if r[1] == "PSY"]
if psy:
    for _, _, og, taro, arab, counts in psy:
        print(f"PSY orthogroup {og}")
        print(f"  taro (Colocasia esculenta) : {taro} gene(s)")
        print(f"  Arabidopsis                : {arab} gene(s)"
              + ("   [expected 1 — sanity check passed]" if arab == 1
                 else "   [expected 1 — SOMETHING IS WRONG UPSTREAM]"))
        print("  full composition:")
        for s in species:
            print(f"    {s:<10} {counts[s]}")
else:
    print("PSY not located. Try resolving the accession manually.")

# ---- save ----------------------------------------------------------------
import csv
outdir = os.path.join(os.path.dirname(__file__), "..", "results", "tables")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "pathway_orthogroups.tsv")
with open(outfile, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pathway_step", "symbol", "orthogroup"] + species)
    for group, sym, og, taro, arab, counts in rows:
        w.writerow([group, sym, og] + [counts[s] for s in species])
print(f"\nwritten: {os.path.normpath(outfile)}")
