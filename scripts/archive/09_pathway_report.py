#!/usr/bin/env python3
"""Report carotenoid pathway orthogroups and taro copy number."""
import os, sys, csv
from collections import defaultdict

REFS, WORK = os.environ.get("TARO_REFS"), os.environ.get("TARO_WORK")
CODE = os.environ.get("TARO_CODE")
if not all([REFS, WORK, CODE]):
    sys.exit("source envs/activate.sh first")

# Arabidopsis locus identifiers. Stable, unambiguous, and independent of
# whatever symbol NCBI happens to prefer.
PATHWAY = [
    ("MEP backbone",    "DXS",     "AT4G15560"),
    ("MEP backbone",    "DXR",     "AT5G62790"),
    ("MEP backbone",    "MCT",     "AT2G02500"),
    ("MEP backbone",    "CMK",     "AT2G26930"),
    ("MEP backbone",    "MDS",     "AT1G63970"),
    ("MEP backbone",    "HDS",     "AT5G60600"),
    ("MEP backbone",    "HDR",     "AT4G34350"),
    ("backbone",        "GGPPS11", "AT4G36810"),
    ("carotenoid core", "PSY",     "AT5G17230"),
    ("carotenoid core", "PDS3",    "AT4G14210"),
    ("carotenoid core", "Z-ISO",   "AT1G10830"),
    ("carotenoid core", "ZDS",     "AT3G04870"),
    ("carotenoid core", "CRTISO",  "AT1G06820"),
    ("cyclisation",     "LCYB",    "AT3G10230"),
    ("cyclisation",     "LCYE",    "AT5G57030"),
    ("xanthophyll",     "CYP97A3", "AT1G31800"),
    ("xanthophyll",     "CYP97C1", "AT3G53130"),
    ("xanthophyll",     "BCH1",    "AT4G25700"),
    ("xanthophyll",     "BCH2",    "AT5G52570"),
    ("xanthophyll",     "ZEP",     "AT5G67030"),
    ("xanthophyll",     "VDE",     "AT1G08550"),
    ("xanthophyll",     "NSX",     "AT1G67080"),
    ("degradation",     "CCD1",    "AT3G63520"),
    ("degradation",     "CCD4",    "AT4G19170"),
    ("degradation",     "NCED3",   "AT3G14440"),
    ("sink/stability",  "OR",      "AT5G61670"),
    ("sink/stability",  "OR-like", "AT5G06130"),
    ("sink/stability",  "FBN1a",   "AT4G04020"),
]

def find(name, root):
    for dp, _, fs in os.walk(root):
        if name in fs:
            return os.path.join(dp, name)

OG = find("Orthogroups.tsv", os.path.join(WORK, "05_orthology"))
MAP = os.path.join(WORK, "05_orthology", "aratha_protein2gene.tsv")

# locus -> protein accessions
locus2prot = defaultdict(list)
with open(MAP) as fh:
    next(fh)
    for ln in fh:
        p, g, sym, loc = ln.rstrip("\n").split("\t")
        if loc:
            locus2prot[loc].append(p)

# protein accession -> orthogroup
acc2og, og_row = {}, {}
with open(OG) as fh:
    species = fh.readline().rstrip("\n").split("\t")[1:]
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        og_row[p[0]] = dict(zip(species, p[1:]))
        for cell in p[1:]:
            for it in cell.split(", "):
                it = it.strip()
                if it:
                    acc2og[it.split("|")[-1]] = p[0]

def n(cell):
    return len([x for x in cell.split(", ") if x.strip()]) if cell else 0

print(f"\nspecies: {', '.join(species)}\n")
print(f"{'step':<16}{'gene':<10}{'locus':<12}{'orthogroup':<14}{'taro':>5}{'arab':>6}")
print("-" * 68)

rows, psy = [], None
for step, gene, locus in PATHWAY:
    prots = locus2prot.get(locus, [])
    ogs = sorted({acc2og[p] for p in prots if p in acc2og})
    if not ogs:
        print(f"{step:<16}{gene:<10}{locus:<12}{'not found':<14}")
        continue
    for og in ogs:
        r = og_row[og]
        c = {s: n(r.get(s, "")) for s in species}
        mark = "" if c["aratha"] == 1 else "  <- check"
        print(f"{step:<16}{gene:<10}{locus:<12}{og:<14}"
              f"{c['colesc']:>5}{c['aratha']:>6}{mark}")
        rows.append([step, gene, locus, og] + [c[s] for s in species])
        if gene == "PSY":
            psy = (og, c, r)

print("\n" + "=" * 68)
if psy:
    og, c, r = psy
    print(f"PSY  —  orthogroup {og}\n")
    for s in species:
        print(f"   {s:<10} {c[s]:>3}")
    print()
    print(f"taro phytoene synthase genes: {c['colesc']}")
    ok = c["aratha"] == 1
    print(f"sanity check: Arabidopsis = {c['aratha']} "
          f"({'pass, expected 1' if ok else 'FAIL — investigate before trusting the taro count'})")
    print("\ntaro members:")
    for m in r["colesc"].split(", "):
        if m.strip():
            print(f"   {m.strip()}")
else:
    print("PSY orthogroup not located")

out = os.path.join(CODE, "results", "tables", "pathway_orthogroups.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["step", "gene", "at_locus", "orthogroup"] + species)
    w.writerows(rows)
print(f"\nwritten: {out}")
