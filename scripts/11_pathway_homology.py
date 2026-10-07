#!/usr/bin/env python3
"""
Carotenoid and MEP pathway inventory in taro, by reciprocal homology search.

For each Arabidopsis pathway gene:
  forward   Arabidopsis protein -> taro proteome
  reverse   each taro hit       -> Arabidopsis proteome

A taro gene is called a confident ortholog when it is a reciprocal best hit
and aligns over at least 70% of the query length. Shorter or weaker matches
are reported as partial, because gene-model fragmentation produces those and
they should not be counted as paralogs without inspection.
"""
import os, subprocess, sys, csv
from collections import defaultdict

REFS = os.environ["TARO_REFS"]; WORK = os.environ["TARO_WORK"]
CODE = os.environ["TARO_CODE"]
D = f"{REFS}/proteomes"
W = f"{WORK}/05_orthology/pathway_homology"
MAP = f"{WORK}/05_orthology/aratha_protein2gene.tsv"

PATHWAY = [
    ("MEP",       "DXS",     "AT4G15560"), ("MEP",       "DXR",     "AT5G62790"),
    ("MEP",       "MCT",     "AT2G02500"), ("MEP",       "CMK",     "AT2G26930"),
    ("MEP",       "MDS",     "AT1G63970"), ("MEP",       "HDS",     "AT5G60600"),
    ("MEP",       "HDR",     "AT4G34350"), ("backbone",  "GGPPS11", "AT4G36810"),
    ("core",      "PSY",     "AT5G17230"), ("core",      "PDS3",    "AT4G14210"),
    ("core",      "Z-ISO",   "AT1G10830"), ("core",      "ZDS",     "AT3G04870"),
    ("core",      "CRTISO",  "AT1G06820"), ("cyclase",   "LCYB",    "AT3G10230"),
    ("cyclase",   "LCYE",    "AT5G57030"), ("xanth",     "CYP97A3", "AT1G31800"),
    ("xanth",     "CYP97C1", "AT3G53130"), ("xanth",     "BCH1",    "AT4G25700"),
    ("xanth",     "BCH2",    "AT5G52570"), ("xanth",     "ZEP",     "AT5G67030"),
    ("xanth",     "VDE",     "AT1G08550"), ("xanth",     "NSX",     "AT1G67080"),
    ("degrade",   "CCD1",    "AT3G63520"), ("degrade",   "CCD4",    "AT4G19170"),
    ("degrade",   "NCED3",   "AT3G14440"), ("sink",      "OR",      "AT5G61670"),
    ("sink",      "OR-like", "AT5G06130"), ("sink",      "FBN1a",   "AT4G04020"),
]

# locus -> accessions, and the subset actually present in the proteome
locus2acc = defaultdict(list)
acc2locus = {}
for ln in open(MAP):
    if ln.startswith("protein\t"): continue
    p, g, sym, loc = ln.rstrip("\n").split("\t")
    if loc:
        locus2acc[loc].append(p); acc2locus[p] = loc

def read_fasta(path):
    seqs, name = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            name = ln[1:].strip(); seqs[name] = []
        else: seqs[name].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}

ara = read_fasta(f"{D}/aratha.primary.faa")
ara_by_acc = {k.split("|")[-1]: (k, v) for k, v in ara.items()}

def dmnd(query_file, db, out):
    subprocess.run(["diamond", "blastp", "-q", query_file, "-d", db, "-o", out,
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "evalue", "bitscore", "--max-target-seqs", "20",
                    "--evalue", "1e-5", "--quiet"], check=True)
    rows = []
    if os.path.exists(out):
        for ln in open(out):
            f = ln.rstrip("\n").split("\t")
            rows.append((f[0], f[1], float(f[2]), int(f[3]), f[4], float(f[5])))
    return rows

os.chdir(W)
results = []
print(f"\n{'step':<9}{'gene':<9}{'AT locus':<12}{'taro gene':<12}"
      f"{'id%':>6}{'cov':>6}  {'evalue':<11} call")
print("-" * 82)

for step, gene, locus in PATHWAY:
    acc = next((a for a in locus2acc.get(locus, []) if a in ara_by_acc), None)
    if not acc:
        print(f"{step:<9}{gene:<9}{locus:<12}{'-':<12}{'':>6}{'':>6}  "
              f"{'':<11} query not in proteome")
        continue
    hdr, seq = ara_by_acc[acc]
    qlen = len(seq)
    with open("q.faa", "w") as fh:
        fh.write(f">{hdr}\n{seq}\n")

    fwd = dmnd("q.faa", "colesc", "fwd.tsv")
    if not fwd:
        print(f"{step:<9}{gene:<9}{locus:<12}{'none':<12}{'':>6}{'':>6}  "
              f"{'':<11} no hit")
        results.append([step, gene, locus, "", 0, 0, "", "absent"])
        continue

    # reverse search each taro hit
    taro = read_fasta(f"{D}/orthofinder_input/colesc.fa")
    for i, (_, sid, pid, alen, ev, bits) in enumerate(fwd):
        with open("t.faa", "w") as fh:
            fh.write(f">{sid}\n{taro[sid]}\n")
        rev = dmnd("t.faa", "aratha", "rev.tsv")
        rbh = bool(rev) and rev[0][1].split("|")[-1] == acc
        cov = alen / qlen
        if rbh and cov >= 0.70:   call = "ortholog (RBH)"
        elif rbh:                 call = "partial (RBH)"
        elif cov >= 0.70:         call = "homolog"
        else:                     call = "fragment"
        g = sid.split("|")[-1]
        print(f"{step if i==0 else '':<9}{gene if i==0 else '':<9}"
              f"{locus if i==0 else '':<12}{g:<12}{pid:>5.1f}{cov*100:>6.0f}  "
              f"{ev:<11} {call}")
        results.append([step, gene, locus, g, pid, round(cov*100), ev, call])
        if i >= 4: break

out = f"{CODE}/results/tables/pathway_homology.tsv"
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["step","gene","at_locus","taro_gene","identity","coverage_pct","evalue","call"])
    w.writerows(results)

print("\n" + "=" * 82)
print("confident orthologs per pathway gene (RBH, coverage >= 70%):")
cnt = defaultdict(list)
for r in results:
    if r[7] == "ortholog (RBH)":
        cnt[r[1]].append(r[3])
for _, gene, _ in PATHWAY:
    if gene in cnt:
        print(f"  {gene:<9} {len(cnt[gene])}   {', '.join(cnt[gene])}")
print(f"\nwritten: {out}")
