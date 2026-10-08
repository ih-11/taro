#!/usr/bin/env python3
"""Reciprocal best hits for the pathway anchors, swept across the grid."""
import csv
import os
import subprocess
from collections import defaultdict
from pathlib import Path

REFS = Path(os.environ["TARO_REFS"])
WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
D = REFS / "proteomes"
W = WORK / "03_homology"
MAP = WORK / "05_orthology" / "aratha_protein2gene.tsv"

MIN_ALN = [100, 150, 200]        # alignment length, aa
MIN_COV = [0.50, 0.70, 0.80]     # query coverage
DEFAULT_ALN, DEFAULT_COV = 150, 0.70

anchors = []
with open(CODE / "results" / "tables" / "pathway_anchors.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        anchors.append(r)

locus2acc = defaultdict(list)
for ln in open(MAP):
    if ln.startswith("protein\t"):
        continue
    p, g, sym, loc = ln.rstrip("\n").split("\t")
    if loc:
        locus2acc[loc.upper()].append(p)


def read_fasta(p):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].strip(); d[k] = []
        else:
            d[k].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


ara = read_fasta(D / "orthofinder_input" / "aratha.fa")
taro = read_fasta(D / "orthofinder_input" / "colesc.fa")
ara_by_acc = {k.split("|")[-1]: (k, v) for k, v in ara.items()}


def dmnd(q, db, out):
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db), "-o", str(out),
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "evalue", "bitscore", "--max-target-seqs", "25",
                    "--evalue", "1e-5", "--quiet"], check=True, cwd=W)
    rows = []
    p = W / out
    if p.exists():
        for ln in open(p):
            f = ln.rstrip("\n").split("\t")
            rows.append((f[0], f[1], float(f[2]), int(f[3]), f[4], float(f[5])))
    return rows


print("=" * 86)
print("CANDIDATE RECOVERY  —  reciprocal best hit, default thresholds")
print(f"  alignment >= {DEFAULT_ALN} aa, coverage >= {DEFAULT_COV:.0%}")
print("=" * 86)
print(f"\n{'step':<10}{'gene':<10}{'taro gene':<13}{'id%':>6}{'cov':>6}"
      f"{'aln':>6}  {'evalue':<11} call")
print("-" * 86)

rows = []
raw = defaultdict(list)     # gene -> all forward hits, for the grid

for a in anchors:
    gene, locus = a["gene"], a["at_locus"].upper()
    acc = next((x for x in locus2acc.get(locus, []) if x in ara_by_acc), None)
    if acc is None:
        print(f"{a['step']:<10}{gene:<10}{'-':<13}  anchor not in proteome")
        continue

    hdr, seq = ara_by_acc[acc]
    qlen = len(seq)
    (W / "q.faa").write_text(f">{hdr}\n{seq}\n")
    fwd = dmnd("q.faa", "colesc", "fwd.tsv")

    if not fwd:
        print(f"{a['step']:<10}{gene:<10}{'none':<13}  no hit at e<1e-5")
        rows.append(dict(step=a["step"], gene=gene, at_locus=locus,
                         claim_kind=a["claim_kind"], taro_gene="",
                         identity="", coverage_pct="", aln="", evalue="",
                         call="absent"))
        continue

    first = True
    for _, sid, pid, alen, ev, bits in fwd[:12]:
        (W / "t.faa").write_text(f">{sid}\n{taro[sid]}\n")
        rev = dmnd("t.faa", "aratha", "rev.tsv")
        rbh = bool(rev) and rev[0][1].split("|")[-1] == acc
        cov = alen / qlen
        raw[gene].append((sid, pid, alen, cov, ev, rbh))

        if alen < DEFAULT_ALN and cov < DEFAULT_COV:
            continue
        if rbh and cov >= DEFAULT_COV:
            call = "ortholog (RBH)"
        elif rbh:
            call = "partial (RBH)"
        elif cov >= DEFAULT_COV:
            call = "homolog"
        else:
            call = "fragment"

        g = sid.split("|")[-1]
        print(f"{a['step'] if first else '':<10}{gene if first else '':<10}"
              f"{g:<13}{pid:>5.1f}{cov*100:>6.0f}{alen:>6}  {ev:<11} {call}")
        first = False
        rows.append(dict(step=a["step"], gene=gene, at_locus=locus,
                         claim_kind=a["claim_kind"], taro_gene=g,
                         identity=f"{pid:.1f}", coverage_pct=f"{cov*100:.0f}",
                         aln=alen, evalue=ev, call=call))

# ---------------------------------------------------------------- grid
print()
print("=" * 86)
print("SENSITIVITY  —  is the candidate set stable across the grid?")
print("=" * 86)
print("\nThis asks about candidate recovery, not copy number. Copy number comes")
print("from 04_ (phylogeny) and 05_ (genomic locus).\n")

header = "  ".join(f"{a}aa/{int(c*100)}%" for a in MIN_ALN for c in MIN_COV)
print(f"{'gene':<10}  {header}   stable")
print("-" * 86)

unstable = []
for a in anchors:
    gene = a["gene"]
    hits = raw.get(gene, [])
    counts = []
    for mina in MIN_ALN:
        for minc in MIN_COV:
            n = sum(1 for _, _, alen, cov, _, rbh in hits
                    if rbh and alen >= mina and cov >= minc)
            counts.append(n)
    cells = "  ".join(f"{n:>7d}" for n in counts)
    stable = len(set(counts)) == 1
    print(f"{gene:<10}  {cells}   {'yes' if stable else 'NO'}")
    if not stable:
        unstable.append((gene, counts))

print()
if unstable:
    print(f"  {len(unstable)} of {len(anchors)} candidate sets move across the grid:")
    for gene, c in unstable:
        print(f"    {gene:<10} {min(c)} to {max(c)}")
    print()
    print("  These are reported as unstable in the master table. A candidate")
    print("  count that depends on an arbitrary threshold is a property of the")
    print("  threshold and belongs in the paper as such.")
else:
    print("  Every candidate set is identical across all nine threshold")
    print("  combinations. The default thresholds are not doing the work.")

out = CODE / "results" / "tables" / "homology_candidates.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["step", "gene", "at_locus", "claim_kind",
                                   "taro_gene", "identity", "coverage_pct",
                                   "aln", "evalue", "call"])
    w.writeheader()
    w.writerows(rows)

grid_out = CODE / "results" / "tables" / "homology_sensitivity.tsv"
with open(grid_out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene"] + [f"aln{a}_cov{int(c*100)}"
                           for a in MIN_ALN for c in MIN_COV] + ["stable"])
    for a in anchors:
        hits = raw.get(a["gene"], [])
        counts = [sum(1 for _, _, alen, cov, _, rbh in hits
                      if rbh and alen >= mina and cov >= minc)
                  for mina in MIN_ALN for minc in MIN_COV]
        w.writerow([a["gene"]] + counts + [len(set(counts)) == 1])

print(f"\nwritten: {out}")
print(f"written: {grid_out}")
print("next:    bash scripts/04_families.sh")
