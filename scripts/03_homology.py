#!/usr/bin/env python3
"""
03_homology.py — candidate recovery by reciprocal best hit, with a sensitivity
grid over the detection thresholds.

REWRITTEN. The previous version is in scripts/archive/. Four things were wrong
with it, and all four were visible in its own output.

1. COVERAGE WAS NOT COVERAGE.

   It computed alignment_length / query_length. Alignment length counts gap
   columns, so the ratio exceeds 1 whenever the alignment contains insertions
   relative to the query. The output reported LCYB at 102%, LCYE at 103%,
   ORANGE at 107% and fibrillin at 106%, which is the arithmetic announcing
   itself. Coverage is now DIAMOND's qcovhsp: the percentage of query residues
   inside the aligned block. Both are kept in the output table so the change is
   auditable, and the grid now sweeps real coverage.

2. RECIPROCAL BEST HIT WAS TESTED AGAINST THE ANCHOR ALONE.

   A taro gene qualified only if its best Arabidopsis hit was the exact anchor
   locus. For a target that stands for a family, that is too strict: a taro
   NCED whose closest Arabidopsis relative is NCED5 fails the test even though
   NCED5 is in the family the target measures. Two levels are now recorded:

     rbh_anchor  reverse best hit is the anchor locus
     rbh_family  reverse best hit is the anchor OR a declared family member

   rbh_family is what the grid counts, because the target is the family. Both
   columns are written out, so the stricter count is recoverable.

3. SUB-THRESHOLD HITS WERE DISCARDED IN MEMORY.

   Hits failing the default thresholds were dropped before anything was
   written. The two PSY fragments Ces12496 and Ces12497, which are a recorded
   finding in LOGBOOK.md, are short enough to fall under them, and they do not
   appear anywhere in the previous output. A finding cannot be allowed to
   disappear because a threshold moved. Every forward hit at e < 1e-5 is now
   written to homology_hits_all.tsv regardless of whether it passes anything,
   and the reverse best hit of each is recorded next to it.

4. IT CRASHED, AND THAT ONE IS MINE.

   patch_schema.sh renamed the anchor columns by string replacement and missed
   raw.get(a["gene"]) at line 179, so the run died after printing both tables
   but before writing homology_sensitivity.tsv. Patching a script by replacing
   strings in it is how that happens. This version reads the schema in one
   place and passes a target object around instead.

METHOD, and what it does not establish.

  Reciprocal best hit identifies likely taro counterparts of Arabidopsis
  pathway genes. That is DETECTION. It is a lower bound on family size, because
  a paralog divergent enough that its closest Arabidopsis relative lies outside
  the declared family set never appears.

  Copy number is established later, by phylogeny (04_) and genomic coordinates
  (05_). Nothing in this file is a copy-number claim.

  The grid sweeps alignment length and query coverage to ask whether the
  CANDIDATE SET is stably recovered. None of these thresholds has a principled
  justification; they are conventional, which is the reason for sweeping them
  rather than defending them.
"""
import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REFS = Path(os.environ["TARO_REFS"])
WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
D = REFS / "proteomes"
W = WORK / "03_homology"
MAP = WORK / "05_orthology" / "aratha_protein2gene.tsv"
TBL = CODE / "results" / "tables"

MIN_ALN = [100, 150, 200]        # alignment length, aa
MIN_COV = [50, 70, 80]           # query coverage, percent
DEFAULT_ALN, DEFAULT_COV = 150, 70
MAX_TARGETS = 50                 # forward hits kept per anchor

W.mkdir(parents=True, exist_ok=True)
TBL.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------------ inputs
def read_fasta(p):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].strip()
            d[k] = []
        else:
            d[k].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


targets = list(csv.DictReader(open(TBL / "pathway_anchors.tsv"), delimiter="\t"))
if "family_members" not in targets[0]:
    sys.exit("pathway_anchors.tsv has no family_members column; "
             "run 02e_group_families.py first")

locus2acc = defaultdict(list)
acc2locus = {}
acc2sym = {}
for ln in open(MAP):
    if ln.startswith("protein\t"):
        continue
    p, g, sym, loc = ln.rstrip("\n").split("\t")
    if loc:
        locus2acc[loc.upper()].append(p)
        acc2locus[p] = loc.upper()
        acc2sym[p] = sym

ara = read_fasta(D / "orthofinder_input" / "aratha.fa")
taro = read_fasta(D / "orthofinder_input" / "colesc.fa")
ara_by_acc = {k.split("|")[-1]: (k, v) for k, v in ara.items()}


def dmnd(qfile, db, out, maxt):
    subprocess.run(
        ["diamond", "blastp", "-q", str(qfile), "-d", str(db), "-o", str(out),
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qcovhsp",
         "scovhsp", "evalue", "bitscore",
         "--max-target-seqs", str(maxt), "--evalue", "1e-5",
         "--threads", os.environ.get("TARO_THREADS", "4"), "--quiet"],
        check=True, cwd=W)
    rows = []
    p = W / out
    if p.exists():
        for ln in open(p):
            f = ln.rstrip("\n").split("\t")
            rows.append(dict(q=f[0], s=f[1], pid=float(f[2]), aln=int(f[3]),
                             qcov=float(f[4]), scov=float(f[5]),
                             ev=f[6], bits=float(f[7])))
    return rows


# ------------------------------------------------- pass 1, anchors into taro
resolved, unresolved = [], []
with open(W / "anchors.faa", "w") as fh:
    for t in targets:
        loc = t["anchor_locus"].upper()
        acc = next((x for x in locus2acc.get(loc, []) if x in ara_by_acc), None)
        if acc is None:
            unresolved.append(t)
            continue
        hdr, seq = ara_by_acc[acc]
        t["_acc"] = acc
        t["_qlen"] = len(seq)
        t["_famset"] = {loc} | {x.upper() for x in t["family_members"].split(";") if x}
        resolved.append(t)
        fh.write(f">{acc}\n{seq}\n")

if unresolved:
    print("anchors absent from the Arabidopsis proteome, no search possible:")
    for t in unresolved:
        print(f"  {t['target']:<16} {t['anchor_locus']}")
    print()

for db, src in (("colesc", "colesc.fa"), ("aratha", "aratha.fa")):
    if not (W / f"{db}.dmnd").exists():
        subprocess.run(["diamond", "makedb", "--in",
                        str(D / "orthofinder_input" / src),
                        "-d", db, "--quiet"], check=True, cwd=W)

print(f"forward search: {len(resolved)} anchors into the taro proteome")
fwd = dmnd("anchors.faa", "colesc", "fwd.tsv", MAX_TARGETS)
by_anchor = defaultdict(list)
for h in fwd:
    by_anchor[h["q"]].append(h)

# ------------------------------------------- pass 2, those taro genes back
taro_hits = sorted({h["s"] for h in fwd})
with open(W / "taro_hits.faa", "w") as fh:
    for sid in taro_hits:
        fh.write(f">{sid}\n{taro[sid]}\n")

print(f"reverse search: {len(taro_hits)} unique taro genes into Arabidopsis\n")
rev = dmnd("taro_hits.faa", "aratha", "rev.tsv", 5)
rev_best = {}
for h in rev:
    cur = rev_best.get(h["q"])
    if cur is None or h["bits"] > cur["bits"]:
        rev_best[h["q"]] = h

rev_locus, rev_sym = {}, {}
for sid, h in rev_best.items():
    acc = h["s"].split("|")[-1]
    rev_locus[sid] = acc2locus.get(acc, "")
    rev_sym[sid] = acc2sym.get(acc, "")


# ------------------------------------------------------------------- report
def call_of(h, rbh_fam, rbh_anc):
    if rbh_fam and h["qcov"] >= DEFAULT_COV:
        return "ortholog (RBH)" if rbh_anc else "ortholog (RBH, family)"
    if rbh_fam:
        return "partial (RBH)"
    if h["qcov"] >= DEFAULT_COV:
        return "homolog"
    return "fragment"


W_T, W_G = 16, 12
print("=" * 118)
print("CANDIDATE RECOVERY  —  reciprocal best hit")
print(f"  default thresholds: alignment >= {DEFAULT_ALN} aa, query coverage "
      f">= {DEFAULT_COV}%")
print("  coverage is DIAMOND qcovhsp, the percent of query residues in the "
      "aligned block")
print("=" * 118)
print(f"\n{'pathway':<10}{'target':<{W_T}}{'taro gene':<{W_G}}{'id%':>6}"
      f"{'qcov':>6}{'aln':>6}  {'evalue':<11} {'rev best hit':<22} call")
print("-" * 118)

all_rows, pass_rows = [], []
raw = defaultdict(list)

for t in resolved:
    hits = sorted(by_anchor.get(t["_acc"], []), key=lambda h: -h["bits"])
    if not hits:
        print(f"{t['pathway']:<10}{t['target']:<{W_T}}{'none':<{W_G}}"
              f"  no hit at e < 1e-5")
        continue

    first = True
    for h in hits:
        g = h["s"].split("|")[-1]
        rl = rev_locus.get(h["s"], "")
        rs = rev_sym.get(h["s"], "")
        rbh_anc = rl == t["anchor_locus"].upper()
        rbh_fam = rl in t["_famset"]
        c = call_of(h, rbh_fam, rbh_anc)

        rec = dict(pathway=t["pathway"], target=t["target"],
                   anchor_locus=t["anchor_locus"], claim_kind=t["claim_kind"],
                   family_group=t["family_group"], taro_gene=g,
                   identity=f"{h['pid']:.1f}", qcov_pct=f"{h['qcov']:.0f}",
                   scov_pct=f"{h['scov']:.0f}", aln_aa=h["aln"],
                   aln_over_qlen_pct=f"{100.0 * h['aln'] / t['_qlen']:.0f}",
                   evalue=h["ev"], rev_best_locus=rl, rev_best_symbol=rs,
                   rbh_anchor=int(rbh_anc), rbh_family=int(rbh_fam), call=c)
        all_rows.append(rec)
        raw[t["target"]].append((h["aln"], h["qcov"], rbh_fam))

        if h["aln"] < DEFAULT_ALN and h["qcov"] < DEFAULT_COV:
            continue
        rb = f"{rs or rl or '-'}"[:21]
        print(f"{(t['pathway'] if first else ''):<10}"
              f"{(t['target'] if first else ''):<{W_T}}"
              f"{g:<{W_G}}{h['pid']:>5.1f}{h['qcov']:>6.0f}{h['aln']:>6}  "
              f"{h['ev']:<11} {rb:<22} {c}")
        first = False
        pass_rows.append(rec)
    if first:
        print(f"{t['pathway']:<10}{t['target']:<{W_T}}"
              f"{'(all sub-threshold)':<{W_G}}")

# ------------------------------------------------------------------- grid
print()
print("=" * 118)
print("SENSITIVITY  —  is the candidate set stable across the grid?")
print("=" * 118)
print("\nCounts are reciprocal best hits to the family, not copy number. Copy")
print("number comes from 04_ (phylogeny) and 05_ (genomic locus).\n")

cells = [(a, c) for a in MIN_ALN for c in MIN_COV]
print(f"{'target':<{W_T}}" + "".join(f"{a}/{c}%".rjust(10) for a, c in cells)
      + "   stable")
print("-" * 118)

unstable = []
grid = {}
for t in resolved:
    hits = raw.get(t["target"], [])
    counts = [sum(1 for aln, cov, rf in hits if rf and aln >= a and cov >= c)
              for a, c in cells]
    grid[t["target"]] = counts
    ok = len(set(counts)) == 1
    print(f"{t['target']:<{W_T}}" + "".join(f"{n:>10d}" for n in counts)
          + f"   {'yes' if ok else 'NO'}")
    if not ok:
        unstable.append((t["target"], counts))

print()
if unstable:
    print(f"  {len(unstable)} of {len(resolved)} candidate sets move across the grid:")
    for name, c in unstable:
        print(f"    {name:<16} {min(c)} to {max(c)}")
    print()
    print("  These are marked unstable in the master table. A candidate count")
    print("  that depends on an arbitrary cutoff is a property of the cutoff.")
else:
    print("  Every candidate set is identical across all nine combinations.")

# -------------------------------------------------- family group summary
print()
print("=" * 118)
print("SHARED FAMILIES  —  where one pool of taro genes feeds several targets")
print("=" * 118)
groups = defaultdict(list)
for t in resolved:
    if t["family_group"]:
        groups[t["family_group"]].append(t)
for g, ts in sorted(groups.items()):
    pool = set()
    for t in ts:
        pool |= {r["taro_gene"] for r in all_rows if r["target"] == t["target"]}
    print(f"\n  {g}: {len(pool)} distinct taro genes recovered across "
          f"{len(ts)} targets")
    for t in ts:
        n = sum(1 for aln, cov, rf in raw.get(t["target"], [])
                if rf and aln >= DEFAULT_ALN and cov >= DEFAULT_COV)
        print(f"    {t['target']:<16} {n} RBH   ({t['claim_kind']})")
    cn = [t for t in ts if t["claim_kind"] == "copy-number"]
    if len(cn) > 1:
        print(f"    -> 04_ builds ONE tree for this group. The partition above "
              f"is an RBH result, not a premise.")

# ------------------------------------------------------------------ outputs
F_ALL = ["pathway", "target", "anchor_locus", "claim_kind", "family_group",
         "taro_gene", "identity", "qcov_pct", "scov_pct", "aln_aa",
         "aln_over_qlen_pct", "evalue", "rev_best_locus", "rev_best_symbol",
         "rbh_anchor", "rbh_family", "call"]

for path, data, label in (
        (TBL / "homology_hits_all.tsv", all_rows, "every forward hit at e < 1e-5"),
        (TBL / "homology_candidates.tsv", pass_rows, "hits passing the defaults")):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=F_ALL)
        w.writeheader()
        w.writerows(data)
    print(f"\nwritten: {path}   ({len(data)} rows, {label})")

gp = TBL / "homology_sensitivity.tsv"
with open(gp, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["target", "claim_kind"]
               + [f"aln{a}_cov{c}" for a, c in cells] + ["stable"])
    for t in resolved:
        c = grid[t["target"]]
        w.writerow([t["target"], t["claim_kind"]] + c + [len(set(c)) == 1])
print(f"written: {gp}")

print("""
next:  bash scripts/04_families.sh

  Read homology_hits_all.tsv before that. It holds the sub-threshold hits the
  default view suppresses, including anything short enough to be a split gene
  model, and the rev_best_symbol column says which Arabidopsis gene each
  non-reciprocal hit actually belongs to.
""")
