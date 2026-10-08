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

   (Fixed, and then half-undone by the printing. See 5.)

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

5. A RECIPROCAL BEST HIT COULD STILL BE HIDDEN FROM THE PRINTED TABLE.

   homology_hits_all.tsv did keep the two PSY fragments, exactly as intended:

     Ces12497  81.0% identity  qcov 34%  aln 147 aa  e 2.04e-81  RBH to PSY
     Ces12496  57.4% identity  qcov 26%  aln 115 aa  e 4.11e-35  RBH to PSY

   But the print filter suppressed both, because each fails BOTH defaults, and
   so the terminal output showed PSY recovering one gene. Writing a finding to
   a file and then hiding it from the only view anyone reads is not much better
   than discarding it. The filter now never suppresses a reciprocal best hit,
   whatever its coverage, and prints the call beside it.

6. QUERY COORDINATES WERE NOT RECORDED, SO SPLIT MODELS WERE NOT TESTABLE.

   Two partial hits to the same anchor are either one gene split across two
   annotation records or two genuinely separate partial genes, and coverage
   alone cannot tell them apart: 26% and 34% sum to 60% whether the two pieces
   cover different parts of the query or the same part twice. The forward
   search now records qstart, qend, sstart and send, and a section at the end
   lists pairs of reciprocal hits whose query ranges barely overlap. That is
   the sequence-side evidence for a split model. The genome-side evidence,
   adjacency and strand, comes from the GFF3 in 05_.

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
         "scovhsp", "evalue", "bitscore", "qstart", "qend", "sstart", "send",
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
                             ev=f[6], bits=float(f[7]),
                             qs=int(f[8]), qe=int(f[9]),
                             ss=int(f[10]), se=int(f[11])))
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

# Rebuild a DIAMOND database when it is missing OR older than the FASTA it was
# built from. Testing only for existence is the bug that sits in
# 04_families.sh, where a database built from the nine-species set would be
# silently reused against the current seventeen. Same shape as the KEGG query
# returning zero and the single-tip rooting returning all zeros: no error, wrong
# answer. A timestamp comparison is not a checksum, but it catches the case that
# actually occurs, which is a regenerated proteome and a stale database.
for db, src in (("colesc", "colesc.fa"), ("aratha", "aratha.fa")):
    fa = D / "orthofinder_input" / src
    dmnd_path = W / f"{db}.dmnd"
    stale = dmnd_path.exists() and dmnd_path.stat().st_mtime < fa.stat().st_mtime
    if not dmnd_path.exists() or stale:
        if stale:
            print(f"  {db}.dmnd is older than {src}; rebuilding")
        subprocess.run(["diamond", "makedb", "--in", str(fa),
                        "-d", db, "--quiet"], check=True, cwd=W)
    else:
        print(f"  {db}.dmnd reused, newer than {src}")

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
rbh_hits = defaultdict(list)   # target -> reciprocal hits, for split detection

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
                   q_start=h["qs"], q_end=h["qe"],
                   s_start=h["ss"], s_end=h["se"],
                   evalue=h["ev"], rev_best_locus=rl, rev_best_symbol=rs,
                   rbh_anchor=int(rbh_anc), rbh_family=int(rbh_fam), call=c)
        all_rows.append(rec)
        raw[t["target"]].append((h["aln"], h["qcov"], rbh_fam))
        if rbh_fam:
            rbh_hits[t["target"]].append(rec)

        # A reciprocal best hit is always printed. Suppressing one is how the
        # two PSY fragments came to look absent when they were in the file all
        # along.
        if not rbh_fam and h["aln"] < DEFAULT_ALN and h["qcov"] < DEFAULT_COV:
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

# ------------------------------------------- split gene model candidates
print()
print("=" * 118)
print("SPLIT GENE MODEL CANDIDATES  —  sequence-side evidence only")
print("=" * 118)
print("""
Two reciprocal hits to the same anchor are either one gene broken across two
annotation records or two separate partial genes. Coverage cannot tell them
apart, because 26% and 34% sum to 60% whether the pieces cover different parts
of the query or the same part twice. What distinguishes them is WHERE on the
query each piece aligns.

Listed below: pairs of reciprocal hits to the same target that between them
cover much more of the anchor than either does alone, and whose aligned query
ranges barely overlap. That is consistent with a split model and not proof of
one. The genome-side evidence, whether the two models sit adjacent on the same
strand with nothing annotated between them, comes from the GFF3 in 05_.
""")

SPLIT_MAX_OVERLAP = 0.25      # of the shorter aligned range
SPLIT_MIN_GAIN = 15.0         # combined coverage must exceed the better single
                              # hit by this many percentage points

found_any = False
for t in resolved:
    hits = rbh_hits.get(t["target"], [])
    if len(hits) < 2:
        continue
    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            a, b = hits[i], hits[j]
            a_s, a_e = int(a["q_start"]), int(a["q_end"])
            b_s, b_e = int(b["q_start"]), int(b["q_end"])
            la, lb = a_e - a_s + 1, b_e - b_s + 1
            ov = max(0, min(a_e, b_e) - max(a_s, b_s) + 1)
            ov_frac = ov / min(la, lb) if min(la, lb) else 1.0
            ca, cb = float(a["qcov_pct"]), float(b["qcov_pct"])
            union = la + lb - ov
            comb = 100.0 * union / t["_qlen"]
            if ov_frac > SPLIT_MAX_OVERLAP:
                continue
            if comb - max(ca, cb) < SPLIT_MIN_GAIN:
                continue
            found_any = True
            # numerically adjacent gene identifiers are suggestive, not evidence
            def num(g):
                d = "".join(ch for ch in g if ch.isdigit())
                return int(d) if d else -1
            adj = abs(num(a["taro_gene"]) - num(b["taro_gene"])) == 1
            print(f"  {t['target']}   anchor {t['anchor_locus']}  "
                  f"({t['_qlen']} aa)")
            for h, l in ((a, la), (b, lb)):
                print(f"    {h['taro_gene']:<12} query {int(h['q_start']):>4}-"
                      f"{int(h['q_end']):<4} ({l:>3} aa)  "
                      f"id {h['identity']:>5}%  qcov {h['qcov_pct']:>3}%  "
                      f"e {h['evalue']}")
            print(f"    overlap of the shorter range : {ov} aa ({ov_frac*100:.0f}%)")
            print(f"    coverage apart               : {ca:.0f}% and {cb:.0f}%")
            print(f"    coverage together            : {comb:.0f}%")
            print(f"    gene identifiers adjacent    : "
                  f"{'yes' if adj else 'no'}")
            print(f"    -> 05_ must check: same scaffold, same strand, nothing")
            print(f"       annotated between. Until then this is one candidate")
            print(f"       locus reported as two models, not two copies.")
            print()

if not found_any:
    print("  No pair of reciprocal hits splits the anchor between them.")

# ------------------------------------------------------------------ outputs
F_ALL = ["pathway", "target", "anchor_locus", "claim_kind", "family_group",
         "taro_gene", "identity", "qcov_pct", "scov_pct", "aln_aa",
         "aln_over_qlen_pct", "q_start", "q_end", "s_start", "s_end",
         "evalue", "rev_best_locus", "rev_best_symbol",
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
