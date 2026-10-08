#!/usr/bin/env python3
"""
04_families.py — gene family trees for the copy-number targets.

FOUR THINGS WRONG WITH THE FIRST RUN, ALL VISIBLE IN ITS OWN OUTPUT

1. TREE B WAS NOT A FAMILY-ONLY TREE.

   The old collect() put family queries and root queries into one search, kept
   everything above one floor, then built tree B by removing only the
   Arabidopsis root SEQUENCES from that set:

       keep -= set(roots)

   So tree B still held every other species' squalene synthase, with just the
   two Arabidopsis SQS tips taken out. The tip counts said so: PSY tree A 72
   tips, tree B 70. A family-only tree should have lost the whole SQS clade,
   around eighteen tips, not two. GGPPS likewise, 147 and 145, a difference of
   exactly FPS1 and FPS2.

   That also explains the retention running backwards. PSY tree A retained 52%
   and tree B 34%, when removing a divergent outgroup should IMPROVE retention.
   Tree B had kept the divergent sequences and lost the two Arabidopsis anchors
   that helped align them. It was neither tree.

   Now every recovered sequence is assigned to whichever query it hits with the
   highest bitscore: best hit a family locus means family, best hit a root
   locus means root. Same reciprocal logic as 03_, applied inside the job. Tree
   B is the family set. Tree A is family plus root. The monophyly test gets a
   defined tip set instead of an assumed one.

2. THE ROOT GROUP WAS READ FROM THE FIRST ROW OF THE GROUP.

   root = ts[0]["outparalog"], which for the CCD group is CCD1, a detection
   target with nothing declared. So tree A was skipped and the group fell to
   midpoint rooting even though CCD4 and NCED both declare CCD7 and CCD8.
   Taking the first row of a group as the group's property is wrong whenever
   the rows differ. The root spec is now the union over the group, and a
   disagreement between members is printed rather than silently resolved.

3. ONE FLOOR WAS DOING TWO JOBS.

   A reference homolog from another species is scaffolding: it supplies the
   alignment columns everything else is placed against, so it must be
   near-complete in BOTH directions. Half a protein supplies half a column set
   and gaps for the rest. Seven GGPPS-family queries against a 454,114-sequence
   database at qcov >= 50% pulled in cis-prenyltransferases and polyprenyl
   synthases of very different lengths, which is how a family of 300 to 400
   residue proteins produced a 3441-column alignment retaining 8%.

   References now need qcov >= 70% AND scov >= 70%. Taro candidates keep the
   permissive floor declared in 03_, alignment >= 100 aa and qcov >= 50%,
   because they are the question rather than the scaffolding.

4. trimAl -automated1 PRODUCES RETENTION FIGURES THAT CANNOT BE COMPARED.

   -automated1 picks between gappyout and strict from each alignment's own
   statistics, so two alignments can be trimmed by two different methods and
   their retention percentages are not measuring the same thing. Reporting
   retention as a diagnostic while letting the trimming method float defeats
   the diagnostic. A fixed gap threshold is comparable, so trimming is now
   -gt 0.80, keeping columns where at least 80% of sequences have a residue.
   The -automated1 column count is still computed and reported alongside, so
   the change is auditable against the earlier figures.

WHAT A TREE HERE ESTABLISHES, AND WHAT IT DOES NOT

  Reciprocal best hit (03_) says a taro gene's closest Arabidopsis relative is
  inside a declared family. That is detection, a statement about a pairwise
  search, not about a clade.

  A tree says whether a taro gene falls inside the clade the Arabidopsis
  members of the family define. That is family assignment.

  Neither is a copy number. Two sequences can be assigned to one family and
  still be one gene counted twice, as a split annotation model or as two
  haplotigs of a 1.83%-heterozygous assembly. Copy number needs the independent
  genomic locus, which is 05_.

    Copy number = family assignment (phylogeny) + independent genomic locus
                  (coordinates). Presence = reciprocal homology with adequate
                  coverage. A detection count is not a copy number.

TWO TREES PER JOB

  Tree A is family plus root group, and asks one question: is the root group
  monophyletic? If yes, the family boundary can be placed from the tree. If no,
  it cannot, and nothing further is claimed from it.

  Tree B is the family alone, rooted on the most distant species present
  (Amborella where available). It exists for resolution, because a divergent
  root group costs the alignment columns close relatives need.

SUPPORT

  SH-aLRT >= 80 AND UFBoot >= 95, both required. A node meeting one and not the
  other is reported as below threshold.

STOPPING RULE

  If tree A and tree B disagree on a taro gene and neither node clears support,
  the assignment is recorded as not determinable and no third analysis is run.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REFS = Path(os.environ["TARO_REFS"])
WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
D = REFS / "proteomes" / "orthofinder_input"
W = WORK / "04_families"
TBL = CODE / "results" / "tables"
MAP = WORK / "05_orthology" / "aratha_protein2gene.tsv"

# taro candidates: the permissive cell of the grid declared in 03_
TREE_MIN_ALN = 100
TREE_MIN_COV = 50.0
# reference homologs from the other 16 species: near-complete, both directions
REF_MIN_COV = 70.0

TRIM_GT = 0.80
ALRT_MIN, BOOT_MIN = 80.0, 95.0
RETENTION_WARN = 0.40
THREADS = os.environ.get("TARO_THREADS", "4")

W.mkdir(parents=True, exist_ok=True)
TBL.mkdir(parents=True, exist_ok=True)
IQ = "iqtree2" if shutil.which("iqtree2") else "iqtree"


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], check=True, **kw)


def read_fasta(p):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].strip()
            d[k] = []
        elif k is not None:
            d[k].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


def aln_cols(p):
    f = read_fasta(p)
    return len(next(iter(f.values()))) if f else 0


def rule(n=100):
    print("-" * n)


# ======================================================================
# 1. the combined search database, validated rather than assumed
# ======================================================================
SPECIES = sorted(p.stem for p in D.glob("*.fa"))
ALL_FAA, ALL_DMND = W / "all.faa", W / "all.dmnd"
SIDECAR = W / "all.manifest"

want = "\n".join(SPECIES) + "\n"
have = SIDECAR.read_text() if SIDECAR.exists() else ""

rebuild, why = False, ""
if not ALL_FAA.exists() or ALL_FAA.stat().st_size == 0:
    rebuild, why = True, "all.faa missing or empty"
elif not ALL_DMND.exists():
    rebuild, why = True, "all.dmnd missing"
elif ALL_DMND.stat().st_size == 0:
    rebuild, why = True, "all.dmnd is ZERO BYTES, a previous makedb did not finish"
elif have != want:
    rebuild, why = True, "the species set changed since the database was built"
elif ALL_DMND.stat().st_mtime < max(p.stat().st_mtime for p in D.glob("*.fa")):
    rebuild, why = True, "a proteome is newer than the database"

print("=" * 100)
print("SEARCH DATABASE")
print("=" * 100)
print(f"  species in {D.name}: {len(SPECIES)}")
print(f"  {' '.join(SPECIES)}")
if rebuild:
    print(f"\n  rebuilding: {why}")
    with open(ALL_FAA, "w") as out:
        for s in SPECIES:
            out.write((D / f"{s}.fa").read_text())
    run(["diamond", "makedb", "--in", ALL_FAA, "-d", W / "all", "--quiet"])
    SIDECAR.write_text(want)
else:
    print("\n  database reused; species set and timestamps agree")

ALL = read_fasta(ALL_FAA)
print(f"  sequences: {len(ALL)}")
per_sp = defaultdict(int)
for k in ALL:
    per_sp[k.split("|")[0]] += 1
absent = [s for s in SPECIES if per_sp[s] == 0]
if absent:
    sys.exit(f"  species in {D} but absent from all.faa: {absent}")

# ======================================================================
# 2. identifier maps
# ======================================================================
locus2acc, acc2locus, acc2sym = defaultdict(list), {}, {}
for ln in open(MAP):
    if ln.startswith("protein\t"):
        continue
    p, g, sym, loc = ln.rstrip("\n").split("\t")
    if loc:
        locus2acc[loc.upper()].append(p)
        acc2locus[p] = loc.upper()
        acc2sym[p] = sym

ara_ids = {k.split("|")[-1]: k for k in ALL if k.startswith("aratha|")}


def ara_key(locus):
    for acc in locus2acc.get(locus.upper(), []):
        if acc in ara_ids:
            return ara_ids[acc]
    return None


def label(key):
    sp, gid = key.split("|", 1)
    if sp == "aratha":
        sym = acc2sym.get(gid, "")
        loc = acc2locus.get(gid, gid)
        return f"aratha_{(sym or loc).replace(' ', '').replace('/', '-')}_{loc}"
    return f"{sp}_{gid}"


def root_spec(ts, quiet=True):
    """Root loci declared by ANY target in the group, not just the first row."""
    specs, loci = set(), []
    for t in ts:
        sp = t["outparalog"].strip()
        if not sp:
            continue
        specs.add(sp)
        for item in sp.split(";"):
            if item:
                loc = item.split(":")[-1].upper()
                if loc not in loci:
                    loci.append(loc)
    if len(specs) > 1 and not quiet:
        print("    NOTE  members of this group declare different root groups:")
        for sp in sorted(specs):
            print(f"          {sp}")
        print(f"    using the union, {len(loci)} loci. If they disagree because")
        print("    the targets are not one family, the grouping is wrong rather")
        print("    than the rooting.")
    return loci


# ======================================================================
# 3. tree jobs
# ======================================================================
targets = list(csv.DictReader(open(TBL / "pathway_anchors.tsv"), delimiter="\t"))
jobs = {}
for t in targets:
    jobs.setdefault(t["family_group"] or t["target"], []).append(t)
jobs = {g: ts for g, ts in jobs.items()
        if any(t["claim_kind"] == "copy-number" for t in ts)}

print()
print("=" * 100)
print("TREE JOBS")
print("=" * 100)
for g, ts in sorted(jobs.items()):
    cn = [t["target"] for t in ts if t["claim_kind"] == "copy-number"]
    print(f"\n  {g}")
    for t in ts:
        fam = t["family_members"].replace(";", " ")
        print(f"    {t['target']:<16} {t['claim_kind']:<12} anchor {t['anchor_locus']}"
              + (f"  + {fam}" if fam else ""))
    rs = root_spec(ts, quiet=False)
    if rs:
        names = []
        for loc in rs:
            k = ara_key(loc)
            names.append(f"{acc2sym.get(k.split('|')[-1], loc) if k else '?'}:{loc}")
        print(f"    rooted on: {' '.join(names)}")
    else:
        print("    rooted on: MIDPOINT, nothing declared")
    if len(cn) > 1:
        print(f"    one tree; {len(cn)} copy-number targets share this family")

# ======================================================================
# 4. taro candidates from 03_
# ======================================================================
hits = list(csv.DictReader(open(TBL / "homology_hits_all.tsv"), delimiter="\t"))
cand, below = defaultdict(set), defaultdict(dict)
for h in hits:
    if h["rbh_family"] != "1":
        continue
    g = h["family_group"] or h["target"]
    if g not in jobs:
        continue
    if int(h["aln_aa"]) >= TREE_MIN_ALN and float(h["qcov_pct"]) >= TREE_MIN_COV:
        cand[g].add(h["taro_gene"])
    else:
        below[g].setdefault(h["taro_gene"],
                            (h["qcov_pct"], h["aln_aa"], h["target"]))

print()
print("=" * 100)
print("TARO SEQUENCES ENTERING THE TREES")
print("=" * 100)
print(f"\n  taro candidates : alignment >= {TREE_MIN_ALN} aa AND qcov >= "
      f"{TREE_MIN_COV:.0f}%   (the permissive cell of the 03_ grid)")
print(f"  reference tips  : qcov >= {REF_MIN_COV:.0f}% AND scov >= "
      f"{REF_MIN_COV:.0f}%   (near-complete, so they carry columns)\n")
for g in sorted(jobs):
    inc = sorted(cand[g])
    print(f"  {g:<16} in: {len(inc):>2}   {' '.join(inc)}")
    if below[g]:
        print(f"  {'':<16} below the floor, handed to 05_ instead of a tree:")
        for gene, (cv, al, tg) in sorted(below[g].items()):
            print(f"  {'':<18} {gene:<12} qcov {cv:>3}%  aln {al:>4} aa  "
                  f"(recovered by {tg})")


# ======================================================================
# 5. collect, align, infer
# ======================================================================
def collect(group, ts):
    """Partition the database into a family set and a root set."""
    fam_q, root_q = {}, {}
    for t in ts:
        for loc in [t["anchor_locus"]] + [x for x in t["family_members"].split(";") if x]:
            k = ara_key(loc)
            if k:
                fam_q[k] = loc.upper()
            else:
                print(f"    WARNING  {loc} has no sequence in aratha.fa")
    for loc in root_spec(ts):
        k = ara_key(loc)
        if k:
            root_q[k] = loc.upper()
        else:
            print(f"    WARNING  root locus {loc} has no sequence in aratha.fa")

    qf = W / f"{group}_q.faa"
    with open(qf, "w") as fh:
        for k in list(fam_q) + list(root_q):
            fh.write(f">{k}\n{ALL[k]}\n")

    out = W / f"{group}_q_hits.tsv"
    run(["diamond", "blastp", "-q", qf, "-d", W / "all", "-o", out,
         "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qcovhsp",
         "scovhsp", "evalue", "bitscore",
         "--max-target-seqs", "600", "--evalue", "1e-5",
         "--threads", THREADS, "--quiet"])

    best = {}
    for ln in open(out):
        f = ln.rstrip("\n").split("\t")
        q, sid = f[0], f[1]
        qcov, scov, bits = float(f[4]), float(f[5]), float(f[7])
        if qcov < REF_MIN_COV or scov < REF_MIN_COV:
            continue
        if sid not in best or bits > best[sid][0]:
            best[sid] = (bits, q)

    fam, root = set(), set()
    for sid, (_, q) in best.items():
        (root if q in root_q else fam).add(sid)
    # declared queries are whatever they were declared as
    fam |= set(fam_q); fam -= set(root_q)
    root |= set(root_q); root -= set(fam_q)

    forced = set()
    for gene in cand[group]:
        k = f"colesc|{gene}"
        if k in ALL:
            fam.add(k); forced.add(k)
    return fam, root, forced


def tree(group, keep, tag):
    faa = W / f"{group}_{tag}.faa"
    with open(faa, "w") as fh:
        for k in sorted(keep):
            fh.write(f">{label(k)}\n{ALL[k]}\n")
    aln = W / f"{group}_{tag}.aln"
    run(["famsa", "-t", THREADS, faa, aln],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    c0 = aln_cols(aln)

    # the comparable trimming, used for the tree
    trim = W / f"{group}_{tag}.trim.aln"
    run(["trimal", "-in", aln, "-out", trim, "-gt", str(TRIM_GT)],
        stdout=subprocess.DEVNULL)
    c1 = aln_cols(trim)

    # the old heuristic, recorded only so the change is auditable
    auto = W / f"{group}_{tag}.auto.aln"
    run(["trimal", "-in", aln, "-out", auto, "-automated1"],
        stdout=subprocess.DEVNULL)
    c2 = aln_cols(auto)

    if c1 < 50:
        print(f"    alignment trimmed to {c1} columns, too few to infer a tree")
        return dict(n=len(keep), cols=c0, trimmed=c1, auto=c2,
                    retention=(c1 / c0 if c0 else 0), treefile=None)

    pre = W / f"{group}_{tag}_iq"
    run([IQ, "-s", trim, "-m", "MFP", "-B", "1000", "-alrt", "1000",
         "-T", THREADS, "--prefix", pre, "-quiet", "-redo"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return dict(n=len(keep), cols=c0, trimmed=c1, auto=c2,
                retention=(c1 / c0 if c0 else 0),
                treefile=Path(str(pre) + ".treefile"))


print()
print("=" * 100)
print("TREES")
print("=" * 100)
print(f"\n  trimming: trimal -gt {TRIM_GT:.2f}, a fixed threshold so retention is")
print("  comparable between alignments. The -automated1 column count is shown")
print("  beside it for comparison with the earlier runs, and is not used.\n")

results, meta = {}, {}
for group, ts in sorted(jobs.items()):
    print(f"\n{group}")
    rule()
    fam, root, forced = collect(group, ts)
    sp_fam = len({k.split('|')[0] for k in fam})
    sp_root = len({k.split('|')[0] for k in root})
    print(f"  family set  {len(fam):>4} sequences across {sp_fam:>2} species "
          f"({len(forced)} taro candidates forced in)")
    print(f"  root set    {len(root):>4} sequences across {sp_root:>2} species")

    a = None
    if root:
        a = tree(group, fam | root, "A")
        print(f"  tree A  {a['n']:>4} tips  {a['cols']:>5} cols -> "
              f"{a['trimmed']:>4} at -gt {TRIM_GT:.2f} ({a['retention']*100:>3.0f}%)"
              f"   [-automated1 would give {a['auto']}]")
        if a["retention"] < RETENTION_WARN:
            print(f"          retention below {RETENTION_WARN*100:.0f}%: the root group is")
            print("          divergent enough to cost the columns close relatives need.")
            print("          Read tree B for resolution and tree A only for the boundary.")
    else:
        print("  tree A  skipped: no root group declared by any target in this job")

    b = tree(group, fam, "B")
    print(f"  tree B  {b['n']:>4} tips  {b['cols']:>5} cols -> "
          f"{b['trimmed']:>4} at -gt {TRIM_GT:.2f} ({b['retention']*100:>3.0f}%)"
          f"   [-automated1 would give {b['auto']}]")
    if a and b["n"] >= a["n"]:
        print("          WARNING tree B is not smaller than tree A. The root set")
        print("          did not separate. Do not read a boundary off tree A.")

    results[group] = dict(A=a, B=b)
    meta[group] = dict(
        targets=[t["target"] for t in ts],
        anchors={t["target"]: t["anchor_locus"] for t in ts},
        family_members={t["target"]: [x for x in t["family_members"].split(";") if x]
                        for t in ts},
        root_labels=sorted(label(k) for k in root),
        root_loci=root_spec(ts),
        taro=sorted(label(k) for k in forced),
        treeA=str(a["treefile"]) if a and a["treefile"] else "",
        treeB=str(b["treefile"]) if b["treefile"] else "",
        retentionA=(a["retention"] if a else None),
        retentionB=b["retention"],
        support_thresholds=dict(alrt=ALRT_MIN, ufboot=BOOT_MIN))

(W / "jobs.json").write_text(json.dumps(meta, indent=2))

print(f"""
{'=' * 100}

written: {W / 'jobs.json'}
         per job: the anchors, the root group's tip labels, the taro tips whose
         placement is the question, and the support thresholds

next:    python scripts/04b_read_trees.py

  The trees exist; nothing has been read off them. 04b_ does the reading:
  whether the root group is monophyletic in tree A, which anchor-defined clade
  each taro tip falls in, and whether the supporting node clears SH-aLRT >= {ALRT_MIN:.0f}
  and UFBoot >= {BOOT_MIN:.0f}. Keeping the reading in a separate script means a tree that
  took minutes to infer is not rebuilt whenever the reading changes, and the
  reading cannot reach back and alter how the tree was made.
""")
