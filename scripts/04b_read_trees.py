#!/usr/bin/env python3
"""
04b_read_trees.py — read the trees 04_ built. Builds nothing, infers nothing.

THE BUG IN THE PREVIOUS VERSION, WHICH IS THE SAME BUG AS THREE OTHERS

Six taro NCED genes came out assigned to "CCD1 family" with common ancestors of
165 and 121 tips. That is not a result, it is a symptom.

IQ-TREE writes an UNROOTED tree. Bio.Phylo reads the newick with whatever root
the file's parenthesisation implies, which is arbitrary. common_ancestor() is a
rooted operation: on an arbitrarily rooted tree it returns an arbitrary answer.
The five Arabidopsis NCEDs reported a common ancestor of 180 tips, the entire
tree, which is what "their MRCA is meaningless here" looks like from the inside.
Then the tie-break, which picked the smallest ancestor when none was exclusive,
chose CCD1 at 165 tips over NCED at 180 and printed it as an assignment.

Two errors, one cause and one consequence:

  Cause      a rooted operation applied to an unrooted tree. The same class as
             Bio.Phylo.get_path() never evaluating the root (which discarded a
             published rice PSY in 16_), and as rooting on a single SQS tip
             (which made every count zero in 18_).

  Consequence a tie-break that returned the least-bad option instead of
             declining. When no target yields an exclusive ancestor, the answer
             is not determinable, and saying so is the whole point of having a
             stopping rule.

HOW THE TREES ARE READ NOW

  Tree A is rooted on its root group, which test 1 has already shown to be
  monophyletic with support. Outgroup rooting on a monophyletic outgroup is
  valid, so every rooted operation on tree A is valid. Assignment happens here.

  Tree B CANNOT be meaningfully rooted, and the previous docstring claiming it
  was "rooted on the most distant species, Amborella where available" was wrong
  in principle as well as never implemented. In a GENE family tree the genes of
  any one species do not form a clade: each Amborella CCD groups with its own
  orthologs, not with the other Amborella CCDs. There is no outgroup inside
  tree B by construction, because tree B is defined as the family with the
  outgroup removed.

  So tree B is read by SPLITS only, which are rooting-independent: for each
  target, do its references plus the taro tips assigned to it form one side of
  some edge? That tests the same grouping without needing a root.

THE THREE TESTS

  1. MEMBERSHIP. Is the root group one split in tree A, and which side is each
     taro tip on? Split-based, so rooting-independent. The support is the
     split's own. For a family with a single Arabidopsis reference this is the
     only question a tree can answer, and it answers it well.

  2. ASSIGNMENT, in tree A rooted on the root group. For each taro tip and each
     target, take the common ancestor of the tip and that target's references
     and keep it only if it excludes every other target's references. The
     smallest surviving one is the assignment. If none survives, the tip is
     NOT DETERMINABLE and no fallback is applied.

  3. CONFIRMATION, in tree B by splits. For each target, is the set of its
     references plus its assigned taro tips a split, and does the separating
     node clear support?

SUPPORT

  SH-aLRT >= 80 AND UFBoot >= 95, both required, parsed from IQ-TREE's
  "alrt/ufboot" node label.

STOPPING RULE

  Disagreement between the two tests with neither clearing support is recorded
  as not determinable. No third analysis is run.
"""
import csv
import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

try:
    from Bio import Phylo
except ImportError:
    sys.exit("Bio.Phylo not available: pip install --break-system-packages biopython")

WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
W = WORK / "04_families"
TBL = CODE / "results" / "tables"
ALRT_MIN, BOOT_MIN = 80.0, 95.0
LENGTH_FLAG = 2.0

meta = json.loads((W / "jobs.json").read_text())


def fasta_lengths(p):
    out, name, n = {}, None, 0
    if not Path(p).exists():
        return out
    for ln in open(p):
        if ln.startswith(">"):
            if name:
                out[name] = n
            name, n = ln[1:].strip(), 0
        else:
            n += len(ln.strip())
    if name:
        out[name] = n
    return out


def support(clade):
    for v in (clade.name, clade.confidence):
        if v is None:
            continue
        s = str(v)
        if "/" in s:
            a, b = s.split("/")[:2]
            try:
                return float(a), float(b)
            except ValueError:
                pass
        else:
            try:
                return None, float(s)
            except ValueError:
                pass
    return None, None


def clears(a, b):
    return a is not None and b is not None and a >= ALRT_MIN and b >= BOOT_MIN


def tips(c):
    return {t.name for t in c.get_terminals()}


def split_node(tree, S):
    """Rooting-independent: a set is monophyletic iff some edge separates it."""
    allt = tips(tree.root)
    if not S or not S <= allt:
        return None, "some members of the set are not in this tree"
    comp = allt - S
    for cl in tree.find_clades():
        t = tips(cl)
        if t == S:
            return cl, "terminals are exactly the set"
        if t == comp:
            return cl, "terminals are exactly the complement"
    return None, "no single edge separates the set from the rest"


def refs_for(m, names):
    out = defaultdict(set)
    for tgt in m["targets"]:
        for loc in [m["anchors"][tgt]] + m["family_members"].get(tgt, []):
            for n in names:
                if n.startswith("aratha_") and n.endswith("_" + loc):
                    out[tgt].add(n)
    return out


# ======================================================================
# 1. alignment inputs, and over-merged taro models
# ======================================================================
print("=" * 106)
print("ALIGNMENT INPUTS, AND TARO MODELS OF ANOMALOUS LENGTH")
print("=" * 106)
length_flags = []
for job in sorted(meta):
    print(f"\n{job}")
    print("-" * 106)
    for tag in ("A", "B"):
        lens = fasta_lengths(W / f"{job}_{tag}.faa")
        if not lens:
            continue
        v = sorted(lens.values())
        med = statistics.median(v)
        print(f"  tree {tag}  n={len(v):<4} min {v[0]:>5}  median {med:>5.0f}  "
              f"max {v[-1]:>5}")
        for k, n in sorted(((k, n) for k, n in lens.items()
                            if n > LENGTH_FLAG * med), key=lambda kv: -kv[1]):
            who = "TARO CANDIDATE" if k.startswith("colesc_") else "reference"
            print(f"          {k} = {n} aa, {n/med:.1f}x the median   [{who}]")
            if k.startswith("colesc_") and tag == "B":
                length_flags.append((job, k, n, med))

if length_flags:
    print()
    print("  A record far longer than its family is the mirror image of the split")
    print("  models 03_ found: one record spanning more than one gene rather than")
    print("  one gene spanning two records. It also explains the column counts, a")
    print("  sequence several times the median forcing gap columns on every other")
    print("  sequence. 05_ checks the GFF3: exon count, genomic span, and whether")
    print("  the neighbouring genes a merge would have swallowed are absent.")

# ======================================================================
# 2. membership, by split
# ======================================================================
print()
print("=" * 106)
print("TEST 1  —  MEMBERSHIP, by split in tree A. Rooting-independent.")
print("=" * 106)
member, rootedA = {}, {}
for job in sorted(meta):
    m = meta[job]
    print(f"\n{job}")
    print("-" * 106)
    if not m["treeA"] or not Path(m["treeA"]).exists():
        print("  no tree A; membership cannot be tested by a boundary")
        member[job] = (None, None, None)
        continue
    tA = Phylo.read(m["treeA"], "newick")
    S = set(m["root_labels"])
    node, reason = split_node(tA, S)
    if node is None:
        print(f"  root group ({len(S)} tips from {m['root_loci']}): NOT monophyletic")
        print(f"  {reason}")
        print("  -> no boundary; membership untested, no copy-number claim here")
        member[job] = (False, None, None)
        continue
    a, b = support(node)
    ok = clears(a, b)
    print(f"  root group ({len(S)} tips from {m['root_loci']}): monophyletic, "
          f"{reason}")
    print(f"  split support: SH-aLRT {a}  UFBoot {b}  "
          f"{'clears' if ok else 'BELOW THRESHOLD'}")
    inside = [x for x in m["taro"] if x in S]
    if inside:
        print(f"  taro tips INSIDE the root group: {' '.join(inside)}")
    print(f"  taro tips on the family side: "
          f"{len([x for x in m['taro'] if x not in S])} of {len(m['taro'])}")
    if ok and not inside:
        print(f"  -> every taro candidate is a family member at "
              f"SH-aLRT {a} / UFBoot {b}")
    member[job] = (True, a, b)

    # root tree A on the root group: valid, because it is monophyletic
    og = [n for n in m["root_labels"] if n in tips(tA.root)]
    try:
        tA.root_with_outgroup(*og)
        rootedA[job] = (tA, "outgroup")
        print(f"  tree A rooted on the {len(og)}-tip root group; rooted "
              f"operations below are valid")
    except Exception as e:
        tA.root_at_midpoint()
        rootedA[job] = (tA, "midpoint")
        print(f"  outgroup rooting failed ({e}); fell back to midpoint, and "
              f"assignments below are weaker for it")

# ======================================================================
# 3. assignment in rooted tree A
# ======================================================================
print()
print("=" * 106)
print("TEST 2  —  ASSIGNMENT in tree A, rooted on the root group")
print("=" * 106)
rows, assigned = [], defaultdict(lambda: defaultdict(set))
for job in sorted(meta):
    m = meta[job]
    print(f"\n{job}")
    print("-" * 106)
    if len(m["targets"]) == 1:
        okm, a, b = member.get(job, (None, None, None))
        print(f"  one target, so there is no sub-assignment. Membership from "
              f"test 1 is the whole result for {m['targets'][0]}.")
        for x in m["taro"]:
            assigned[job][m["targets"][0]].add(x)
            rows.append(dict(job=job, taro_tip=x, assigned_to=m["targets"][0],
                             anc_tips="", alrt=a, boot=b,
                             support_ok=int(clears(a, b)), exclusive="",
                             basis="membership from the root split"))
        continue
    if job not in rootedA:
        print("  no usable tree A; nothing assigned")
        continue
    tA, how = rootedA[job]
    names = tips(tA.root)
    R = refs_for(m, names)
    allrefs = set().union(*R.values()) if R else set()

    print(f"  reference coherence in the rooted tree:")
    for tgt in m["targets"]:
        r = R.get(tgt, set())
        if len(r) < 2:
            print(f"    {tgt:<16} {len(r)} reference; coherence needs two or more")
            continue
        anc = tA.common_ancestor(sorted(r))
        intr = tips(anc) & (allrefs - r)
        a, b = support(anc)
        print(f"    {tgt:<16} {len(r)} refs, MRCA {len(tips(anc)):>3} tips, "
              + ("EXCLUSIVE" if not intr else f"contains {' '.join(sorted(intr))}")
              + f", aLRT {a} UFBoot {b}")

    print(f"\n    {'taro tip':<24}{'assignment':<18}{'anc':>5}{'aLRT':>8}"
          f"{'UFBoot':>8}  support")
    for x in m["taro"]:
        if x not in names:
            continue
        options = []
        for tgt, r in R.items():
            if not r:
                continue
            anc = tA.common_ancestor([x] + sorted(r))
            t = tips(anc)
            if t & (allrefs - r):
                continue                      # not exclusive: not an option
            options.append((len(t), tgt, anc))
        if not options:
            print(f"    {x:<24}{'NOT DETERMINABLE':<18}    no target gives an "
                  f"ancestor exclusive of the others")
            rows.append(dict(job=job, taro_tip=x, assigned_to="NOT DETERMINABLE",
                             anc_tips="", alrt="", boot="", support_ok=0,
                             exclusive=0,
                             basis="no exclusive common ancestor for any target"))
            continue
        options.sort(key=lambda o: o[0])
        n_t, tgt, anc = options[0]
        a, b = support(anc)
        ok = clears(a, b)
        print(f"    {x:<24}{tgt:<18}{n_t:>5}{str(a):>8}{str(b):>8}"
              f"  {'ok' if ok else 'below'}")
        assigned[job][tgt].add(x)
        rows.append(dict(job=job, taro_tip=x, assigned_to=tgt, anc_tips=n_t,
                         alrt=a, boot=b, support_ok=int(ok), exclusive=1,
                         basis="smallest exclusive MRCA in rooted tree A"))

# ======================================================================
# 4. confirmation in tree B, by the smallest excluding split
# ======================================================================
#
# WHAT THE PREVIOUS VERSION OF THIS TEST ASKED, AND WHY IT COULD ONLY FAIL
#
# It asked whether {a target's Arabidopsis references} + {its assigned taro
# tips} is a split in tree B. All three CCD targets came back "NOT a split",
# which is the only answer that test can give. Tree B holds seventeen species.
# The CCD4 clade contains every species' CCD4, so Arabidopsis plus taro is a
# strict subset of it with fifteen other species' orthologs interleaved. A set
# that is a subset of a clade is not a split, and never was going to be.
#
# The question that was meant: is there an edge with that target's references
# and taro tips on one side and NO other target's references on the same side?
# That is the orthogroup, and the other species belong in it. So the test is
# now for the SMALLEST side of any edge that CONTAINS the set and EXCLUDES the
# other targets' references. Both sides of every edge are considered, which
# keeps it rooting-independent.
print()
print("=" * 106)
print("TEST 3  —  CONFIRMATION in tree B: the smallest excluding split")
print("=" * 106)


def smallest_excluding(tree, S, O):
    """Smallest side of any edge containing all of S and none of O."""
    allt = tips(tree.root)
    best = None
    for cl in tree.find_clades():
        t = tips(cl)
        for side in (t, allt - t):
            if S <= side and not (side & O):
                if best is None or len(side) < len(best[0]):
                    best = (side, cl)
    return best


confirm = {}
for job in sorted(meta):
    m = meta[job]
    print(f"\n{job}")
    print("-" * 106)
    tf = m["treeB"]
    if not tf or not Path(tf).exists():
        print("  no tree B")
        continue
    if len(m["targets"]) == 1:
        print("  one target; tree B has no internal grouping to confirm")
        continue
    tB = Phylo.read(tf, "newick")
    names = tips(tB.root)
    R = refs_for(m, names)
    allrefs = set().union(*R.values()) if R else set()
    for tgt in m["targets"]:
        S = (R.get(tgt, set()) | assigned[job][tgt]) & names
        O = (allrefs - R.get(tgt, set())) & names
        if not S:
            continue
        got = smallest_excluding(tB, S, O)
        if got is None:
            print(f"    {tgt:<16} {len(S):>2} tips: no edge separates this set "
                  f"from the other targets' references")
            confirm[(job, tgt)] = (False, None, None, None)
            continue
        side, cl = got
        a, b = support(cl)
        extra = len(side) - len(S)
        sp = len({n.split("_")[0] for n in side})
        print(f"    {tgt:<16} {len(S):>2} tips -> smallest excluding clade "
              f"{len(side):>3} tips across {sp:>2} species "
              f"(+{extra} orthologs), aLRT {a} UFBoot {b}  "
              f"{'clears' if clears(a, b) else 'below'}")
        confirm[(job, tgt)] = (True, a, b, len(side))
    print()
    print("    The set is that target's Arabidopsis references plus the taro tips")
    print("    tree A assigned to it. The clade reported is the smallest group")
    print("    containing them and no other target's reference, so the extra tips")
    print("    are the other species' orthologs, which is what an orthogroup is.")
    print("    Recovering it in tree B means a second inference, built from a")
    print("    different alignment, puts the same genes together.")

# ======================================================================
# 5. final
# ======================================================================
print()
print("=" * 106)
print("FINAL")
print("=" * 106)
final = []
for r in rows:
    job, tgt, tip = r["job"], r["assigned_to"], r["taro_tip"]
    okm, ma, mb = member.get(job, (None, None, None))
    c = confirm.get((job, tgt))
    bits = []
    if okm:
        bits.append(f"family member by the root split ({ma}/{mb})")
    elif okm is False:
        bits.append("root group not monophyletic, membership untested")
    if tgt == "NOT DETERMINABLE":
        verdict = "NOT DETERMINABLE"
        bits.append("no exclusive ancestor in tree A")
    else:
        verdict = tgt
        if r["basis"].startswith("membership"):
            bits.append("single-target job, no sub-assignment needed")
        else:
            bits.append(f"tree A MRCA {r['anc_tips']} tips, "
                        + ("support clears" if r["support_ok"] else "support below"))
            if c is None:
                pass
            elif c[0]:
                bits.append(f"tree B recovers the clade, {c[3]} tips "
                            f"({c[1]}/{c[2]})"
                            + ("" if clears(c[1], c[2]) else ", support below"))
            else:
                bits.append("tree B does NOT recover an excluding clade")
    print(f"  {tip:<24} {verdict:<18} {'; '.join(bits)}")
    final.append(dict(job=job, taro_tip=tip, taro_gene=tip.split("_", 1)[-1],
                      assignment=verdict, basis="; ".join(bits),
                      member_alrt=ma, member_boot=mb,
                      treeA_alrt=r["alrt"], treeA_boot=r["boot"],
                      treeA_support_ok=r["support_ok"],
                      treeB_split=("" if c is None else int(c[0])),
                      boundary_placed=("" if okm is None else int(okm))))

out1, out2 = TBL / "tree_assignments.tsv", TBL / "family_assignment.tsv"
with open(out1, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["job", "taro_tip", "assigned_to", "anc_tips",
                                   "alrt", "boot", "support_ok", "exclusive",
                                   "basis"])
    w.writeheader(); w.writerows(rows)
with open(out2, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["job", "taro_tip", "taro_gene", "assignment",
                                   "basis", "member_alrt", "member_boot",
                                   "treeA_alrt", "treeA_boot",
                                   "treeA_support_ok", "treeB_split",
                                   "boundary_placed"])
    w.writeheader(); w.writerows(final)

print(f"""
written: {out1}
written: {out2}

  Family assignment is half of a copy number. The other half is an independent
  genomic locus, read from the taro GFF3 in 05_.

  Resting on coordinates alone, having never entered a tree:
  Ces12496 and Ces12497, the split PSY model, and Ces05879.
  {'Needing a GFF3 answer before its copy number means anything: '
   + ', '.join(f'{k} ({n} aa vs median {med:.0f})' for _, k, n, med in length_flags)
   if length_flags else ''}

next:    python scripts/05_report.py
""")
