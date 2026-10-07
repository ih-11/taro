#!/usr/bin/env python3
"""
19_subgroup_call.py

Read the subgroup assignment off the PSY tree, and fix the boundary test that
18_ got wrong.

WHAT 18_ GOT WRONG
    Part 3 rooted the tree on a single squalene synthase sequence
    (-o SQS_aratha_NP_195190.1). That places one outgroup sequence on one side
    of the root and the other seventeen on the other, so no internal clade can
    contain all eighteen except the root. The test reported the whole tree,
    every count came out zero, and the five verdicts were an artifact of the
    rooting rather than evidence about the sequences.

    Fixed here by re-rooting at the midpoint and then asking whether the
    outgroup is monophyletic, reporting honestly if it is not.

WHAT 18_ GOT RIGHT, AND WHY IT MATTERS
    Part 4 mapped our three rice sequences to chromosomes 6, 12 and 9, which
    are exactly the three rice PSY genes in Lisboa et al. 2022 Table 1:

        chr6   LOC_Os06g51290   M2   seedling, leaf, panicle, floret, ovary
        chr12  LOC_Os12g43130   M1   shoot, leaf, blade, central vein, flag leaf
        chr9   LOC_Os09g38320   M2   seedling, coleoptile, leaf, ROOT, radicle,
                                     panicle, ovary; abiotic-stress inducible

    One of those three, XP_015611707.1 on chromosome 9, was among the five
    sequences 16_ excluded as wider superfamily. A published, functionally
    characterised phytoene synthase was therefore discarded by that rule, which
    confirms the rule was unsound and that the naive hit counts were correct.

THE QUESTION THIS SCRIPT ANSWERS
    Monocots carry subgroups M1 and M2 only. Taro has one PSY, so it has lost
    one of the two. Which one it kept is read from the rice sequence it groups
    with, and the support at that node decides whether the assignment can be
    stated or only noted.

    If taro kept M1 it kept the copy that in rice is biased to photosynthetic
    tissue, and lost the copy that reaches root and radicle. For a corm crop
    that bears on whether the storage organ can accumulate carotenoid.
"""
import os
import re
import sys
from Bio import Phylo

W = os.path.join(os.environ["TARO_WORK"], "05_orthology", "trees")
TREE = os.path.join(W, "psy_og_iq.treefile")

# Rice subgroup anchors, from Lisboa et al. 2022 Table 1, resolved to our
# protein accessions by chromosome in 18_.
RICE = {
    "XP_015618170.1": ("M1", "LOC_Os12g43130", "shoot, leaf, blade, central vein, flag leaf"),
    "XP_015611707.1": ("M2", "LOC_Os09g38320", "seedling, coleoptile, leaf, ROOT, radicle, panicle, ovary"),
    "XP_015642514.1": ("M2", "LOC_Os06g51290", "seedling, leaf, panicle, floret, ovary, shoot"),
}

NAMES = {
    "colesc": "Colocasia esculenta (taro)", "manesc": "Manihot esculenta (cassava)",
    "aratha": "Arabidopsis thaliana", "orysat": "Oryza sativa", "zeamay": "Zea mays",
    "musacu": "Musa acuminata", "soltub": "Solanum tuberosum",
    "nelnuc": "Nelumbo nucifera", "ambtri": "Amborella trichopoda",
    "zosmar": "Zostera marina",
}
ORDER = ["colesc", "zosmar", "musacu", "orysat", "zeamay",
         "manesc", "aratha", "soltub", "nelnuc", "ambtri"]

ALRT_MIN, BOOT_MIN = 80.0, 95


def support(node):
    """IQ-TREE writes internal labels as 'SH-aLRT/UFBoot'."""
    raw = node.name or (str(node.confidence) if node.confidence is not None else "")
    m = re.match(r"^([\d.]+)/(\d+)$", str(raw).strip())
    return (float(m.group(1)), int(m.group(2))) if m else (None, None)


def mark(node):
    a, b = support(node)
    if a is None:
        return "    ", "   -   "
    ok = "OK  " if (a >= ALRT_MIN and b >= BOOT_MIN) else "WEAK"
    return ok, f"{a:5.1f}/{b:<3d}"


def species_of(leaf):
    n = leaf.name or ""
    return n.split("_")[1] if n.startswith("SQS_") else n.split("_")[0]


def acc_of(leaf):
    n = leaf.name or ""
    return n[4:].split("_", 1)[1] if n.startswith("SQS_") else n.split("_", 1)[1]


t = Phylo.read(TREE, "newick")
t.root_at_midpoint()
t.ladderize()
tips = t.get_terminals()

sqs = [l for l in tips if l.name and l.name.startswith("SQS_")]
qry = [l for l in tips if l.name and not l.name.startswith("SQS_")]

print("=" * 70)
print("BOUNDARY TEST, corrected rooting")
print("=" * 70)
print(f"  {len(tips)} tips: {len(sqs)} squalene synthase, {len(qry)} query")

sqs_names = {l.name for l in sqs}
best = None
for nd in t.get_nonterminals():
    names = {l.name for l in nd.get_terminals()}
    if sqs_names <= names and (best is None or
                               len(names) < len(best.get_terminals())):
        best = nd

if best is None or len(best.get_terminals()) == len(tips):
    print("\n  The squalene synthase sequences are not monophyletic after")
    print("  midpoint rooting, so this tree cannot place the PSY boundary.")
    print("  That is a statement about the tree, not about the sequences:")
    print("  only 37% of alignment columns survived trimming, because PSY and")
    print("  squalene synthase are divergent enough that much of the alignment")
    print("  is not confidently homologous.")
    print("\n  The subgroup assignment below does not depend on this test,")
    print("  because it is read from rice sequences of published identity.")
    boundary_ok = False
else:
    inside = {l.name for l in best.get_terminals()}
    intruders = sorted(n for n in inside if n not in sqs_names)
    a, b = support(best)
    print(f"\n  smallest clade holding every outgroup sequence: "
          f"{len(inside)} tips, support {a}/{b}" if a else
          f"\n  smallest clade holding every outgroup sequence: {len(inside)} tips")
    print(f"  query sequences falling inside it: {len(intruders)}")
    for n in intruders:
        print(f"    {n}")
    boundary_ok = True

    counts = {}
    for l in qry:
        if l.name not in inside:
            counts.setdefault(species_of(l), []).append(acc_of(l))
    print("\n  PSY count per species (query sequences outside the outgroup clade):")
    for s in ORDER:
        g = counts.get(s, [])
        print(f"    {NAMES[s]:<34}{len(g):>3}  {', '.join(g) if g else '-'}")

# ----------------------------------------------------------------- subgroup
print()
print("=" * 70)
print("SUBGROUP ASSIGNMENT  —  which monocot copy did taro keep?")
print("=" * 70)

taro = [l for l in tips if l.name and l.name.startswith("colesc")]
if not taro:
    sys.exit("  taro not in this tree")

print("\n  rice anchors (Lisboa et al. 2022 Table 1, mapped by chromosome):")
for acc, (sg, loc, expr) in RICE.items():
    present = any(acc in (l.name or "") for l in tips)
    print(f"    {acc:<18} {sg:<3} {loc:<16} {'in tree' if present else 'ABSENT'}")
    print(f"    {'':<18} {expr}")

for tip in taro:
    print(f"\n  walking outward from {tip.name}:\n")
    print(f"  {'':4} {'support':>9}  {'tips':>4}  composition")
    path = t.get_path(tip)
    first_rice = None
    for node in reversed(path[:-1]):
        leaves = node.get_terminals()
        ok, sup = mark(node)
        rice_here = sorted({RICE[a][0] for l in leaves
                            for a in RICE if a in (l.name or "")})
        sp = sorted({species_of(l) for l in leaves})
        tag = f"   <- rice {'/'.join(rice_here)}" if rice_here else ""
        print(f"  {ok} {sup:>9}  {len(leaves):>4}  "
              f"{' '.join(sp)}{tag}")
        if rice_here and first_rice is None:
            first_rice = (rice_here, node, len(leaves))
        if len(leaves) > 20:
            break

    print()
    if first_rice is None:
        print("  Taro does not group with any labelled rice sequence before the")
        print("  clade grows past 20 tips. No subgroup assignment from this tree.")
    else:
        sgs, node, n = first_rice
        a, b = support(node)
        strong = a is not None and a >= ALRT_MIN and b >= BOOT_MIN
        print(f"  First clade containing taro and a labelled rice sequence:")
        print(f"    {n} tips, subgroup {'/'.join(sgs)}, support "
              f"{a}/{b}" if a is not None else f"    {n} tips, no support label")
        print()
        if len(sgs) == 1 and strong:
            sg = sgs[0]
            print(f"  ASSIGNMENT: taro's PSY is subgroup {sg}, well supported.")
            if sg == "M1":
                print("  Taro kept the copy that in rice is biased to shoot, leaf")
                print("  and blade, and lost M2, which in rice reaches root and")
                print("  radicle and is abiotic-stress inducible.")
            else:
                print("  Taro kept the broadly expressed copy, which in rice")
                print("  includes root and radicle, and lost M1.")
        elif len(sgs) == 1:
            print(f"  Taro falls with subgroup {sgs[0]}, but the node is below")
            print("  SH-aLRT 80 / UFBoot 95. Note it; do not state it.")
        else:
            print("  Taro's nearest clade contains both M1 and M2 rice genes, so")
            print("  the tree does not separate them at this depth. No assignment.")

print()
print("  Caveat for either outcome: 37% of alignment columns survived trimming")
print("  and the PSY backbone was already weakly supported in 17_. A subgroup")
print("  assignment from this tree is a working hypothesis, and confirming it")
print("  would need the Lisboa et al. reference sequences themselves rather")
print("  than rice standing in for them.")
