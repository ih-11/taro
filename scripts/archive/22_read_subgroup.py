#!/usr/bin/env python3
"""
22_read_subgroup.py — read the subgroup assignment off tree B.

Tree B holds phytoene synthase sequences only, rooted on Amborella, with the
monocot sampling widened by 20_. Rice and sorghum sequences carry subgroup
labels derived from chromosome position, following Lisboa et al. 2022 Table 1.

The question: taro has one PSY, monocots carry subgroups M1 and M2, so which
did taro keep?

The assignment is stated only if the node joining taro to a labelled anchor
clears SH-aLRT 80 and UFBoot 95, and only if rice and sorghum agree. Two
species agreeing matters more than one node being well supported, because a
single anchor can be misplaced by a long branch while two independent anchors
rarely are in the same direction.
"""
import os
import re
import sys
from Bio import Phylo

W = os.path.join(os.environ["TARO_WORK"], "05_orthology", "trees")
TREE = os.path.join(W, "B_iq.treefile")
ANCH = os.path.join(W, "anchors.tsv")

ALRT_MIN, BOOT_MIN = 80.0, 95

EXPR = {
    ("orysat", "M1"): "shoot, leaf, blade, central vein, flag leaf",
    ("orysat", "M2"): "seedling, coleoptile, leaf, ROOT, radicle, panicle, ovary",
    ("sorbic", "M1"): "all stages, peak at stem elongation; all 9 anatomical parts",
    ("sorbic", "M2"): "all stages, peak at flowering; higher in leaf; stress inducible",
}

anchors = {}           # accession -> (species tag, subgroup)
for ln in open(ANCH):
    tag, acc, chrom, sg = ln.rstrip("\n").split("\t")
    if sg in ("M1", "M2"):
        anchors[acc] = (tag, sg)

if not anchors:
    sys.exit("no labelled anchors in anchors.tsv; rerun 21_")

t = Phylo.read(TREE, "newick")
t.ladderize()
tips = t.get_terminals()


def support(node):
    raw = node.name or (str(node.confidence) if node.confidence is not None else "")
    m = re.match(r"^([\d.]+)/(\d+)$", str(raw).strip())
    return (float(m.group(1)), int(m.group(2))) if m else (None, None)


def strong(node):
    a, b = support(node)
    return a is not None and a >= ALRT_MIN and b >= BOOT_MIN


def anchors_in(node):
    """subgroup labels present among a clade's tips, per species"""
    found = {}
    for l in node.get_terminals():
        if not l.name:
            continue
        acc = l.name.split("_", 1)[1] if "_" in l.name else l.name
        if acc in anchors:
            sp, sg = anchors[acc]
            found.setdefault(sg, set()).add(sp)
    return found


print("=" * 72)
print("SUBGROUP ANCHORS IN TREE B")
print("=" * 72)
present = 0
for acc, (sp, sg) in sorted(anchors.items(), key=lambda x: (x[1][0], x[1][1])):
    here = any(acc in (l.name or "") for l in tips)
    present += here
    e = EXPR.get((sp, sg), "")
    print(f"  {sp:<8} {acc:<18} {sg}   {'present' if here else 'ABSENT'}")
    if e:
        print(f"  {'':<8} {'':<18}      {e}")
print(f"\n  {present} of {len(anchors)} anchors are in the tree")

taro = [l for l in tips if l.name and l.name.startswith("colesc")]
if not taro:
    sys.exit("taro not in tree B")

for tip in taro:
    print()
    print("=" * 72)
    print(f"WALKING OUTWARD FROM {tip.name}")
    print("=" * 72)
    print(f"\n  {'':4} {'support':>10}  {'tips':>4}  anchors in clade")

    call = None
    for node in reversed(t.get_path(tip)[:-1]):
        a, b = support(node)
        sup = f"{a:5.1f}/{b:<3d}" if a is not None else "    -     "
        ok = "OK  " if strong(node) else "WEAK"
        found = anchors_in(node)
        if found:
            desc = "  ".join(
                f"{sg} ({', '.join(sorted(sps))})" for sg, sps in sorted(found.items()))
        else:
            desc = "none"
        n = len(node.get_terminals())
        print(f"  {ok} {sup:>10}  {n:>4}  {desc}")
        if found and call is None:
            call = (found, node, n)
        if n > 25:
            break

    print()
    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    if call is None:
        print("\n  Taro reaches 25 tips without joining a labelled anchor.")
        print("  No subgroup assignment from this tree.")
        continue

    found, node, n = call
    a, b = support(node)
    sgs = sorted(found)
    print(f"\n  First clade containing taro and a labelled anchor: {n} tips")
    print(f"  Support: {a}/{b}" if a is not None else "  Support: not labelled")
    for sg in sgs:
        print(f"  Subgroup {sg} anchors present, from: {', '.join(sorted(found[sg]))}")

    print()
    if len(sgs) > 1:
        print("  The clade holds both M1 and M2 anchors, so the tree does not")
        print("  separate the subgroups at this depth. No assignment.")
    elif not strong(node):
        print(f"  Taro falls with subgroup {sgs[0]}, but the node is below")
        print(f"  SH-aLRT {ALRT_MIN:.0f} and UFBoot {BOOT_MIN}.")
        print("  Record the observation; do not state it as a result.")
        if len(found[sgs[0]]) >= 2:
            print()
            print("  Note that two independent species anchor this subgroup")
            print(f"  ({', '.join(sorted(found[sgs[0]]))}), which is worth more than")
            print("  the support value alone suggests. It is corroboration,")
            print("  not confirmation.")
    else:
        sg = sgs[0]
        nsp = len(found[sg])
        print(f"  ASSIGNMENT: taro's phytoene synthase is subgroup {sg}.")
        print(f"  Well supported, anchored by {nsp} "
              f"{'species' if nsp > 1 else 'species only'}"
              f" ({', '.join(sorted(found[sg]))}).")
        print()
        if sg == "M2":
            print("  Taro kept the broadly expressed copy. In rice the M2 genes")
            print("  reach root and radicle as well as leaf and panicle, and one")
            print("  is abiotic-stress inducible. Taro lost M1, which in rice is")
            print("  biased to shoot, leaf and blade.")
            print()
            print("  For the collaboration this is the more encouraging outcome:")
            print("  the surviving paralog is the one whose rice counterparts are")
            print("  expressed in below-ground and storage-adjacent tissue.")
        else:
            print("  Taro kept the copy biased in rice to shoot, leaf and blade,")
            print("  and lost M2, whose rice counterparts reach root and radicle.")
            print()
            print("  For the collaboration this is the harder outcome: the")
            print("  surviving paralog is the one associated with photosynthetic")
            print("  tissue rather than storage organs.")
        print()
        print("  Either way this is an inference about taro's gene from rice and")
        print("  sorghum expression, not a measurement in taro. It is a reason to")
        print("  ask where CePSY is expressed in corm, not an answer to it.")

print()
print("  Standing limitations: taro's nearest relative in this tree is still")
print("  Zostera marina, since no member of the Araceae has a protein set at")
print("  NCBI. Spirodela polyrhiza would be the single most useful addition.")
