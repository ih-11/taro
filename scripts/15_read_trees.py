#!/usr/bin/env python3
"""
15_read_trees.py — extract the answers from the PSY and CCD gene trees.

Two questions:

  PSY  Taro has one full-length phytoene synthase. Cassava has the PSY1/PSY2
       pair that motivates the engineering plan. Where does taro's single copy
       sit relative to that pair, and did taro simply never undergo the
       duplication?

  CCD  Three taro genes were called CCD4 orthologs by reciprocal best hit, but
       two of them also matched NCED3. CCD and NCED are the same superfamily.
       The tree separates the clades, so we can see which taro genes are
       genuinely CCD4.
"""
import os, sys
from Bio import Phylo

W = os.environ["TARO_WORK"] + "/05_orthology/trees"
CODE = os.environ["TARO_CODE"]

SPNAME = {
    "colesc": "Colocasia esculenta (taro)",
    "manesc": "Manihot esculenta (cassava)",
    "aratha": "Arabidopsis thaliana",
    "orysat": "Oryza sativa",
    "zeamay": "Zea mays",
    "musacu": "Musa acuminata",
    "soltub": "Solanum tuberosum",
    "nelnuc": "Nelumbo nucifera",
    "ambtri": "Amborella trichopoda",
    "zosmar": "Zostera marina",
}

def load(name):
    t = Phylo.read(f"{W}/{name}.tree", "newick")
    t.root_at_midpoint()
    return t

def sp(leaf):
    return leaf.name.split("_")[0] if leaf.name else "?"

def neighbours(tree, target, k=6):
    """Leaves nearest to target by patristic distance."""
    leaves = tree.get_terminals()
    tgt = next((l for l in leaves if target in (l.name or "")), None)
    if not tgt:
        return None, []
    d = []
    for l in leaves:
        if l is tgt: continue
        d.append((tree.distance(tgt, l), l.name))
    d.sort()
    return tgt, d[:k]

# ---------------------------------------------------------------- PSY
print("=" * 72)
print("PHYTOENE SYNTHASE")
print("=" * 72)

psy = load("psy")
leaves = psy.get_terminals()
print(f"\n{len(leaves)} sequences in tree\n")

by_sp = {}
for l in leaves:
    by_sp.setdefault(sp(l), []).append(l.name.split("_", 1)[1])
for s in ["colesc", "manesc", "aratha", "orysat", "zeamay",
          "musacu", "soltub", "nelnuc", "ambtri", "zosmar"]:
    if s in by_sp:
        print(f"  {SPNAME[s]:<34} {len(by_sp[s])}  {', '.join(by_sp[s])}")

tgt, nb = neighbours(psy, "Ces24605")
if tgt:
    print(f"\nnearest relatives of taro Ces24605:\n")
    for dist, name in nb:
        s, g = name.split("_", 1)
        print(f"  {dist:7.3f}  {g:<18} {SPNAME.get(s, s)}")

print("\nnewick (paste into iTOL or ETE for a picture):")
print(open(f"{W}/psy.tree").read().strip()[:400] + " ...")

# ---------------------------------------------------------------- CCD
print("\n" + "=" * 72)
print("CAROTENOID CLEAVAGE DIOXYGENASE / NCED FAMILY")
print("=" * 72)

ccd = load("ccd")
print(f"\n{len(ccd.get_terminals())} sequences in tree")

# anchor on the two Arabidopsis references we know the identity of:
#   AT4G19170 = CCD4, AT3G14440 = NCED3
anchors = {}
for ln in open(os.environ["TARO_WORK"] + "/05_orthology/aratha_protein2gene.tsv"):
    p, g, sym, loc = ln.rstrip("\n").split("\t")
    if loc in ("AT4G19170", "AT3G14440"):
        anchors.setdefault(loc, []).append(p)

print("\ntaro genes ranked by distance to each Arabidopsis anchor:")
for loc, label in [("AT4G19170", "CCD4"), ("AT3G14440", "NCED3")]:
    acc = None
    for l in ccd.get_terminals():
        if l.name and any(a in l.name for a in anchors.get(loc, [])):
            acc = l; break
    if not acc:
        print(f"\n  {label} ({loc}) anchor not in tree")
        continue
    print(f"\n  distance to {label} ({loc}):")
    d = []
    for l in ccd.get_terminals():
        if l is acc or not l.name or not l.name.startswith("colesc"): continue
        d.append((ccd.distance(acc, l), l.name.split("_", 1)[1]))
    d.sort()
    for dist, g in d:
        print(f"    {dist:7.3f}  {g}")

out = f"{CODE}/results/tables/tree_summary.txt"
print(f"\n(save this output: {out})")
