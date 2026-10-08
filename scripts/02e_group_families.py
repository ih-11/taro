#!/usr/bin/env python3
"""
02e_group_families.py — de-duplicate the target table and record which targets
share a gene family.

Two problems in pathway_anchors.tsv that only became visible once 03_ ran.

1. BCH IS SEEDED TWICE.

   The hand list carried both BCH1 (AT4G25700) and BCH2 (AT5G52570). They are
   the two Arabidopsis paralogs of ONE family, the non-heme di-iron beta-carotene
   hydroxylases, so 02_ produced two targets, "BCH family" and "BCH family 2",
   for one family.

   The consequence in the 03_ output is a phantom result. Ces06731 is recovered
   by both, but its reverse best hit is BCH1, so it counts as a reciprocal best
   hit for "BCH family" and not for "BCH family 2". The grid therefore records
   "BCH family 2" as a stable zero across all nine threshold combinations, which
   would read in the master table as a family absent from taro. It is not
   absent. It is the same gene, counted once and missed once.

   Why this slipped through: 02c_ collapsed multi-member families only for the
   KEGG-only loci. NCED2/5/6/9, GGPS2/3/4/6, DXPS1/3 and IPP2 were all demoted
   to family members by that pass. BCH2 was never seen by it, because BCH2 was
   already in the hand list of 28 and 02c_ only triaged what KEGG added.

   This is a de-duplication, not a scope change. The set of Arabidopsis loci in
   scope is identical before and after. Recorded here, and in LOGBOOK.md, as
   having been found after 03_ ran.

2. SOME TARGETS DRAW FROM A SHARED FAMILY.

   CCD1, CCD4 and the NCEDs are all carotenoid cleavage dioxygenases. In the
   03_ output they recover from one pool of ten taro genes, which reciprocal
   best hit then partitions: one to CCD1, three to CCD4, six to NCED.

   That partition is a result, and a tree built from only one target's
   candidates has already assumed it. Asking whether taro has three CCD4 copies
   requires the NCEDs in the same tree, otherwise the tree cannot distinguish a
   CCD4 paralog from a divergent NCED. So the CCD targets get one tree, not
   three.

   PDS/ZDS/CRTISO, LCYB/LCYE and CYP97A/CYP97C are also same-family pairs or
   triples. All of them are detection targets, so no tree is built and no
   copy-number claim is made. The grouping is recorded anyway, because it says
   their candidate counts are not independent of each other.

   GGPPS and PSY are in no group. They keep the tree design already written.

This writes two new columns and changes no claim kind.
"""
import csv
import os
import sys
from pathlib import Path

CODE = Path(os.environ["TARO_CODE"])
ANCH = CODE / "results" / "tables" / "pathway_anchors.tsv"
REVIEW = CODE / "results" / "tables" / "kegg_scope_review.tsv"

# target name -> the Arabidopsis loci that are members of the same family as
# its anchor. These are the loci 02c_ and 02d_ demoted to family-member, plus
# BCH2, which 02c_ never saw.
FAMILY_MEMBERS = {
    "NCED family":  ["AT1G30100", "AT1G78390", "AT3G24220", "AT4G18350"],
    "GGPPS family": ["AT1G49530", "AT2G18640", "AT2G23800", "AT3G14550"],
    "DXS family":   ["AT3G21500", "AT5G11380"],
    "IDI family":   ["AT3G02780"],
    "BCH family":   ["AT5G52570"],
}

# targets that share one gene family. A group with more than one copy-number
# target gets ONE tree containing every member's candidates.
GROUPS = {
    "CCD1 family":   "CCD",
    "CCD4 family":   "CCD",
    "NCED family":   "CCD",
    "PDS family":    "CRT_desaturase",
    "ZDS family":    "CRT_desaturase",
    "CRTISO family": "CRT_desaturase",
    "LCYB family":   "LCY",
    "LCYE family":   "LCY",
    "CYP97A family": "CYP97",
    "CYP97C family": "CYP97",
}

if not ANCH.exists():
    sys.exit(f"{ANCH} not found")

rows = list(csv.DictReader(open(ANCH), delimiter="\t"))
n_before = len(rows)

# ---------------------------------------------------------------- 1. collapse
dropped = [r for r in rows if r["target"] == "BCH family 2"]
rows = [r for r in rows if r["target"] != "BCH family 2"]

print("=" * 78)
print("DE-DUPLICATION")
print("=" * 78)
if dropped:
    d = dropped[0]
    print(f"\n  dropped target : {d['target']}")
    print(f"  its anchor     : {d['anchor_gene']}  {d['anchor_locus']}")
    print(f"  now recorded as: a family member of BCH family (anchor AT4G25700)")
    print(f"\n  rows: {n_before} -> {len(rows)}")
else:
    print("\n  no 'BCH family 2' row; already collapsed or never created")

# ---------------------------------------------------------------- 2. annotate
for r in rows:
    t = r["target"]
    r["family_group"] = GROUPS.get(t, "")
    r["family_members"] = ";".join(FAMILY_MEMBERS.get(t, []))

print()
print("=" * 78)
print("FAMILY MEMBER SETS  —  what counts as inside the family for RBH")
print("=" * 78)
print("\n  A taro gene whose reverse best hit is any locus in this set is a")
print("  reciprocal best hit for the family, not only for the anchor. Without")
print("  this, a taro NCED whose closest Arabidopsis relative is NCED5 rather")
print("  than NCED3 is silently dropped.\n")
for r in rows:
    if r["family_members"]:
        print(f"  {r['target']:<16} {r['anchor_locus']}  + "
              f"{r['family_members'].replace(';', ' ')}")

# cross-check against the frozen scope file
if REVIEW.exists():
    frozen = {x["at_locus"].upper() for x in csv.DictReader(open(REVIEW), delimiter="\t")
              if x["decision"] == "family-member"}
    declared = set()
    for r in rows:
        declared |= {x.upper() for x in r["family_members"].split(";") if x}
    # BCH2 is expected to be in declared but not in frozen: it came from the
    # hand list, so 02c_ never triaged it.
    missing = frozen - declared
    extra = declared - frozen - {"AT5G52570"}
    print()
    if missing:
        print(f"  WARNING  frozen scope calls these family-member but no target "
              f"claims them: {' '.join(sorted(missing))}")
    if extra:
        print(f"  WARNING  claimed as family members but not family-member in the "
              f"frozen scope: {' '.join(sorted(extra))}")
    if not missing and not extra:
        print(f"  cross-check against kegg_scope_review.tsv: the {len(frozen)} "
              f"family-member loci all accounted for, plus BCH2 from the hand list")

print()
print("=" * 78)
print("TREE GROUPS")
print("=" * 78)
groups = {}
for r in rows:
    g = r["family_group"]
    if g:
        groups.setdefault(g, []).append((r["target"], r["claim_kind"]))

for g, members in sorted(groups.items()):
    cn = [t for t, k in members if k == "copy-number"]
    print(f"\n  {g}")
    for t, k in members:
        print(f"    {t:<16} {k}")
    if len(cn) > 1:
        print(f"    -> ONE tree for this group. {len(cn)} copy-number targets "
              f"cannot be rooted apart.")
    elif len(cn) == 1:
        print(f"    -> one tree, seeded by {cn[0]}, with the whole group's "
              f"candidates in it.")
    else:
        print(f"    -> no tree. All detection. Grouping recorded because the "
              f"counts are not independent.")

solo = [r["target"] for r in rows if not r["family_group"] and r["claim_kind"] == "copy-number"]
print(f"\n  ungrouped copy-number targets, tree design unchanged: {', '.join(solo)}")

FIELDS = ["pathway", "target", "anchor_gene", "anchor_locus", "target_level",
          "claim_kind", "outparalog", "in_kegg", "scope_reason",
          "family_group", "family_members"]

with open(ANCH, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=FIELDS, extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

print(f"""
written: {ANCH}   ({len(rows)} targets, 2 new columns)

  Nothing about which Arabidopsis loci are in scope changed. One duplicate
  target was removed and the family structure that was implicit in the 03_
  output is now written down.

next:  bash scripts/03_homology.sh 2>&1 | tee results/tables/log_03_homology.txt
""")
