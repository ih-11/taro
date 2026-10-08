#!/usr/bin/env python3
"""
02_pathway_set.py — define the pathway anchor set, and check it for completeness.

The anchor list in the exploratory pass was assembled by hand and never
verified. That matters because this project has already missed one published
phytoene synthase, and an unverified gene set is the obvious place for the same
failure to be sitting undetected.

This script cross-checks the hand list against KEGG:

    ath00906  carotenoid biosynthesis
    ath00900  terpenoid backbone biosynthesis

KEGG's Arabidopsis gene identifiers are locus tags in the same format as ours,
so the comparison is direct and needs no identifier mapping.

A first version of this script queried `ko00906`, KEGG's KO-based reference
pathway. That returns no organism genes. The organism-specific identifier
`ath00906` is required. Worse than the wrong identifier was what the script did
with the result: it reported "0 Arabidopsis loci", then "only in our list: 28",
which reads as KEGG confirming all 28 anchors. It confirmed nothing, because the
check had not run.

So this version fails loudly. A pathway cannot contain zero genes, so an empty
response means the query failed, and a completeness check that cannot be
performed exits non-zero rather than passing silently. A check that quietly
succeeds when it did not run is worse than no check, because it buys false
confidence.
"""
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
MAP = WORK / "05_orthology" / "aratha_protein2gene.tsv"

# ---------------------------------------------------------------- hand list
# One Arabidopsis locus per pathway step. Locus identifiers rather than gene
# symbols: symbols are ambiguous, change, and the NCBI symbol query failed on
# this dataset.
ANCHORS = [
    # step,       gene,      locus,       claim kind for this family
    ("MEP",       "DXS",     "AT4G15560", "copy-number"),
    ("MEP",       "DXR",     "AT5G62790", "detection"),
    ("MEP",       "MCT",     "AT2G02500", "detection"),
    ("MEP",       "CMK",     "AT2G26930", "detection"),
    ("MEP",       "MDS",     "AT1G63970", "detection"),
    ("MEP",       "HDS",     "AT5G60600", "detection"),
    ("MEP",       "HDR",     "AT4G34350", "copy-number"),
    ("backbone",  "GGPPS11", "AT4G36810", "detection"),
    ("core",      "PSY",     "AT5G17230", "copy-number"),
    ("core",      "PDS3",    "AT4G14210", "detection"),
    ("core",      "Z-ISO",   "AT1G10830", "detection"),
    ("core",      "ZDS",     "AT3G04870", "detection"),
    ("core",      "CRTISO",  "AT1G06820", "detection"),
    ("cyclase",   "LCYB",    "AT3G10230", "detection"),
    ("cyclase",   "LCYE",    "AT5G57030", "detection"),
    ("xanth",     "CYP97A3", "AT1G31800", "detection"),
    ("xanth",     "CYP97C1", "AT3G53130", "detection"),
    ("xanth",     "BCH1",    "AT4G25700", "detection"),
    ("xanth",     "BCH2",    "AT5G52570", "detection"),
    ("xanth",     "ZEP",     "AT5G67030", "detection"),
    ("xanth",     "VDE",     "AT1G08550", "detection"),
    ("xanth",     "NSX",     "AT1G67080", "detection"),
    ("degrade",   "CCD1",    "AT3G63520", "detection"),
    ("degrade",   "CCD4",    "AT4G19170", "copy-number"),
    ("degrade",   "NCED3",   "AT3G14440", "copy-number"),
    ("sink",      "OR",      "AT5G61670", "presence"),
    ("sink",      "OR-like", "AT5G06130", "presence"),
    ("sink",      "FBN1a",   "AT4G04020", "presence"),
]

# Outparalogs for tree A. Declared here, before any tree is built, so the
# boundary test is not chosen after seeing where it falls.
#
# DXS and HDR are deliberately absent. Picking an outparalog from memory for a
# transketolase-family or reductase-family gene is the kind of guess that
# produced an earlier error, so 04_ skips them and neither gets a copy-number
# claim until one is looked up and recorded here.
OUTPARALOG = {
    "PSY":   [("SQS1", "AT4G34640"), ("SQS2", "AT4G34650")],
    "CCD4":  [("CCD7", "AT2G44990"), ("CCD8", "AT4G32810")],
    "NCED3": [("CCD7", "AT2G44990"), ("CCD8", "AT4G32810")],
}

# Organism-specific pathway identifiers. `ko00906` is the KO-based reference
# pathway and returns no organism genes.
KEGG = {"ath00906": "carotenoid biosynthesis",
        "ath00900": "terpenoid backbone biosynthesis"}


def kegg_loci(pathway):
    """
    Arabidopsis loci KEGG assigns to a pathway.

    Returns (status, set). Status is 'ok', 'unreachable' or 'empty'.
    An empty response is treated as a failure, not as a pathway with no genes,
    because the latter does not exist.
    """
    url = f"https://rest.kegg.jp/link/ath/{pathway}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            txt = r.read().decode()
    except (urllib.error.URLError, OSError, TimeoutError) as exc:
        return "unreachable", {"_reason": str(exc)}

    out = set()
    for ln in txt.splitlines():
        parts = ln.split("\t")
        if len(parts) == 2 and parts[1].startswith("ath:"):
            out.add(parts[1][4:].upper())
    return ("ok", out) if out else ("empty", set())


print("=" * 72)
print("PATHWAY ANCHOR SET")
print("=" * 72)
print(f"\nhand-curated list: {len(ANCHORS)} anchors")

# every anchor must resolve to a protein present in the proteome
loci_present = set()
for ln in open(MAP):
    if ln.startswith("protein\t"):
        continue
    loc = ln.rstrip("\n").split("\t")[3]
    if loc:
        loci_present.add(loc.upper())

missing = [a for a in ANCHORS if a[2].upper() not in loci_present]
if missing:
    print("\n  anchors that do not resolve in the Arabidopsis proteome:")
    for step, gene, loc, kind in missing:
        print(f"    {gene:<10} {loc}")
    print("\n  these are errors in the hand list and must be corrected")
    sys.exit(1)
print("  all anchors resolve in the Arabidopsis proteome")

print()
print("=" * 72)
print("CROSS-CHECK AGAINST KEGG")
print("=" * 72)

ours = {a[2].upper() for a in ANCHORS}
by_locus = {a[2].upper(): a for a in ANCHORS}
kegg_all = set()
failed = []

for pid, name in KEGG.items():
    status, loci = kegg_loci(pid)
    if status == "ok":
        kegg_all |= loci
        print(f"\n  {pid}  {name}: {len(loci)} Arabidopsis loci")
    elif status == "unreachable":
        print(f"\n  {pid}  {name}: UNREACHABLE — {loci.get('_reason', '')}")
        failed.append(pid)
    else:
        print(f"\n  {pid}  {name}: query returned no loci")
        failed.append(pid)

if failed:
    print(f"""
========================================================================
CROSS-CHECK FAILED
========================================================================

  Could not retrieve: {', '.join(failed)}

  A KEGG pathway cannot contain zero genes, so an empty response means the
  query failed rather than that the pathway is empty.

  This exits non-zero rather than reporting the hand list as unchecked,
  because the earlier version printed "only in our list: 28" after a failed
  query, which reads as KEGG confirming all 28 anchors.

  Check the endpoint by hand:
    curl -sS 'https://rest.kegg.jp/link/ath/ath00906' | head

  The organism-specific identifier (ath00906) is required; the KO-based
  reference identifier (ko00906) returns no organism genes.
""")
    sys.exit(1)

shared = ours & kegg_all
only_ours = ours - kegg_all
only_kegg = kegg_all - ours

print(f"\n  KEGG loci across both pathways : {len(kegg_all)}")
print(f"  in both                        : {len(shared)}")
print(f"  only in our list               : {len(only_ours)}")
print(f"  only in KEGG                   : {len(only_kegg)}")

if shared:
    print("\n  confirmed by KEGG:")
    for loc in sorted(shared):
        step, gene, _, kind = by_locus[loc]
        print(f"    {gene:<10} {loc}   ({step})")

if only_ours:
    print("\n  ours but not KEGG — each needs a reason to be here:")
    for loc in sorted(only_ours):
        step, gene, _, kind = by_locus[loc]
        print(f"    {gene:<10} {loc}   ({step})")
    print("""
    ORANGE, ORANGE-like and fibrillin are expected in this list. They act on
    carotenoid accumulation without catalysing a pathway step, so KEGG does
    not place them in ko00906. Anything else here should be justified or
    dropped.""")

if only_kegg:
    print(f"\n  KEGG but not ours — {len(only_kegg)} loci. This is where a gene")
    print("  the hand list missed would be:")
    for loc in sorted(only_kegg):
        print(f"    {loc}")
    print("""
    Not all of these belong. ath00900 includes the cytosolic mevalonate
    pathway, which is not the plastidial MEP pathway, and ath00906 continues
    into abscisic acid catabolism downstream of NCED. Each locus is judged on
    that basis rather than added wholesale.

    Any locus added here gets a claim kind assigned before the search runs,
    the same as the rest.""")
else:
    print("\n  Nothing in KEGG is absent from the hand list.")

# ---------------------------------------------------------------- write
out_dir = CODE / "results" / "tables"
out_dir.mkdir(parents=True, exist_ok=True)
out = out_dir / "pathway_anchors.tsv"
with open(out, "w") as fh:
    fh.write("step\tgene\tat_locus\tclaim_kind\tin_kegg\toutparalogs\n")
    for step, gene, loc, kind in ANCHORS:
        op = ";".join(f"{n}:{l}" for n, l in OUTPARALOG.get(gene, []))
        fh.write(f"{step}\t{gene}\t{loc}\t{kind}\t"
                 f"{'yes' if loc.upper() in kegg_all else 'no'}\t{op}\n")

kegg_out = out_dir / "kegg_only_loci.tsv"
with open(kegg_out, "w") as fh:
    fh.write("at_locus\tdecision\treason\n")
    for loc in sorted(only_kegg):
        fh.write(f"{loc}\t\t\n")

print()
print("=" * 72)
print("CLAIM KINDS, fixed before any search")
print("=" * 72)
kinds = {}
for step, gene, loc, kind in ANCHORS:
    kinds.setdefault(kind, []).append(gene)
for kind in ("copy-number", "detection", "presence"):
    g = kinds.get(kind, [])
    print(f"\n  {kind}  ({len(g)})")
    print(f"    {', '.join(g)}")

print(f"""
  copy-number families get tree A, tree B and the independent-locus check.
  detection families are reported as candidates recovered by reciprocal best
  hit, which is a lower bound on family size and is labelled as such.
  presence families need only a reciprocal best hit with adequate coverage.

  Assigning these before the search is what keeps the evidence standard from
  being set by which findings turn out to be interesting.

  DXS and HDR are copy-number families with no outparalog declared, so 04_
  skips them and neither gets a copy-number claim until one is looked up.

written: {out}
written: {kegg_out}   (fill in decision and reason for any locus added)
next:    bash scripts/03_homology.sh
""")
