#!/usr/bin/env python3
"""
02d_finalise_scope.py — clear the remaining REVIEW rows and freeze the scope.

02c_ left 19 rows unresolved. Most were unresolved for a mechanical reason: the
rules matched on KEGG's description text, while for these loci the information
is in the gene symbol before the semicolon. "ABA2; NAD(P)-binding Rossmann-fold
superfamily protein" carries no word my ABA rule looked for, but the symbol says
what it is.

This adds symbol-based rules for those, and records an explicit decision for the
two rows that were genuine scope questions.

THE SCOPE BOUNDARY, stated once

  In:  carotenoid pathway enzymes, the plastidial MEP chain supplying them, and
       the accumulation regulators (ORANGE, ORANGE-like, fibrillin).
  Out: everything else, including pathways that compete for the same precursor.

GGR, geranylgeranyl reductase, is the test of that boundary. It reduces GGPP to
phytyl-PP for chlorophyll and tocopherol, so it is the main competing sink for
the pool PSY draws on, and in a biofortification context that is interesting.

It is excluded anyway. GGPP also feeds gibberellins and protein
geranylgeranylation; if a competing sink is in scope because it competes, those
are in too, and the inventory stops being a carotenoid inventory. Competition
for GGPP is a flux question for Part 2, not a gene to inventory in Part 1.

GPS1 is the same shape: geranyl diphosphate synthase makes C10, not the C20
GGPP that PSY condenses. If the GGPPS family search recovers its orthologs
anyway, tree A places them and the tree reports it.
"""
import csv
import os
import sys
from pathlib import Path

CODE = Path(os.environ["TARO_CODE"])
F = CODE / "results" / "tables" / "kegg_scope_review.tsv"

# gene symbol -> (decision, reason).  Symbols are read from the text before the
# first semicolon in the KEGG name.
BY_SYMBOL = {
    "ABA2":   ("exclude", "xanthoxin dehydrogenase; ABA biosynthesis, downstream of NCED"),
    "FLDH":   ("exclude", "ABA-related short-chain dehydrogenase; downstream of NCED"),
    "ERA1":   ("exclude", "protein farnesyltransferase beta subunit; protein prenylation"),
    "ICME-LIKE1": ("exclude", "prenylcysteine methylesterase-like; protein prenylation"),
    "ICME-LIKE2": ("exclude", "prenylcysteine methylesterase-like; protein prenylation"),
    "FOLK":   ("exclude", "farnesol kinase; isoprenoid salvage, not the MEP chain"),
    "GGR":    ("exclude", "geranylgeranyl reductase; competing GGPP sink for chlorophyll "
                          "and tocopherol. Excluded because GGPP also feeds gibberellins "
                          "and protein geranylgeranylation, and a scope that admits one "
                          "competing sink admits all. Noted for the Part 2 flux discussion"),
    "GPS1":   ("exclude", "geranyl diphosphate synthase; makes C10 GPP, not the C20 GGPP "
                          "PSY condenses. If the GGPPS family search recovers its "
                          "orthologs, tree A places them"),
}

# loci whose KEGG name gives neither a usable symbol nor a usable description
BY_LOCUS = {
    "AT1G26640": ("exclude", "amino acid kinase family; in ath00900 but not the MEP chain"),
    "AT1G31910": ("exclude", "GHMP kinase family; mevalonate or phosphomevalonate kinase, "
                             "cytosolic MVA"),
    "AT3G54250": ("exclude", "GHMP kinase family; cytosolic MVA"),
    "AT1G74470": ("exclude", "pyridine nucleotide-disulfide oxidoreductase; not a "
                             "carotenoid or MEP enzyme"),
    "AT4G36470": ("exclude", "SAM-dependent methyltransferase superfamily; protein "
                             "prenylation processing"),
}

# unnamed loci. Not targets under any reading: a homology search seeded by the
# GGPPS anchor recovers them if they are family members, and tree A places them.
UNNAMED = ("AT3G14510", "AT2G18620", "AT3G14530", "AT3G20160",
           "AT3G29430", "AT3G32040")
UNNAMED_REASON = ("unnamed polyprenyl or terpenoid synthase. Not added as a target: "
                  "if it is a GGPPS family member the GGPPS family search recovers it "
                  "and tree A places it; if it is not, it is out of scope")

if not F.exists():
    sys.exit(f"{F} not found; run 02c_propose_scope.py first")

rows = list(csv.DictReader(open(F), delimiter="\t"))
before = sum(1 for r in rows if r["decision"] == "REVIEW")

for r in rows:
    if r["decision"] != "REVIEW":
        continue
    loc, name = r["at_locus"], r["kegg_name"]
    sym = name.split(";")[0].strip() if ";" in name else ""

    if sym in BY_SYMBOL:
        r["decision"], r["reason"] = BY_SYMBOL[sym]
    elif loc in BY_LOCUS:
        r["decision"], r["reason"] = BY_LOCUS[loc]
    elif loc in UNNAMED:
        r["decision"], r["reason"] = "exclude", UNNAMED_REASON

after = sum(1 for r in rows if r["decision"] == "REVIEW")

counts = {}
for r in rows:
    counts[r["decision"]] = counts.get(r["decision"], 0) + 1

print(f"REVIEW rows: {before} -> {after}\n")
print("=" * 70)
print("FINAL SCOPE DISPOSITION")
print("=" * 70)
for d in sorted(counts):
    print(f"  {counts[d]:>3}  {d}")

if after:
    print("\n  still unresolved:")
    for r in rows:
        if r["decision"] == "REVIEW":
            print(f"    {r['at_locus']}  {r['kegg_name'][:56]}")
    print("\n  Resolve these by hand before committing.")

with open(F, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["at_locus", "kegg_name", "bucket",
                                   "decision", "reason"])
    w.writeheader()
    w.writerows(rows)

print(f"""
written: {F}

  Read it before committing. Every decision is a judgement I made and you can
  overturn, and the reason column says what each one rests on.

  The two that were genuine scope questions rather than mechanical gaps:

    GGR   excluded. The main competing sink for GGPP, so interesting for flux,
          but GGPP also feeds gibberellins and protein geranylgeranylation. A
          scope that admits one competing sink admits all of them.

    GPS1  excluded. Different product, C10 rather than C20.

  If you disagree with either, change it now. After 03_ has seen taro results a
  scope change cannot be told apart from a result-driven one.
""")
