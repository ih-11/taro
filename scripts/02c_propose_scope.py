#!/usr/bin/env python3
"""
02c_propose_scope.py — propose a disposition for each KEGG-only locus.

Reads  : results/tables/kegg_scope_review.tsv   (from 02b_, decision column empty)
Writes : results/tables/kegg_scope_review.tsv   (decision column proposed)

WHY PROPOSE RATHER THAN DECIDE

Sixty-nine blank rows invite rubber-stamping, and rubber-stamping a scope table
is worse than not having one. But a script cannot make a biological scope
decision either.

So every proposal here is derived from the name KEGG itself returned, and the
rule that produced it is printed beside it. Where the name does not settle the
question the proposal is REVIEW, with the specific question stated rather than
a guess. Those are the rows that need a person.

THE PRINCIPLE THAT DOES MOST OF THE WORK

A locus belonging to a family already targeted is not a new target. GGPS2, 3, 4
and 6 are members of the GGPPS family; NCED2, 5, 6 and 9 are members of the NCED
family; DXPS1 and DXPS3 are members of the DXS family. A homology search seeded
by any member recovers the rest, so these are recorded as family members rather
than added as targets. Adding them would double-count the same search.

Decisions used:
  exclude        different pathway; not part of this inventory
  family-member  recovered by an existing target's family search
  outparalog     excluded from scope, used as a phylogenetic reference
  REVIEW         the name does not settle it; a person decides
"""
import csv
import os
import sys
from pathlib import Path

CODE = Path(os.environ["TARO_CODE"])
F = CODE / "results" / "tables" / "kegg_scope_review.tsv"

# (substring in the KEGG name, decision, reason)
# Matched in order; first match wins. Substrings come from the names KEGG
# returned, not from recollection of what a locus does.
RULES = [
    # --- members of families we already target --------------------------
    ("geranylgeranyl pyrophosphate synthase", "family-member",
     "member of the GGPPS family; recovered by the GGPPS family search"),
    ("epoxycarotenoid dioxygenase", "family-member",
     "member of the NCED family; recovered by the NCED family search"),
    ("1-deoxy-d-xylulose 5-phosphate synthase", "family-member",
     "member of the DXS family; recovered by the DXS family search"),
    ("isopentenyl diphosphate isomerase", "family-member",
     "member of the IDI family; recovered by the IDI family search"),
    ("isopentenyl pyrophosphate:dimethylallyl", "family-member",
     "member of the IDI family; recovered by the IDI family search"),

    # --- excluded from scope, retained as phylogenetic reference ---------
    ("carotenoid cleavage dioxygenase 7", "outparalog",
     "strigolactone branch; excluded from scope, used as the CCD4 and NCED outparalog"),
    ("carotenoid cleavage dioxygenase 8", "outparalog",
     "strigolactone branch; excluded from scope, used as the CCD4 and NCED outparalog"),
    ("farnesyl diphosphate synthase", "outparalog",
     "cytosolic MVA; excluded from scope, proposed as the GGPPS outparalog"),

    # --- strigolactone branch -------------------------------------------
    ("beta-carotene isomerase d27", "exclude",
     "first committed step of strigolactone biosynthesis, not carotenoid accumulation"),
    ("cytochrome p450, family 711", "exclude",
     "MAX1; strigolactone biosynthesis"),
    ("2-oxoglutarate (2og) and fe(ii)-dependent oxygenase", "exclude",
     "LBO1; strigolactone biosynthesis"),

    # --- cytosolic mevalonate -------------------------------------------
    ("mevalonate", "exclude", "cytosolic MVA pathway, not the plastidial MEP pathway"),
    ("hydroxymethylglutaryl", "exclude", "cytosolic MVA pathway"),
    ("hydroxy methylglutaryl", "exclude", "cytosolic MVA pathway"),
    ("3-hydroxy-3-methylglutaryl", "exclude", "cytosolic MVA pathway"),
    ("thiolase", "exclude", "cytosolic MVA pathway"),
    ("acetoacetyl-coa", "exclude", "cytosolic MVA pathway"),

    # --- abscisic acid, downstream of NCED -------------------------------
    ("cytochrome p450, family 707", "exclude",
     "ABA 8'-hydroxylase; ABA catabolism, downstream of the inventory"),
    ("abscisic aldehyde oxidase", "exclude", "ABA biosynthesis, downstream of NCED"),
    ("beta glucosidase", "exclude", "ABA-glucose ester hydrolysis, downstream of NCED"),

    # --- protein prenylation and processing ------------------------------
    ("farnesyltransferase", "exclude", "protein prenylation, not isoprenoid biosynthesis"),
    ("farnesylated protein-converting", "exclude", "protein prenylation"),
    ("farnesylcysteine", "exclude", "protein prenylation"),
    ("isoprenylcysteine", "exclude", "protein prenylation"),
    ("prenylcysteine", "exclude", "protein prenylation"),
    ("peptidase family m48", "exclude", "CAAX protease, protein prenylation"),

    # --- long-chain polyprenyl, other products ---------------------------
    ("undecaprenyl pyrophosphate synthetase", "exclude",
     "cis-prenyltransferase; dolichol synthesis, not GGPP"),
    ("cis-prenyltransferase", "exclude", "dolichol synthesis, not GGPP"),
    ("solanesyl diphosphate synthase", "exclude",
     "plastoquinone side chain, not GGPP"),
]

# Rows the name does not settle. The question is stated rather than answered.
REVIEW = {
    "AT4G38460": ("GGR; geranylgeranyl reductase",
                  "Reduces GGPP to phytyl-PP for chlorophyll and tocopherol, so it "
                  "is the main competing sink for the GGPP pool. Relevant to flux "
                  "in a biofortification context but not a carotenoid pathway "
                  "enzyme. Include as a target, or note for Part 2 discussion?"),
    "AT2G34630": ("GPS1; geranyl diphosphate synthase 1",
                  "Makes GPP (C10), not GGPP (C20), so a different product. Same "
                  "trans-prenyltransferase superfamily, so the GGPPS family search "
                  "may recover its orthologs anyway. Family member, or out of scope?"),
    "AT3G14510": ("Polyprenyl synthetase family protein",
                  "Unnamed. Could be a GGPS family member or a distinct "
                  "polyprenyl synthase. Tree A for GGPPS will place it."),
    "AT2G18620": ("Terpenoid synthases superfamily protein", "Unnamed; see AT3G14510."),
    "AT3G14530": ("Terpenoid synthases superfamily protein", "Unnamed; see AT3G14510."),
    "AT3G20160": ("Terpenoid synthases superfamily protein", "Unnamed; see AT3G14510."),
    "AT3G29430": ("Terpenoid synthases superfamily protein", "Unnamed; see AT3G14510."),
    "AT3G32040": ("Terpenoid synthases superfamily protein", "Unnamed; see AT3G14510."),
    "AT5G40280": ("ERA1; Prenyltransferase family protein",
                  "Beta subunit of protein farnesyltransferase, so protein "
                  "prenylation, but the name alone does not say so."),
}

if not F.exists():
    sys.exit(f"{F} not found; run 02b_annotate_kegg.py first")

rows = list(csv.DictReader(open(F), delimiter="\t"))
print(f"{len(rows)} loci\n")

counts, review_rows = {}, []
for r in rows:
    loc, name = r["at_locus"], r["kegg_name"]
    low = name.lower()

    if loc in REVIEW:
        r["decision"] = "REVIEW"
        r["reason"] = REVIEW[loc][1]
        review_rows.append((loc, name, REVIEW[loc][1]))
    else:
        for key, dec, why in RULES:
            if key in low:
                r["decision"], r["reason"] = dec, why
                break
        else:
            r["decision"] = "REVIEW"
            r["reason"] = "no rule matched this KEGG name; decide by hand"
            review_rows.append((loc, name, r["reason"]))
    counts[r["decision"]] = counts.get(r["decision"], 0) + 1

print("=" * 76)
print("PROPOSED")
print("=" * 76)
for d in ("exclude", "family-member", "outparalog", "REVIEW"):
    if d in counts:
        print(f"  {counts[d]:>3}  {d}")

for d, title in [("family-member", "FAMILY MEMBERS of targets already defined"),
                 ("outparalog", "EXCLUDED FROM SCOPE, USED AS PHYLOGENETIC REFERENCE")]:
    sel = [r for r in rows if r["decision"] == d]
    if not sel:
        continue
    print(f"\n{title}")
    print("-" * 76)
    for r in sel:
        print(f"  {r['at_locus']:<12} {r['kegg_name'][:60]}")

print()
print("=" * 76)
print("NEEDS A PERSON")
print("=" * 76)
print("\nThe KEGG name does not settle these. Each proposal is REVIEW and the")
print("question is stated rather than answered.\n")
for loc, name, q in review_rows:
    print(f"  {loc}  {name[:62]}")
    for line in [q[i:i + 68] for i in range(0, len(q), 68)]:
        print(f"      {line}")
    print()

with open(F, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["at_locus", "kegg_name", "bucket",
                                   "decision", "reason"])
    w.writeheader()
    w.writerows(rows)

print("=" * 76)
print(f"""
written: {F}

  Every row now carries a proposed decision and the reason behind it. Change
  what you disagree with, resolve the REVIEW rows, then commit the file. That
  commit is the frozen scope.

  Nothing in the REVIEW set blocks 03_: every one of them is either already
  covered by an existing family search or excluded. They need resolving before
  the scope is called final, not before the homology step runs.
""")
