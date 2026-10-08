#!/usr/bin/env python3
"""
02b_annotate_kegg.py — name the loci KEGG has that we do not, and check the
loci we have that KEGG does not.

The cross-check in 02_ found 70 KEGG loci absent from the hand list. Deciding
which belong by recognising locus identifiers would be guessing from memory,
which is the habit that produced the squalene synthase error and the clade-rule
error. So KEGG names them instead.

Two questions, not one:

  1. Which of the 70 KEGG-only loci belong in a carotenoid and MEP inventory?
     Most will not. ath00900 carries the cytosolic mevalonate pathway, which is
     not the plastidial MEP pathway, and ath00906 continues into abscisic acid
     catabolism downstream of NCED.

  2. Why are five of our anchors absent from KEGG? ORANGE, ORANGE-like and
     fibrillin are expected, because they act on carotenoid accumulation
     without catalysing a step. CCD1 and NSX are not: CCD1 is a carotenoid
     cleavage dioxygenase and NSX is neoxanthin synthase, and both should
     appear in ath00906. Either KEGG assigns them different loci, or our locus
     identifiers are wrong. The second would be an error in the hand list.

Nothing is added or removed automatically. This writes a table for a decision.
"""
import csv
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

CODE = Path(os.environ["TARO_CODE"])
TAB = CODE / "results" / "tables"

OURS_ONLY = {
    "AT1G67080": "NSX",
    "AT3G63520": "CCD1",
    "AT4G04020": "FBN1a",
    "AT5G06130": "OR-like",
    "AT5G61670": "OR",
}


def kegg_get(entries):
    """KEGG /list accepts up to 10 entries per call."""
    out = {}
    for i in range(0, len(entries), 10):
        chunk = entries[i:i + 10]
        url = "https://rest.kegg.jp/list/" + "+".join(f"ath:{e}" for e in chunk)
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                txt = r.read().decode()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            print(f"  lookup failed for {chunk[0]}...: {exc}", file=sys.stderr)
            continue
        for ln in txt.splitlines():
            parts = ln.split("\t")
            if len(parts) >= 2:
                out[parts[0].replace("ath:", "").upper()] = parts[-1]
        time.sleep(0.34)          # KEGG asks for no more than 3 requests/second
    return out


# ------------------------------------------------------------- KEGG-only
src = TAB / "kegg_only_loci.tsv"
loci = [r["at_locus"] for r in csv.DictReader(open(src), delimiter="\t")
        if r["at_locus"]]
print(f"naming {len(loci)} KEGG-only loci\n")

names = kegg_get(loci)
if not names:
    sys.exit("KEGG lookup returned nothing; the decision cannot be made on data")

# Group by what the description says, so the decision is made on KEGG's own
# words rather than on recognising an identifier.
BUCKETS = [
    ("mevalonate / cytosolic",
     ("hydroxymethylglutaryl", "hmg-coa", "mevalonate", "mevalonate kinase",
      "diphosphomevalonate", "farnesyl")),
    ("abscisic acid, downstream of NCED",
     ("abscisic", "aba ", "8'-hydroxylase", "xanthoxin", "aldehyde oxidase",
      "beta-glucosidase")),
    ("carotenoid cleavage / apocarotenoid",
     ("carotenoid cleavage", "ccd", "9-cis-epoxycarotenoid", "nced")),
    ("prenyl / polyprenyl diphosphate",
     ("prenyl", "polyprenyl", "solanesyl", "geranyl", "dehydrodolichyl",
      "isopentenyl", "undecaprenyl")),
    ("carotenoid biosynthesis proper",
     ("phytoene", "carotene", "lycopene", "zeaxanthin", "violaxanthin",
      "neoxanthin", "carotenoid")),
]

grouped = {b[0]: [] for b in BUCKETS}
grouped["unclassified"] = []
for loc in sorted(names):
    d = names[loc].lower()
    for label, keys in BUCKETS:
        if any(k in d for k in keys):
            grouped[label].append((loc, names[loc]))
            break
    else:
        grouped["unclassified"].append((loc, names[loc]))

print("=" * 78)
print("KEGG LOCI ABSENT FROM THE HAND LIST, grouped by KEGG's own description")
print("=" * 78)

SUGGEST = {
    "mevalonate / cytosolic": "exclude — cytosolic MVA, not the plastidial MEP pathway",
    "abscisic acid, downstream of NCED": "exclude — downstream of the inventory's scope",
    "carotenoid cleavage / apocarotenoid": "REVIEW — same family as CCD4 and NCED",
    "prenyl / polyprenyl diphosphate": "REVIEW — only those feeding GGPP belong",
    "carotenoid biosynthesis proper": "REVIEW — likely belongs",
    "unclassified": "REVIEW — read the description",
}

rows = []
for label in [b[0] for b in BUCKETS] + ["unclassified"]:
    items = grouped[label]
    if not items:
        continue
    print(f"\n{label}  ({len(items)})")
    print(f"  suggested: {SUGGEST[label]}")
    for loc, desc in items:
        print(f"    {loc:<12} {desc[:78]}")
        rows.append(dict(at_locus=loc, kegg_name=desc, bucket=label,
                         suggestion=SUGGEST[label], decision="", reason=""))

# ------------------------------------------------------------- ours-only
print()
print("=" * 78)
print("OUR ANCHORS THAT KEGG DID NOT RETURN")
print("=" * 78)
print("\nIf KEGG knows the locus but places it outside these two pathways, the")
print("anchor is fine. If KEGG does not know the locus at all, the hand list")
print("has a wrong identifier, which is an error rather than a difference.\n")

ours = kegg_get(list(OURS_ONLY))
for loc, gene in sorted(OURS_ONLY.items(), key=lambda x: x[1]):
    d = ours.get(loc)
    if d:
        print(f"  {gene:<9} {loc}   KEGG knows it: {d[:60]}")
        print(f"  {'':<9} {'':<12}  -> outside ath00906/ath00900, anchor stands")
    else:
        print(f"  {gene:<9} {loc}   KEGG DOES NOT RECOGNISE THIS LOCUS")
        print(f"  {'':<9} {'':<12}  -> check the identifier in the hand list")

out = TAB / "kegg_only_loci.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["at_locus", "kegg_name", "bucket",
                                   "suggestion", "decision", "reason"])
    w.writeheader()
    w.writerows(rows)

print(f"""
written: {out}

  Fill in `decision` (add / exclude) and `reason` for anything in a REVIEW
  bucket. Loci marked exclude need no further work; loci added get a claim kind
  assigned in 02_ before 03_ re-runs, the same as every other anchor.
""")
