#!/usr/bin/env python3
"""
02b_annotate_kegg.py — name the KEGG-only loci so scope decisions are made on
descriptions rather than on recognising identifiers.

Reads  : results/tables/kegg_only_loci_raw.tsv   (written by 02_, never altered)
Writes : results/tables/kegg_scope_review.tsv    (for a human to fill in)

The split matters. An earlier version read and overwrote the same file, which
destroyed the raw list the moment the annotation ran and left no record of what
02_ actually produced.

THE BUCKETS ARE NAVIGATION, NOT EVIDENCE.

Keyword grouping cannot make a biological decision and is not meant to.
"isopentenyl" catches IDI, which belongs, alongside prenyltransferases that do
not. The carotenoid-cleavage bucket mixes CCD1, CCD4, NCED, CCD7 and CCD8,
whose roles in this inventory differ completely: CCD4 and NCED are targets,
CCD7 and CCD8 are excluded from scope but used as phylogenetic outparalogs.

So every bucket is a sorting aid for reading 70 descriptions quickly. Nothing is
added or excluded automatically, and the decision column starts empty.
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
RAW = TAB / "kegg_only_loci_raw.tsv"
OUT = TAB / "kegg_scope_review.tsv"

# anchors KEGG did not return; checked to distinguish "KEGG files it elsewhere"
# from "our locus identifier is wrong"
OURS_ONLY = {
    "AT1G67080": "NSX", "AT3G63520": "CCD1", "AT4G04020": "FBN1a",
    "AT5G06130": "OR-like", "AT5G61670": "OR",
}


def kegg_names(entries):
    out = {}
    for i in range(0, len(entries), 10):
        chunk = entries[i:i + 10]
        url = "https://rest.kegg.jp/list/" + "+".join(f"ath:{e}" for e in chunk)
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                txt = r.read().decode()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            print(f"  lookup failed near {chunk[0]}: {exc}", file=sys.stderr)
            continue
        for ln in txt.splitlines():
            p = ln.split("\t")
            if len(p) >= 2:
                out[p[0].replace("ath:", "").upper()] = p[-1]
        time.sleep(0.34)          # KEGG asks for at most 3 requests per second
    return out


if not RAW.exists():
    sys.exit(f"{RAW} not found; run 02_pathway_set.py first")

loci = [ln.strip() for ln in open(RAW) if ln.strip() and ln.strip() != "at_locus"]
print(f"naming {len(loci)} KEGG-only loci\n")
names = kegg_names(loci)
if not names:
    sys.exit("KEGG returned nothing; the review cannot be prepared from data")

BUCKETS = [
    ("mevalonate / cytosolic",
     ("hydroxymethylglutaryl", "hmg-coa", "mevalonate", "farnesyl", "thiolase",
      "acetoacetyl")),
    ("abscisic acid, downstream of NCED",
     ("abscisic", "707a", "xanthoxin", "aldehyde oxidase", "glucosidase")),
    ("carotenoid cleavage / apocarotenoid",
     ("carotenoid cleavage", "epoxycarotenoid", "carotene isomerase")),
    ("protein prenylation",
     ("prenyltransferase a", "farnesyltransferase", "isoprenylcysteine",
      "prenylcysteine", "peptidase family m48", "methyltransferase")),
    ("prenyl / polyprenyl diphosphate",
     ("prenyl", "polyprenyl", "solanesyl", "geranyl", "undecaprenyl",
      "isopentenyl")),
    ("carotenoid biosynthesis proper",
     ("phytoene", "lycopene", "zeaxanthin", "violaxanthin", "neoxanthin",
      "carotenoid")),
]
NOTE = {
    "mevalonate / cytosolic": "likely exclude: cytosolic MVA, not plastidial MEP",
    "abscisic acid, downstream of NCED": "likely exclude: downstream of scope",
    "carotenoid cleavage / apocarotenoid":
        "READ EACH: CCD4 and NCED are targets; CCD7, CCD8 and D27 are the "
        "strigolactone branch, excluded from scope but used as outparalogs",
    "protein prenylation": "likely exclude: protein modification",
    "prenyl / polyprenyl diphosphate":
        "READ EACH: GGPS members belong, IDI belongs, dolichol and solanesyl "
        "chains do not",
    "carotenoid biosynthesis proper": "READ EACH: likely belongs",
    "unclassified": "READ EACH",
}

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
print("KEGG LOCI ABSENT FROM THE SCOPE")
print("grouped by KEGG's own description; grouping is navigation, not evidence")
print("=" * 78)

rows = []
for label in [b[0] for b in BUCKETS] + ["unclassified"]:
    items = grouped[label]
    if not items:
        continue
    print(f"\n{label}  ({len(items)})\n  {NOTE[label]}")
    for loc, desc in items:
        print(f"    {loc:<12} {desc[:76]}")
        rows.append(dict(at_locus=loc, kegg_name=desc, bucket=label,
                         decision="", reason=""))

print()
print("=" * 78)
print("ANCHORS KEGG DID NOT RETURN")
print("=" * 78)
print("\nKEGG knowing the locus but filing it elsewhere means the anchor is")
print("fine. KEGG not knowing it at all means our identifier is wrong.\n")
ours = kegg_names(list(OURS_ONLY))
for loc, gene in sorted(OURS_ONLY.items(), key=lambda x: x[1]):
    d = ours.get(loc)
    if d:
        print(f"  {gene:<9} {loc}  ok, filed elsewhere: {d[:52]}")
    else:
        print(f"  {gene:<9} {loc}  NOT RECOGNISED — check the identifier")

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t",
                       fieldnames=["at_locus", "kegg_name", "bucket",
                                   "decision", "reason"])
    w.writeheader()
    w.writerows(rows)

print(f"""
written: {OUT}

  Fill in decision (add / exclude) and reason for every row, then commit it.
  That file is the frozen record of what the scope includes and why.

  Anything added gets a target, an anchor, a target_level and a claim_kind in
  02_ before 03_ runs. Changing scope after 03_ has seen taro results cannot be
  distinguished from a result-driven decision, which is the reason this is
  settled first.
""")
