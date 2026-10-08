#!/usr/bin/env python3
"""
06_readme_block.py — regenerate the README's findings section from the master
table, so prose cannot drift from data.

The README previously carried hand-typed counts. Over the course of this
analysis the CCD4 count moved from 3 to 2 and back to 3, the GGPPS count moved
from 1 to withheld, and the target count moved from 29 to 28. Every one of
those was a correction to the analysis, and every one would have left a stale
number in the README if the number were typed there by hand.

This writes the section between two markers in README.md:

    <!-- BEGIN part1-findings -->
    ... generated, do not edit by hand ...
    <!-- END part1-findings -->

If the markers are absent the block is printed to stdout instead, so you can
paste it where you want it the first time and run this afterwards.

Nothing here decides anything. Every statement is read out of
results/tables/part1_inventory.tsv and the support columns 04b_ wrote.
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

CODE = Path(os.environ["TARO_CODE"])
TBL = CODE / "results" / "tables"
README = CODE / "README.md"
BEGIN, END = "<!-- BEGIN part1-findings -->", "<!-- END part1-findings -->"

inv = list(csv.DictReader(open(TBL / "part1_inventory.tsv"), delimiter="\t"))
anchors = list(csv.DictReader(open(TBL / "pathway_anchors.tsv"), delimiter="\t"))
by_target = defaultdict(list)
for r in inv:
    by_target[r["target"]].append(r)

kind = {a["target"]: a["claim_kind"] for a in anchors}
pathway = {a["target"]: a["pathway"] for a in anchors}
L = []
w = L.append

n_cn = sum(1 for a in anchors if a["claim_kind"] == "copy-number")
n_det = sum(1 for a in anchors if a["claim_kind"] == "detection")
n_pre = sum(1 for a in anchors if a["claim_kind"] == "presence")

w("## Part 1 findings")
w("")
w("Generated from `results/tables/part1_inventory.tsv` by")
w("`scripts/06_readme_block.py`. Do not edit by hand.")
w("")
w(f"{len(anchors)} targets: {n_cn} copy-number, {n_det} detection, "
  f"{n_pre} presence. {len(inv)} target-and-gene rows.")
w("")
w("A copy number is reported only where a family assignment clearing "
  "SH-aLRT ≥ 80")
w("and UFBoot ≥ 95 meets an independent locus on one of the fourteen anchored")
w("sequences. A detection count is a lower bound, not a copy number.")
w("")

# ------------------------------------------------------------ copy-number
w("### Copy-number targets")
w("")
w("| target | copies verified | genes | basis |")
w("|---|---|---|---|")
detail = []
for a in anchors:
    if a["claim_kind"] != "copy-number":
        continue
    t = a["target"]
    rr = by_target.get(t, [])
    ver = [r for r in rr if r["evidence"] == "copy-number verified"]
    oth = [r for r in rr if r not in ver]
    if ver:
        basis = (f"tree A/B agree, UFBoot {ver[0]['tree_boot']}, "
                 f"{len({r['seq'] for r in ver})} anchored sequence(s)")
    else:
        basis = "see note below"
    w(f"| {t} | **{len(ver)}** | {' '.join(r['taro_gene'] for r in ver) or '—'} "
      f"| {basis} |")
    for r in oth:
        detail.append((t, r))
w("")
if detail:
    w("Genes in a copy-number family that do not contribute a verified copy:")
    w("")
    w("| target | gene | sequence | evidence |")
    w("|---|---|---|---|")
    for t, r in detail:
        w(f"| {t} | {r['taro_gene']} | {r['seq']} | {r['evidence']} |")
    w("")

# ------------------------------------------------------------ what the trees said
w("### Family membership, from the root-group split")
w("")
w("| job | root group | SH-aLRT | UFBoot | taro tips placed |")
w("|---|---|---|---|---|")
seen = set()
for r in inv:
    g = r["family_group"] or r["target"]
    if not r["member_alrt"] or g in seen:
        continue
    seen.add(g)
    tips = len({x["taro_gene"] for x in inv
                if (x["family_group"] or x["target"]) == g and x["member_alrt"]})
    w(f"| {g} | monophyletic | {r['member_alrt']} | {r['member_boot']} | {tips} |")
w("")

# ------------------------------------------------------------ detection
w("### Detection targets")
w("")
w("Lower bounds on family size. Reciprocal best hit cannot see a paralog")
w("whose closest Arabidopsis relative lies outside the declared family.")
w("")
w("| pathway | target | counted | independent loci | below the floor |")
w("|---|---|---|---|---|")
for a in anchors:
    if a["claim_kind"] != "detection":
        continue
    t = a["target"]
    rr = by_target.get(t, [])
    c = [r for r in rr if r["counted"] == "1"]
    s = [r for r in rr if r["counted"] != "1"]
    loci = len({r["taro_gene"] for r in c if r["independent_locus"] == "1"})
    w(f"| {a['pathway']} | {t} | {len(c)} | {loci} | {len(s) or '—'} |")
w("")

# ------------------------------------------------------------ presence
w("### Presence targets")
w("")
for a in anchors:
    if a["claim_kind"] != "presence":
        continue
    rr = by_target.get(a["target"], [])
    g = ", ".join(f"`{r['taro_gene']}` on {r['seq']}" for r in rr)
    w(f"- **{a['target']}**: {'present' if rr else 'ABSENT'}"
      + (f" — {g}" if g else ""))
w("")

# ------------------------------------------------------------ anomalies
split = [r for r in inv if "split model" in (r["notes"] or "")]
flag = [r for r in inv if r["length_flag"] == "1"]
if split or flag:
    w("### Annotation anomalies")
    w("")
if split:
    w("One gene across two records. Called only on non-overlapping collinear")
    w("anchor coverage plus same strand and no annotated gene between; genomic")
    w("distance is not used, because mean gene spacing here is 82 kb.")
    w("")
    w("| target | genes | sequence | anchor residues | combined span |")
    w("|---|---|---|---|---|")
    done = set()
    for r in split:
        partner = r["notes"].split("same locus as ")[-1].split(",")[0]
        key = tuple(sorted((r["taro_gene"], partner)))
        if key in done:
            continue
        done.add(key)
        o = next((x for x in inv if x["taro_gene"] == partner
                  and x["target"] == r["target"]), None)
        if not o:
            continue
        span = (max(int(r["end"]), int(o["end"]))
                - min(int(r["start"]), int(o["start"])) + 1)
        w(f"| {r['target']} | `{key[0]}` + `{key[1]}` | {r['seq']} "
          f"| {o['anchor_from']}-{o['anchor_to']} and "
          f"{r['anchor_from']}-{r['anchor_to']} | {span:,} bp |")
    w("")
if flag:
    w("One record spanning more than one gene, or a long lineage-specific")
    w("extension. Copy number withheld either way.")
    w("")
    w("| target | gene | protein | span | exons | sequence |")
    w("|---|---|---|---|---|---|")
    for r in flag:
        w(f"| {r['target']} | `{r['taro_gene']}` | {r['protein_aa']} aa "
          f"| {int(r['span_bp']):,} bp | {r['exons']} | {r['seq']} |")
    w("")

block = "\n".join(L)

if README.exists():
    txt = README.read_text()
    if BEGIN in txt and END in txt:
        pre = txt.split(BEGIN)[0]
        post = txt.split(END)[1]
        README.write_text(pre + BEGIN + "\n" + block + "\n" + END + post)
        print(f"updated {README} between the markers "
              f"({len(block.splitlines())} lines)")
        sys.exit(0)

print(f"{BEGIN} ... {END} not found in {README}.")
print("Paste the block below into README.md between those two markers, then")
print("re-run this script to keep it in step with the table.\n")
print(BEGIN)
print(block)
print(END)
