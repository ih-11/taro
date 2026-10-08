#!/usr/bin/env bash
# 03_preflight.sh — check that 03_ and 04_ speak the current anchors schema,
# and patch them if they do not.
#
# Why this exists. pathway_anchors.tsv changed shape when the anchor was
# separated from the target it seeds:
#
#   old:  step  gene  at_locus  claim_kind  outparalog
#   new:  pathway  target  anchor_gene  anchor_locus  target_level
#         claim_kind  outparalog  in_kegg  scope_reason
#
# 03_homology.py reads the table by column NAME and 04_families.sh reads it by
# column NUMBER, so both break silently or loudly depending on which. A missing
# dict key raises, which is fine. A wrong awk field number does not raise, it
# just selects the wrong text, which is the failure mode this project has been
# bitten by twice. Checking is cheap and the check is the point.
#
# Safe to run more than once: it reports what is already correct and changes
# nothing in that case.
set -euo pipefail
: "${TARO_CODE:?source envs/activate.sh first}"
cd "$TARO_CODE"

ANCH="results/tables/pathway_anchors.tsv"
[ -f "$ANCH" ] || { echo "$ANCH not found; run 02_pathway_set.py first"; exit 1; }

echo "=============================================================="
echo "ANCHORS TABLE"
echo "=============================================================="
head -1 "$ANCH" | tr '\t' '\n' | nl -w3 -s'  '
echo
echo "  rows: $(( $(wc -l < "$ANCH") - 1 ))"
echo
awk -F'\t' 'NR>1{c[$6]++} END{for(k in c) printf "  %3d  %s\n", c[k], k}' "$ANCH"
echo
echo "  copy-number targets with an outparalog declared, which is what 04_ runs:"
awk -F'\t' 'NR>1 && $6=="copy-number" && $7!="" {printf "    %-14s anchor %-10s outparalog %s\n", $2, $4, $7}' "$ANCH"
echo
echo "  copy-number targets with no outparalog, which 04_ skips:"
awk -F'\t' 'NR>1 && $6=="copy-number" && $7=="" {printf "    %s\n", $2}' "$ANCH"
echo

echo "=============================================================="
echo "SCRIPT / SCHEMA AGREEMENT"
echo "=============================================================="

py_ok=0
grep -q 'a\["anchor_locus"\]' scripts/03_homology.py && py_ok=1
if [ "$py_ok" = 1 ]; then
  echo "  03_homology.py   reads anchor_locus / target        OK"
else
  echo "  03_homology.py   still reads at_locus / gene        PATCHING"
fi

awk_ok=0
grep -q '\$6=="copy-number"' scripts/04_families.sh && awk_ok=1
if [ "$awk_ok" = 1 ]; then
  echo "  04_families.sh   awk reads \$2 \$4 \$7                  OK"
else
  echo "  04_families.sh   awk still reads \$2 \$3 \$5            PATCHING"
fi

if [ "$py_ok" = 0 ] || [ "$awk_ok" = 0 ]; then
  echo
  bash scripts/patch_schema.sh
fi

echo
echo "=============================================================="
echo "INPUTS 03_ NEEDS"
echo "=============================================================="
need=(
  "$TARO_REFS/proteomes/orthofinder_input/colesc.fa"
  "$TARO_REFS/proteomes/orthofinder_input/aratha.fa"
  "$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
)
missing=0
for f in "${need[@]}"; do
  if [ -s "$f" ]; then
    printf "  OK       %-58s %s\n" "$(basename "$f")" "$(du -h "$f" | cut -f1)"
  else
    printf "  MISSING  %s\n" "$f"
    missing=1
  fi
done

command -v diamond >/dev/null && echo "  OK       diamond  $(diamond --version 2>&1 | head -1)" \
  || { echo "  MISSING  diamond"; missing=1; }

echo
if [ "$missing" = 1 ]; then
  echo "  Something 03_ reads is absent. Do not run it yet."
  exit 1
fi
echo "  Ready. Next:"
echo
echo "    bash scripts/03_homology.sh 2>&1 | tee results/tables/log_03_homology.txt"
