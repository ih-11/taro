#!/usr/bin/env bash
# patch_schema.sh — update 03_ and 04_ for the new pathway_anchors.tsv schema.
#
# The table now separates the anchor from the target it seeds:
#   1 pathway  2 target  3 anchor_gene  4 anchor_locus  5 target_level
#   6 claim_kind  7 outparalog  8 in_kegg  9 scope_reason
set -euo pipefail
cd ~/Code/taro/scripts

python3 - <<'PY'
import re
from pathlib import Path

# ---- 03_homology.py -------------------------------------------------------
p = Path("03_homology.py"); s = p.read_text()
s = s.replace('gene, locus = a["gene"], a["at_locus"].upper()',
              'gene, locus = a["target"], a["anchor_locus"].upper()')
s = s.replace('rows.append(dict(step=a["step"], gene=gene, at_locus=locus,',
              'rows.append(dict(step=a["pathway"], gene=gene, at_locus=locus,')
s = s.replace("print(f\"{a['step']:<10}", "print(f\"{a['pathway']:<10}")
s = s.replace("{a['step'] if first else '':<10}", "{a['pathway'] if first else '':<10}")
s = s.replace('gene = a["gene"]', 'gene = a["target"]')
s = s.replace('w.writerow([a["gene"]] + counts', 'w.writerow([a["target"]] + counts')
p.write_text(s)
print("03_homology.py: columns updated")

# ---- 04_families.sh -------------------------------------------------------
p = Path("04_families.sh"); s = p.read_text()
s = s.replace(
  "awk -F'\\t' 'NR>1 && $4==\"copy-number\" {print $2\"\\t\"$3\"\\t\"$5}' \"$ANCH\"",
  "awk -F'\\t' 'NR>1 && $6==\"copy-number\" {print $2\"\\t\"$4\"\\t\"$7}' \"$ANCH\"")
p.write_text(s)
print("04_families.sh: awk columns updated")
PY

echo
echo "verify:"
grep -n 'anchor_locus\|a\["target"\]' 03_homology.py | head -4
grep -n "awk -F'\\\\t' 'NR>1" 04_families.sh
