#!/usr/bin/env bash
# 03_probe_outgroups.sh — which Araceae proteomes can NCBI give us?
# Reports only. Writes nothing but temp JSON.
set -uo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"
T="${TMPDIR:-/tmp}"

for e in "Spirodela polyrhiza|spipol" \
         "Pistia stratiotes|pisstr" \
         "Amorphophallus konjac|amokon" \
         "Zantedeschia elliottiana|zanell" \
         "Colocasia esculenta|colesc"
do
  IFS='|' read -r sp tag <<< "$e"
  printf "\n%-26s (%s)\n" "$sp" "$tag"
  J="$T/${tag}.json"
  datasets summary genome taxon "$sp" --assembly-source all > "$J" 2>/dev/null \
    || { echo "    query failed"; continue; }

  python - "$J" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
n = d.get("total_count", 0)
if not n:
    print("    none at NCBI under this name"); raise SystemExit
reps = d.get("reports", [])
reps.sort(key=lambda r: (not r["accession"].startswith("GCF"),
                         "annotation_info" not in r))
print(f"    {n} assembl{'y' if n==1 else 'ies'}")
for r in reps[:5]:
    ai  = r.get("annotation_info")
    lvl = r.get("assembly_info", {}).get("assembly_level", "?")
    cnt = ai.get("stats", {}).get("gene_counts", {}).get("protein_coding", "?") if ai else "-"
    print(f"      {r['accession']:20s} {lvl:14s} "
          f"{'PROTEINS' if ai else 'no annot':10s} {str(cnt):>7s}  "
          f"{r.get('organism',{}).get('organism_name','?')}")
PY
done

cat <<'TXT'

  PROTEINS  -> datasets download genome accession <acc> --include protein
  no annot  -> assembly only; proteome must come from the paper
  none      -> try NGDC/CNCB (ngdc.cncb.ac.cn) or journal supplement
TXT
