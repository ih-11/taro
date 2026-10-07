#!/usr/bin/env bash
# 04_probe_yin_set.sh — outgroup proteomes, Yin et al. 2021 species set.
#
# No Araceae relative has a protein set at NCBI (03_probe_outgroups.sh).
# Yin et al. 2021 (Mol Ecol Resour 21:68) clustered taro gene families against
# ten species and published the accessions. Most are RefSeq, so annotated.
# Reusing their set anchors our orthology to a published taro analysis.
set -uo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"
T="${TMPDIR:-/tmp}"

# accession | tag | species  (as listed in Yin et al. 2021 methods)
SET=(
"GCA_008360905.1|spipol|Spirodela polyrhiza|closest relative, ~73 MYA"
"GCA_001185155.1|zosmar|Zostera marina|monocot, seagrass"
"GCF_001433935.1|orysat|Oryza sativa|monocot reference"
"GCF_000005005.2|zeamay|Zea mays|monocot reference"
"GCF_000313855.2|musacu|Musa acuminata|monocot, starchy"
"GCF_000226075.1|soltub|Solanum tuberosum|eudicot, tuber — carotenoid literature"
"GCF_000365185.1|nelnuc|Nelumbo nucifera|basal eudicot"
"GCF_000471905.2|ambtri|Amborella trichopoda|basal angiosperm outgroup"
"GCF_000001735.4|aratha|Arabidopsis thaliana|PSY/OR functional reference"
)

printf "\n%-18s %-8s %-10s %8s  %s\n" ACCESSION TAG STATUS PROTEINS SPECIES
printf '%.0s-' {1..78}; printf "\n"

for e in "${SET[@]}"; do
  IFS='|' read -r acc tag sp note <<< "$e"
  J="$T/${tag}.json"
  datasets summary genome accession "$acc" > "$J" 2>/dev/null || { 
      printf "%-18s %-8s %-10s %8s  %s\n" "$acc" "$tag" "QUERYFAIL" "-" "$sp"; continue; }

  python - "$J" "$acc" "$tag" "$sp" <<'PY'
import json, sys
_, j, acc, tag, sp = sys.argv
d = json.load(open(j))
r = (d.get("reports") or [{}])[0]
ai = r.get("annotation_info")
cnt = ai.get("stats", {}).get("gene_counts", {}).get("protein_coding", "?") if ai else "-"
print(f"{acc:<18} {tag:<8} {'PROTEINS' if ai else 'no annot':<10} {str(cnt):>8}  {sp}")
PY
done

cat <<'TXT'

PROTEINS -> fetchable now
no annot -> substitute, or take the proteome from Phytozome / the paper

Minimum viable set for counting taro PSY paralogs:
  Arabidopsis (functional anchor) + rice or maize (monocot) + one more.
More outgroups sharpen orthogroup boundaries; they do not change the count.
TXT
