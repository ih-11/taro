#!/usr/bin/env bash
# 20_add_monocots.sh — widen the monocot sampling.
#
# WHY
#   Taro's PSY sits in a four-tip clade with banana, rice and maize at
#   SH-aLRT 74.8 and UFBoot 59, which is below the threshold for stating a
#   subgroup assignment. Low support on a short internal branch usually means
#   too few informative sites and too few taxa breaking up long branches,
#   rather than a genuine polytomy. More monocots is the cheapest remedy.
#
# WHICH, AND WHY EACH
#   Sorghum bicolor is the most valuable addition, because Lisboa et al. 2022
#   Table 1 assigns its genes to subgroups just as they do for rice:
#       Sobic.008G180800  chr8  M1
#       Sobic.002G292600  chr2  M2
#   That gives a second labelled anchor pair, independent of rice, so the
#   assignment no longer rests on one species.
#
#   Asparagus officinalis is an Asparagales, which is phylogenetically closer
#   to taro than any grass. The nearest available relative to taro remains
#   Zostera marina, and nothing fixes that, but Asparagales sits between the
#   grasses and the Alismatales and should break up the branch leading to taro.
#
#   Ananas comosus is a non-grass commelinid, Setaria italica and Brachypodium
#   distachyon are grasses that subdivide the Poales branch, and the two palms
#   are Arecales. None is essential; each adds a tip where there is currently
#   a long edge.
#
# This script probes first and downloads only what NCBI actually annotates,
# because several monocot genomes are assembly-only. Nothing is assumed.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"; mkdir -p "$D"; cd "$D"
T="${TMPDIR:-/tmp}"

# accession | tag | species | why
CAND=(
"GCF_000003195.3|sorbic|Sorghum bicolor|Poales; published M1 and M2 anchors"
"GCF_001876935.1|aspoff|Asparagus officinalis|Asparagales; closest order to taro here"
"GCF_001540865.1|anacom|Ananas comosus|Poales, non-grass commelinid"
"GCF_000263155.2|setita|Setaria italica|Poales; subdivides the grass branch"
"GCF_000005505.3|bradis|Brachypodium distachyon|Poales; subdivides the grass branch"
"GCF_009389715.1|phodac|Phoenix dactylifera|Arecales"
"GCF_000442705.1|elagui|Elaeis guineensis|Arecales"
)

echo "=========================================================="
echo "PROBE  which of these does NCBI annotate?"
echo "=========================================================="
printf "\n%-18s %-8s %-10s %9s  %s\n" ACCESSION TAG STATUS PROTEINS SPECIES
printf '%.0s-' {1..86}; printf "\n"

AVAIL=()
for e in "${CAND[@]}"; do
    IFS='|' read -r acc tag sp why <<< "$e"
    J="$T/probe_${tag}.json"
    if ! datasets summary genome accession "$acc" > "$J" 2>/dev/null; then
        printf "%-18s %-8s %-10s %9s  %s\n" "$acc" "$tag" "QUERYFAIL" "-" "$sp"
        continue
    fi
    OUT=$(python3 - "$J" "$acc" "$tag" "$sp" <<'PY'
import json, sys
_, j, acc, tag, sp = sys.argv
d = json.load(open(j))
r = (d.get("reports") or [{}])[0]
ai = r.get("annotation_info")
cnt = ai.get("stats", {}).get("gene_counts", {}).get("protein_coding", "?") if ai else "-"
print(f"{acc:<18} {tag:<8} {'PROTEINS' if ai else 'no annot':<10} {str(cnt):>9}  {sp}")
print("YES" if ai else "NO")
PY
)
    echo "$OUT" | head -1
    [ "$(echo "$OUT" | tail -1)" = "YES" ] && AVAIL+=("$acc|$tag|$sp")
done

echo
echo "  annotated and usable: ${#AVAIL[@]} of ${#CAND[@]}"
[ "${#AVAIL[@]}" -eq 0 ] && { echo "  nothing to add"; exit 0; }

echo
echo "=========================================================="
echo "FETCH  one protein per gene, same rule as 06_"
echo "=========================================================="

for e in "${AVAIL[@]}"; do
    IFS='|' read -r acc tag sp <<< "$e"
    if [ -s "$tag.primary.faa" ]; then
        echo "have  $tag  ($(grep -c '^>' "$tag.primary.faa") genes)"
        cp "$tag.primary.faa" "orthofinder_input/$tag.fa"
        continue
    fi
    echo "get   $tag  ($sp)"
    datasets download genome accession "$acc" --include protein,gff3 \
        --filename "$tag.zip" --no-progressbar
    unzip -oq "$tag.zip" -d "${tag}_x"
    FAA=$(find "${tag}_x" -name 'protein.faa' | head -1)
    GFF=$(find "${tag}_x" -name 'genomic.gff' | head -1)

    # keep the GFF3: chromosome positions are needed to label the sorghum
    # subgroup anchors, exactly as rice was labelled in 18_
    mkdir -p "$TARO_WORK/05_orthology/${tag}_gff"
    cp "$GFF" "$TARO_WORK/05_orthology/${tag}_gff/genomic.gff"

    python3 - "$FAA" "$GFF" "$tag.primary.faa" "$tag" <<'PY'
import sys, re
from collections import defaultdict
faa, gff, out, tag = sys.argv[1:5]

p2g = {}
for ln in open(gff):
    if ln.startswith('#'): continue
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'CDS': continue
    p = re.search(r'protein_id=([^;\n]+)', f[8])
    g = (re.search(r'Dbxref=[^;\n]*GeneID:(\d+)', f[8])
         or re.search(r'locus_tag=([^;\n]+)', f[8]))
    if p and g: p2g[p.group(1)] = g.group(1)

seqs, name = {}, None
for ln in open(faa):
    if ln.startswith('>'):
        name = ln[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(ln.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

grp = defaultdict(list)
for acc, s in seqs.items():
    grp[p2g.get(acc, acc)].append((len(s), acc, s))
with open(out, 'w') as fh:
    for g in sorted(grp):
        _, acc, s = max(grp[g])
        fh.write(f">{tag}|{acc}\n")
        for i in range(0, len(s), 60): fh.write(s[i:i+60] + "\n")
print(f"      {len(seqs)} proteins -> {len(grp)} genes")
PY
    cp "$tag.primary.faa" "orthofinder_input/$tag.fa"
    rm -rf "${tag}_x" "$tag.zip"
done

echo
echo "=========================================================="
echo "REBUILD  combined search database"
echo "=========================================================="

W="$TARO_WORK/05_orthology/trees"; mkdir -p "$W"; cd "$W"
: > all.faa
BASE="colesc manesc aratha orysat zeamay musacu soltub nelnuc ambtri zosmar"
for s in $BASE; do cat "$D/orthofinder_input/$s.fa" >> all.faa; done
NEW=""
for e in "${AVAIL[@]}"; do
    IFS='|' read -r _ tag _ <<< "$e"
    cat "$D/orthofinder_input/$tag.fa" >> all.faa
    NEW="$NEW $tag"
done
diamond makedb --in all.faa -d all --quiet

NSP=$(( $(echo $BASE | wc -w) + ${#AVAIL[@]} ))
echo "  $(grep -c '^>' all.faa) proteins, $NSP species"
echo "  added:$NEW"
echo
echo "  next: bash scripts/21_subgroup_resolve.sh"
