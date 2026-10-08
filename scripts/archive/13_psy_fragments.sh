#!/usr/bin/env bash
# 13_psy_fragments.sh — are Ces12496 and Ces12497 one gene or two?
#
# Both matched Arabidopsis PSY but covered only 26% and 34% of the query, and
# both fell below the 150 aa threshold for the gene tree. Their gene IDs are
# adjacent and their coverage sums to roughly one full protein, which is the
# signature of a single gene split across two models by the annotation.
#
# Genome coordinates settle it. Adjacent, same strand, small gap, non
# overlapping alignment regions means one gene. Distant or overlapping means
# two genes.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

GFF="$TARO_REFS/proteomes/Colocasia_esculenta.Genome.V1.gff3"

echo "coordinates:"
for g in Ces24605 Ces12496 Ces12497; do
    awk -v g="$g" -F'\t' '$3=="gene" && $9 ~ g {
        printf "  %-10s %-14s %10d %10d  %s  %6.1f kb\n",
               g, $1, $4, $5, $7, ($5-$4)/1000 }' "$GFF"
done

echo
echo "anything else annotated between them:"
read -r CHR S1 E1 < <(awk -F'\t' '$3=="gene" && $9 ~ /Ces12496/ {print $1, $4, $5}' "$GFF")
read -r _   S2 E2 < <(awk -F'\t' '$3=="gene" && $9 ~ /Ces12497/ {print $1, $4, $5}' "$GFF")
LO=$(( S1 < S2 ? S1 : S2 )); HI=$(( E1 > E2 ? E1 : E2 ))
awk -v c="$CHR" -v lo="$LO" -v hi="$HI" -F'\t' \
    '$1==c && $3=="gene" && $4>=lo-5000 && $5<=hi+5000 {
        match($9,/ID=[^;]+/); printf "  %10d %10d %s  %s\n", $4,$5,$7,substr($9,RSTART+3,RLENGTH-3)}' "$GFF"

echo
echo "protein lengths:"
for g in Ces24605 Ces12496 Ces12497; do
    L=$(awk -v g="$g" '/^>/{p=($0=="colesc|"g)} p&&!/^>/{n+=length($0)} END{print n+0}' \
        "$TARO_REFS/proteomes/orthofinder_input/colesc.fa")
    printf "  %-10s %4d aa\n" "$g" "$L"
done
echo "  Arabidopsis PSY   437 aa  (for comparison)"
