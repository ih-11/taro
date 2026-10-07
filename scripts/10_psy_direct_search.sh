#!/usr/bin/env bash
# 10_psy_direct_search.sh
#
# OrthoFinder placed zero taro genes in the PSY orthogroup (OG0008455), which
# cannot be right: taro accumulates carotenoids and must have a phytoene
# synthase. Zostera marina is also absent from that orthogroup, and 6.6% of
# genes went unassigned in the run, so the likely explanation is that MCL
# split the family rather than that taro lacks the gene.
#
# A direct homology search settles it, because it does not depend on
# clustering at all.
#
# Note on accessions: 06_ kept the longest isoform per gene, so the protein
# present in aratha.primary.faa is not necessarily the first accession listed
# for that locus. We therefore take every protein accession belonging to
# AT5G17230 and use whichever one survived into the file.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
W="$TARO_WORK/05_orthology/psy_check"
mkdir -p "$W" && cd "$W"

MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
mapfile -t ACCS < <(awk -F'\t' '$4=="AT5G17230"{print $1}' "$MAP")
echo "AT5G17230 protein accessions in GFF3: ${ACCS[*]}"

python - "$D/aratha.primary.faa" query_psy.faa "${ACCS[@]}" <<'PY'
import sys
faa, out = sys.argv[1], sys.argv[2]
accs = set(sys.argv[3:])
keep, buf, found = False, [], None
for ln in open(faa):
    if ln.startswith('>'):
        hdr = ln[1:].strip().split('|')[-1]
        keep = hdr in accs
        if keep: found = hdr
    if keep:
        buf.append(ln)
open(out, 'w').writelines(buf)
if found:
    print(f"  query: {found}  ({sum(len(l.strip()) for l in buf[1:])} aa)")
else:
    print("  none of those accessions is in aratha.primary.faa")
    sys.exit(1)
PY

echo
echo "searching taro proteome (28,253 genes)"
diamond makedb --in "$D/orthofinder_input/colesc.fa" -d colesc --quiet
diamond blastp -q query_psy.faa -d colesc -o psy_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 25 --evalue 1e-5 --quiet

echo
if [ ! -s psy_hits.tsv ]; then
    echo "  NO HITS at e<1e-5. That would be extraordinary — investigate the"
    echo "  taro proteome itself before believing it."
    exit 0
fi

printf "%-26s %7s %6s %11s %7s\n" TARO_GENE IDENT ALN EVALUE BITS
awk -F'\t' '{printf "%-26s %6.1f%% %6s %11s %7s\n", $2,$3,$4,$5,$6}' psy_hits.tsv

echo
echo "where OrthoFinder put each hit:"
OGF=$(find "$TARO_WORK/05_orthology" -name 'Orthogroups.tsv' | head -1)
UNA=$(find "$TARO_WORK/05_orthology" -name 'Orthogroups_UnassignedGenes.tsv' | head -1)
cut -f2 psy_hits.tsv | sed 's/^colesc|//' | sort -u | while read -r g; do
    og=$(grep -m1 -F "colesc|$g" "$OGF" 2>/dev/null | cut -f1)
    if [ -n "$og" ]; then
        printf "  %-22s %s\n" "$g" "$og"
    elif grep -q -F "colesc|$g" "$UNA" 2>/dev/null; then
        printf "  %-22s UNASSIGNED\n" "$g"
    else
        printf "  %-22s not located\n" "$g"
    fi
done

echo
echo "reciprocal check — top taro hit back against Arabidopsis:"
TOP=$(head -1 psy_hits.tsv | cut -f2)
python - "$D/orthofinder_input/colesc.fa" "$TOP" top_taro.faa <<'PY'
import sys
faa, acc, out = sys.argv[1:4]
keep, buf = False, []
for ln in open(faa):
    if ln.startswith('>'): keep = ln[1:].strip() == acc
    if keep: buf.append(ln)
open(out,'w').writelines(buf)
PY
diamond makedb --in "$D/aratha.primary.faa" -d aratha --quiet
diamond blastp -q top_taro.faa -d aratha -o recip.tsv \
    --outfmt 6 qseqid sseqid pident evalue bitscore \
    --max-target-seqs 5 --evalue 1e-5 --quiet
awk -F'\t' -v m="$MAP" 'BEGIN{while((getline l < m)>0){split(l,a,"\t"); loc[a[1]]=a[4]}}
    {split($2,s,"|"); printf "  %-18s %-12s %6.1f%% %11s\n", s[2], loc[s[2]], $3, $4}' recip.tsv
