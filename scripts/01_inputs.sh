#!/usr/bin/env bash
# 01_inputs.sh — proteomes for 17 species, one protein per gene.
#
# Verifies before fetching. Everything was downloaded during the exploratory
# pass and is still on disk; re-downloading 454,114 proteins to prove the
# pipeline runs would waste an hour and change nothing.
#
# One protein per gene, the longest isoform. The gene identifier comes from the
# GFF3 CDS attributes (protein_id= paired with Dbxref=GeneID:), never from the
# protein header. An earlier attempt grouped isoforms by description text, which
# merged unrelated genes sharing a description — "F-box protein",
# "uncharacterized protein" — and produced 14,297 Arabidopsis genes against
# 27,562 expected.
#
# Every count is checked against what NCBI reports. A mismatch stops the
# pipeline rather than being noted, because a silently wrong proteome produces
# silently wrong copy numbers downstream.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
IN="$D/orthofinder_input"
mkdir -p "$IN"

# tag | accession | expected genes | role
SPECIES=(
"colesc|figshare|28253|Colocasia esculenta, the subject"
"manesc|GCF_001659605.2|29735|Manihot esculenta, source of the proposed transgene"
"aratha|GCF_000001735.4|27562|Arabidopsis thaliana, functional reference"
"orysat|GCF_001433935.1|28740|Oryza sativa, published PSY subgroup assignments"
"sorbic|GCF_000003195.3|28236|Sorghum bicolor, published PSY subgroup assignments"
"zeamay|GCF_000005005.2|33461|Zea mays"
"setita|GCF_000263155.2|27505|Setaria italica"
"bradis|GCF_000005505.3|25535|Brachypodium distachyon"
"anacom|GCF_001540865.1|22251|Ananas comosus, non-grass commelinid"
"musacu|GCF_000313855.2|30737|Musa acuminata"
"aspoff|GCF_001876935.1|26460|Asparagus officinalis, closest order to taro here"
"phodac|GCF_009389715.1|29239|Phoenix dactylifera"
"elagui|GCF_000442705.1|26381|Elaeis guineensis"
"zosmar|GCA_001185155.1|20436|Zostera marina, the only other Alismatid available"
"soltub|GCF_000226075.1|28404|Solanum tuberosum"
"nelnuc|GCF_000365185.1|24073|Nelumbo nucifera"
"ambtri|GCF_000471905.2|17106|Amborella trichopoda, outgroup for tree B"
)

printf "%-8s %-20s %9s %9s  %s\n" TAG ACCESSION EXPECTED FOUND STATUS
printf '%.0s-' {1..72}; printf "\n"

FAIL=0
for e in "${SPECIES[@]}"; do
    IFS='|' read -r tag acc want role <<< "$e"
    f="$IN/$tag.fa"

    if [ ! -s "$f" ]; then
        printf "%-8s %-20s %9s %9s  %s\n" "$tag" "$acc" "$want" "-" "MISSING"
        FAIL=1
        continue
    fi

    got=$(grep -c '^>' "$f")
    if [ "$got" -eq "$want" ]; then
        printf "%-8s %-20s %9s %9s  ok\n" "$tag" "$acc" "$want" "$got"
    else
        printf "%-8s %-20s %9s %9s  MISMATCH\n" "$tag" "$acc" "$want" "$got"
        FAIL=1
    fi
done

echo
if [ "$FAIL" -ne 0 ]; then
    cat <<'EOF'
One or more proteomes is missing or has an unexpected gene count.

The pipeline stops here deliberately. A proteome with the wrong number of
sequences produces copy numbers that look plausible and are wrong, and that
failure mode has already occurred once in this project.

To rebuild a species from scratch:
  datasets download genome accession <ACC> --include protein,gff3
and collapse to one protein per gene using the GFF3 CDS attributes. The
archived script scripts/archive/06_rebuild_primary.sh does exactly this.

For taro, the proteome comes from Figshare doi 10.6084/m9.figshare.29917034
(Sun et al. 2026) and collapses from 34,340 transcripts to 28,253 genes on the
>rna-Ces#####.N header pattern.
EOF
    exit 1
fi

TOT=$(cat "$IN"/*.fa | grep -c '^>')
echo "17 species, $TOT genes, all counts match NCBI"

# the Arabidopsis protein -> gene map, needed by every later step
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
if [ ! -s "$MAP" ]; then
    echo
    echo "building the Arabidopsis protein -> locus map"
    G="$TARO_WORK/05_orthology/aratha_gff"
    mkdir -p "$G" && ( cd "$G" \
        && datasets download genome accession GCF_000001735.4 --include gff3 \
             --filename a.zip --no-progressbar \
        && unzip -oq a.zip && find . -name genomic.gff -exec mv {} ./genomic.gff \; \
        && rm -rf ncbi_dataset a.zip README.md md5sum.txt 2>/dev/null || true )
    python3 - "$G/genomic.gff" "$MAP" <<'PY'
import sys, re
gff, out = sys.argv[1:3]
seen = {}
for ln in open(gff):
    if ln.startswith('#'): continue
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'CDS': continue
    p = re.search(r'protein_id=([^;\n]+)', f[8])
    g = re.search(r'Dbxref=[^;\n]*GeneID:(\d+)', f[8])
    s = re.search(r';gene=([^;\n]+)', f[8])
    l = re.search(r'locus_tag=([^;\n]+)', f[8])
    if p and g:
        seen[p.group(1)] = (g.group(1), s.group(1) if s else '', l.group(1) if l else '')
with open(out, 'w') as fh:
    fh.write("protein\tgeneid\tsymbol\tlocus\n")
    for p, (g, s, l) in sorted(seen.items()):
        fh.write(f"{p}\t{g}\t{s}\t{l}\n")
print(f"  {len(seen)} proteins mapped to loci")
PY
else
    echo "Arabidopsis map present: $(($(wc -l < "$MAP") - 1)) proteins"
fi

echo
echo "next: python scripts/02_pathway_set.py"
