#!/usr/bin/env bash
# 09_pathway_by_geneid.sh
#
# Locate carotenoid pathway orthogroups by joining on Arabidopsis NCBI GeneID.
#
# Why this route. The orthogroup table identifies sequences by protein
# accession, because that is what 06_ wrote into the FASTA headers. The NCBI
# gene query returns a GeneID for a symbol but does not list protein
# accessions. The bridge between the two lives in the Arabidopsis GFF3, where
# each CDS feature carries both protein_id= and Dbxref=GeneID:.
#
# So: download the Arabidopsis GFF3 once, build protein -> gene, resolve each
# pathway symbol to a GeneID, and find which orthogroup holds its proteins.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

W="$TARO_WORK/05_orthology"
GFFDIR="$W/aratha_gff"
MAP="$W/aratha_protein2gene.tsv"

if [ ! -s "$MAP" ]; then
    echo "building Arabidopsis protein -> gene map"
    mkdir -p "$GFFDIR" && cd "$GFFDIR"
    datasets download genome accession GCF_000001735.4 --include gff3 \
        --filename a.zip --no-progressbar
    unzip -oq a.zip
    GFF=$(find . -name 'genomic.gff' | head -1)
    python - "$GFF" "$MAP" <<'PY'
import sys, re
gff, out = sys.argv[1], sys.argv[2]
seen = {}
for ln in open(gff):
    if ln.startswith('#'): continue
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'CDS': continue
    a = f[8]
    p = re.search(r'protein_id=([^;\n]+)', a)
    g = re.search(r'Dbxref=[^;\n]*GeneID:(\d+)', a)
    s = re.search(r';gene=([^;\n]+)', a)
    l = re.search(r'locus_tag=([^;\n]+)', a)
    if p and g:
        seen[p.group(1)] = (g.group(1),
                            s.group(1) if s else '',
                            l.group(1) if l else '')
with open(out, 'w') as fh:
    fh.write("protein\tgeneid\tsymbol\tlocus\n")
    for p, (g, s, l) in sorted(seen.items()):
        fh.write(f"{p}\t{g}\t{s}\t{l}\n")
print(f"  {len(seen)} proteins mapped")
PY
    rm -rf "$GFFDIR/ncbi_dataset" "$GFFDIR/a.zip" "$GFFDIR/README.md" 2>/dev/null || true
fi

wc -l < "$MAP" | xargs echo "map rows:"
python "$TARO_CODE/scripts/09_pathway_report.py"
