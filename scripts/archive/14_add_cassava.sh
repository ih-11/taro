#!/usr/bin/env bash
# 14_add_cassava.sh — add cassava to the proteome set.
#
# Cassava is the source of the transgene this collaboration proposes to move
# into taro, and it carries the PSY1/PSY2 pair that motivates the copy number
# question. Without it the PSY tree has no anchor to the gene being
# transferred.
#
# Not added earlier because the species set came from Yin et al. 2021, who
# were studying taro gene families rather than carotenoid engineering.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"; cd "$D"
ACC=GCF_001659605.2; TAG=manesc

if [ ! -s "$TAG.primary.faa" ]; then
    datasets download genome accession "$ACC" --include protein,gff3 \
        --filename "$TAG.zip" --no-progressbar
    unzip -oq "$TAG.zip" -d "${TAG}_x"
    FAA=$(find "${TAG}_x" -name 'protein.faa' | head -1)
    GFF=$(find "${TAG}_x" -name 'genomic.gff' | head -1)

    python - "$FAA" "$GFF" "$TAG.primary.faa" "$TAG" <<'PY'
import sys, re
from collections import defaultdict
faa, gff, out, tag = sys.argv[1:5]
p2g = {}
for ln in open(gff):
    if ln.startswith('#'): continue
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'CDS': continue
    p = re.search(r'protein_id=([^;\n]+)', f[8])
    g = re.search(r'Dbxref=[^;\n]*GeneID:(\d+)', f[8])
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
print(f"  {len(seqs)} proteins -> {len(grp)} genes  (expect ~29735)")
PY
    rm -rf "${TAG}_x" "$TAG.zip"
fi

cp "$TAG.primary.faa" "orthofinder_input/$TAG.fa"

# rebuild the combined search database so the trees can see cassava
W="$TARO_WORK/05_orthology/trees"; mkdir -p "$W"; cd "$W"
: > all.faa
for s in colesc manesc aratha orysat zeamay musacu soltub nelnuc ambtri zosmar; do
    cat "$D/orthofinder_input/$s.fa" >> all.faa
done
diamond makedb --in all.faa -d all --quiet
echo "combined database rebuilt: $(grep -c '^>' all.faa) proteins, 10 species"
