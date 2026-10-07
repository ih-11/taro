#!/usr/bin/env bash
# 02_primary_isoform.sh — one protein per gene.
#
# The Figshare .pep carries 34,340 proteins for 28,253 genes: alternative
# transcripts from PASA. OrthoFinder assumes one sequence per gene; isoforms
# inflate orthogroups and would overstate taro PSY paralogy, which is the
# number this project exists to establish.
#
# Headers are >rna-Ces#####.N — gene Ces#####, isoform .N.
# Keep the longest isoform per gene, rename to the gene ID.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
IN="$D/Colocasia_esculenta.Genome.V1.pep"
OUT="$D/colesc.primary.faa"

python - "$IN" "$OUT" <<'PY'
import sys, re
from collections import defaultdict

inp, out = sys.argv[1], sys.argv[2]
seqs, name = {}, None
for line in open(inp):
    if line.startswith('>'):
        name = line[1:].split()[0]
        seqs[name] = []
    else:
        seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

# rna-Ces00003.2 -> gene Ces00003
pat = re.compile(r'^(?:rna-)?(.+?)\.(\d+)$')
by_gene = defaultdict(list)
for tid, s in seqs.items():
    m = pat.match(tid)
    gene = m.group(1) if m else tid
    by_gene[gene].append((len(s), tid, s))

with open(out, 'w') as fh:
    for gene in sorted(by_gene):
        _, tid, s = max(by_gene[gene])          # longest isoform
        fh.write(f">{gene}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")

multi = sum(1 for v in by_gene.values() if len(v) > 1)
print(f"  transcripts in : {len(seqs)}")
print(f"  genes out      : {len(by_gene)}")
print(f"  multi-isoform  : {multi} genes")
print(f"  max isoforms   : {max(len(v) for v in by_gene.values())}")
PY

echo
echo "  stop codons    : $(grep -c '\*' "$OUT" || true)  (OrthoFinder dislikes these)"
sed -i 's/\*//g' "$OUT"
echo "  -> $OUT  ($(grep -c '^>' "$OUT") sequences)"
