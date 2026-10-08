#!/usr/bin/env bash
# 05_fetch_outgroups.sh — outgroup proteomes, one protein per gene.
#
# Species set from Yin et al. 2021 (Mol Ecol Resour 21:68), who used it for
# taro gene-family clustering. Spirodela polyrhiza — the closest relative —
# has no protein set at NCBI and is omitted; note as a limitation.
#
# NCBI protein.faa carries isoforms. OrthoFinder wants one per gene, so we
# keep the longest per locus, same rule applied to taro in 02_.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"; mkdir -p "$D"; cd "$D"

SET=(
"GCF_000001735.4|aratha|Arabidopsis thaliana"
"GCF_001433935.1|orysat|Oryza sativa"
"GCF_000005005.2|zeamay|Zea mays"
"GCF_000313855.2|musacu|Musa acuminata"
"GCF_000226075.1|soltub|Solanum tuberosum"
"GCF_000365185.1|nelnuc|Nelumbo nucifera"
"GCF_000471905.2|ambtri|Amborella trichopoda"
"GCA_001185155.1|zosmar|Zostera marina"
)

for e in "${SET[@]}"; do
  IFS='|' read -r acc tag sp <<< "$e"
  [ -s "$tag.primary.faa" ] && { echo "have  $tag"; continue; }
  echo "get   $tag  ($sp)"
  datasets download genome accession "$acc" --include protein \
      --filename "$tag.zip" --no-progressbar
  unzip -oq "$tag.zip" -d "${tag}_x"
  find "${tag}_x" -name 'protein.faa' -exec mv {} "$tag.raw.faa" \;
  rm -rf "${tag}_x" "$tag.zip"

  python - "$tag.raw.faa" "$tag.primary.faa" "$tag" <<'PY'
import sys, re
from collections import defaultdict
inp, out, tag = sys.argv[1], sys.argv[2], sys.argv[3]

seqs, name = {}, None
for ln in open(inp):
    if ln.startswith('>'):
        name = ln[1:].rstrip('\n'); seqs[name] = []
    else:
        seqs[name].append(ln.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

# NCBI headers: >XP_123456.1 description [Species]
# Isoforms share a gene but the header does not carry the gene ID, so collapse
# on the description, which is isoform-stable for RefSeq.
grp = defaultdict(list)
for h, s in seqs.items():
    acc = h.split()[0]
    desc = re.sub(r'\s*isoform\s+\S+', '', h.split(' ', 1)[1] if ' ' in h else acc)
    desc = re.sub(r'\s*\[.*\]$', '', desc).strip()
    grp[desc or acc].append((len(s), acc, s))

with open(out, 'w') as fh:
    for k in sorted(grp):
        _, acc, s = max(grp[k])
        fh.write(f">{tag}|{acc}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")
print(f"      {len(seqs)} proteins -> {len(grp)} loci")
PY
  rm -f "$tag.raw.faa"
done

# taro, already collapsed in 02_ — just tag the headers to match
if [ ! -s colesc.tagged.faa ]; then
  sed 's/^>/>colesc|/' colesc.primary.faa > colesc.tagged.faa
fi

mkdir -p "$D/orthofinder_input"
cp colesc.tagged.faa "$D/orthofinder_input/colesc.fa"
for e in "${SET[@]}"; do
  IFS='|' read -r _ tag _ <<< "$e"
  [ -s "$tag.primary.faa" ] && cp "$tag.primary.faa" "$D/orthofinder_input/$tag.fa"
done

echo
echo "OrthoFinder input:"
for f in "$D"/orthofinder_input/*.fa; do
  printf "  %-14s %7s seqs\n" "$(basename "$f")" "$(grep -c '^>' "$f")"
done
