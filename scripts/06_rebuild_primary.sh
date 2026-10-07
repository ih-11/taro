#!/usr/bin/env bash
# 06_rebuild_primary.sh — rebuild outgroup proteomes using real gene IDs.
#
# 05_ grouped isoforms by description text, which merged unrelated genes
# sharing a description ("F-box protein", "uncharacterized protein"...).
# Arabidopsis came out 14,297 vs 27,562 expected — about half.
#
# Protein FASTA headers carry no gene ID, so take it from the GFF3:
# CDS features have protein_id= and Dbxref=GeneID:#####.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"; cd "$D"

SET=(
"GCF_000001735.4|aratha|27562"
"GCF_001433935.1|orysat|28738"
"GCF_000005005.2|zeamay|33461"
"GCF_000313855.2|musacu|30737"
"GCF_000226075.1|soltub|28404"
"GCF_000365185.1|nelnuc|24073"
"GCF_000471905.2|ambtri|17106"
"GCA_001185155.1|zosmar|20436"
)

for e in "${SET[@]}"; do
  IFS='|' read -r acc tag want <<< "$e"
  echo "== $tag (expect ~$want genes)"

  datasets download genome accession "$acc" --include protein,gff3 \
      --filename "$tag.zip" --no-progressbar
  unzip -oq "$tag.zip" -d "${tag}_x"
  FAA=$(find "${tag}_x" -name 'protein.faa'  | head -1)
  GFF=$(find "${tag}_x" -name 'genomic.gff'  | head -1)

  python - "$FAA" "$GFF" "$tag.primary.faa" "$tag" <<'PY'
import sys, re
from collections import defaultdict
faa, gff, out, tag = sys.argv[1:5]

# protein accession -> gene id, straight from the GFF3
p2g = {}
for ln in open(gff):
    if ln.startswith('#'): continue
    f = ln.split('\t')
    if len(f) < 9 or f[2] != 'CDS': continue
    a = f[8]
    pid = re.search(r'protein_id=([^;\n]+)', a)
    if not pid: continue
    gid = (re.search(r'Dbxref=[^;\n]*GeneID:(\d+)', a)
           or re.search(r'locus_tag=([^;\n]+)', a)
           or re.search(r'gene=([^;\n]+)', a))
    if gid:
        p2g[pid.group(1)] = gid.group(1)

seqs, name = {}, None
for ln in open(faa):
    if ln.startswith('>'):
        name = ln[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(ln.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

grp = defaultdict(list)
unmapped = 0
for acc, s in seqs.items():
    g = p2g.get(acc)
    if g is None:
        g = acc; unmapped += 1          # keep it; worst case one gene per protein
    grp[g].append((len(s), acc, s))

with open(out, 'w') as fh:
    for g in sorted(grp):
        _, acc, s = max(grp[g])
        fh.write(f">{tag}|{acc}\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i+60] + "\n")

print(f"   {len(seqs):6d} proteins -> {len(grp):6d} genes"
      + (f"   ({unmapped} unmapped)" if unmapped else ""))
PY

  rm -rf "${tag}_x" "$tag.zip"
done

# refresh the OrthoFinder input set
rm -rf orthofinder_input && mkdir -p orthofinder_input
cp colesc.tagged.faa orthofinder_input/colesc.fa
for e in "${SET[@]}"; do
  IFS='|' read -r _ tag _ <<< "$e"
  cp "$tag.primary.faa" "orthofinder_input/$tag.fa"
done

echo
printf "%-14s %8s\n" FILE GENES
for f in orthofinder_input/*.fa; do
  printf "%-14s %8s\n" "$(basename "$f")" "$(grep -c '^>' "$f")"
done
