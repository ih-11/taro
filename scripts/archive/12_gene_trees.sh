#!/usr/bin/env bash
# 12_gene_trees.sh — phylogenies for PSY and the CCD/NCED family.
#
# Two questions that best-hit searching cannot answer:
#
#   1. Taro has one full-length PSY (Ces24605) plus two short matches
#      (Ces12497, Ces12496) covering 34% and 26% of the query. Adjacent gene
#      IDs, coverage summing to roughly one protein. Are they a single gene
#      split by the annotation, or a divergent paralog? Only topology decides.
#
#   2. Three taro genes were called CCD4 orthologs, but Ces11733 and Ces27429
#      also matched NCED3. CCD and NCED are the same superfamily and cross-hit
#      readily, so the three-copy CCD4 claim is unsafe until the clades are
#      separated.
#
# Method: collect every family member across all nine species by homology to
# the Arabidopsis anchor, align with FAMSA, build a tree with FastTree.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
W="$TARO_WORK/05_orthology/trees"
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
mkdir -p "$W" && cd "$W"

SPECIES="colesc aratha orysat zeamay musacu soltub nelnuc ambtri zosmar"

# one combined database, so family members are collected from every species
if [ ! -f all.dmnd ]; then
    : > all.faa
    for s in $SPECIES; do cat "$D/orthofinder_input/$s.fa" >> all.faa; done
    diamond makedb --in all.faa -d all --quiet
    echo "combined database: $(grep -c '^>' all.faa) proteins, 9 species"
fi

# family | arabidopsis locus | evalue cutoff
build_tree () {
    local fam="$1" locus="$2" ev="$3"
    echo
    echo "=== $fam  (anchor $locus) ==="

    local acc
    acc=$(awk -F'\t' -v l="$locus" '$4==l{print $1}' "$MAP" \
          | while read -r a; do grep -q "|$a\$" "$D/aratha.primary.faa" && echo "$a" && break; done)
    [ -z "$acc" ] && { echo "  anchor not in proteome"; return; }

    python - "$D/aratha.primary.faa" "$acc" "${fam}_q.faa" <<'PY'
import sys
faa, acc, out = sys.argv[1:4]
keep, buf = False, []
for ln in open(faa):
    if ln.startswith('>'): keep = ln.rstrip().endswith('|' + acc)
    if keep: buf.append(ln)
open(out, 'w').writelines(buf)
PY

    diamond blastp -q "${fam}_q.faa" -d all -o "${fam}_hits.tsv" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 200 --evalue "$ev" --quiet

    echo "  hits: $(wc -l < "${fam}_hits.tsv")"
    cut -f2 "${fam}_hits.tsv" | cut -d'|' -f1 | sort | uniq -c \
        | awk '{printf "    %-8s %s\n", $2, $1}'

    # pull sequences, keeping only reasonably complete matches
    python - "all.faa" "${fam}_hits.tsv" "${fam}.faa" <<'PY'
import sys
allf, hits, out = sys.argv[1:4]
seqs, name = {}, None
for ln in open(allf):
    if ln.startswith('>'):
        name = ln[1:].strip(); seqs[name] = []
    else: seqs[name].append(ln.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}

keep = []
for ln in open(hits):
    f = ln.split('\t')
    if int(f[3]) >= 150:          # alignment at least 150 aa
        keep.append(f[1])
seen = set()
with open(out, 'w') as fh:
    for k in keep:
        if k in seen: continue
        seen.add(k)
        fh.write(f">{k.replace('|','_')}\n{seqs[k]}\n")
print(f"  kept {len(seen)} sequences (alignment >= 150 aa)")
PY

    famsa -t "${TARO_THREADS:-24}" "${fam}.faa" "${fam}.aln" 2>/dev/null
    fasttree -quiet -nosupport < "${fam}.aln" > "${fam}.tree" 2>/dev/null
    echo "  tree: $W/${fam}.tree"
}

build_tree psy AT5G17230 1e-20
build_tree ccd AT4G19170 1e-30      # CCD4; will also pull NCED

echo
echo "taro members in each tree:"
for f in psy ccd; do
    echo "  $f:"
    grep '^>' "${f}.faa" | grep colesc | sed 's/>colesc_/    /'
done
