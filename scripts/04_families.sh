#!/usr/bin/env bash
# 04_families.sh — family assignment by phylogeny, for copy-number claims only.
#
# Two trees per family, because one tree cannot do both jobs.
#
#   Tree A, identity.   The family plus a declared outparalog of known identity.
#                       Establishes where the family boundary falls.
#   Tree B, resolution. What tree A placed inside the family, without the
#                       outparalog, rooted on Amborella.
#
# The split is not fussiness. A divergent outparalog is what makes the boundary
# testable, and it is also what destroys the alignment columns needed to resolve
# close relatives: including squalene synthase left 50% of columns after
# trimming against 55% without it.
#
# Outparalogs are declared in 02_ before any tree is built, so the boundary test
# is not chosen after seeing where it falls. The CCD4 and NCED outgroup is CCD7
# and CCD8 rather than each other: using NCED to define the CCD4 boundary while
# using CCD4 to define NCED's is close to circular and works only if the two are
# reciprocally monophyletic, which is the thing being tested.
#
# Only families whose claim kind is copy-number get this treatment. Every such
# family gets the same treatment; none is exempted for being less interesting.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

IQ=$(command -v iqtree2 || command -v iqtree)
[ -z "$IQ" ] && { echo "iqtree not found"; exit 1; }

D="$TARO_REFS/proteomes"
W="$TARO_WORK/04_families"
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
ANCH="$TARO_CODE/results/tables/pathway_anchors.tsv"
mkdir -p "$W" && cd "$W"

SPECIES="colesc manesc aratha orysat sorbic zeamay setita bradis anacom \
musacu aspoff phodac elagui zosmar soltub nelnuc ambtri"

if [ ! -f all.dmnd ]; then
    : > all.faa
    for s in $SPECIES; do cat "$D/orthofinder_input/$s.fa" >> all.faa; done
    diamond makedb --in all.faa -d all --quiet
fi
echo "search database: $(grep -c '^>' all.faa) proteins, $(echo $SPECIES | wc -w) species"

# locus -> the protein accession actually present in the proteome
acc_for () {
    local loc="$1"
    for a in $(awk -F'\t' -v l="$loc" 'toupper($4)==toupper(l){print $1}' "$MAP"); do
        grep -q "|${a}\$" "$D/orthofinder_input/aratha.fa" && { echo "$a"; return; }
    done
}

extract () {   # accession tag -> fasta on stdout
    python3 - "$D/orthofinder_input/aratha.fa" "$1" "$2" <<'PY'
import sys
faa, acc, tag = sys.argv[1:4]
keep, buf = False, []
for ln in open(faa):
    if ln.startswith('>'):
        keep = ln.rstrip().endswith('|' + acc)
        if keep: buf.append(f">{tag}aratha_{acc}\n")
    elif keep:
        buf.append(ln)
sys.stdout.writelines(buf)
PY
}

build () {
    local gene="$1" locus="$2" outpar="$3"
    echo
    echo "=========================================================="
    echo "$gene  (anchor $locus)"
    echo "=========================================================="

    local acc; acc=$(acc_for "$locus")
    [ -z "$acc" ] && { echo "  anchor not in proteome"; return; }

    extract "$acc" "" > "${gene}_q.faa"
    : > "${gene}_og_q.faa"
    local have_og=0
    IFS=';' read -ra OPS <<< "$outpar"
    for op in "${OPS[@]}"; do
        [ -z "$op" ] && continue
        local oloc="${op#*:}" oname="${op%%:*}"
        local oacc; oacc=$(acc_for "$oloc")
        [ -z "$oacc" ] && { echo "  outparalog $oname ($oloc) not in proteome"; continue; }
        extract "$oacc" "OG_" >> "${gene}_og_q.faa"
        echo "  outparalog: $oname $oloc -> $oacc"
        have_og=1
    done
    [ "$have_og" -eq 0 ] && { echo "  no usable outparalog; family boundary untestable"; return; }

    for part in q og_q; do
        diamond blastp -q "${gene}_${part}.faa" -d all -o "${gene}_${part}_hits.tsv" \
            --outfmt 6 qseqid sseqid pident length evalue bitscore \
            --max-target-seqs 400 --evalue 1e-20 --quiet
    done
    echo "  family hits: $(cut -f2 "${gene}_q_hits.tsv" | sort -u | wc -l)" \
         " outparalog hits: $(cut -f2 "${gene}_og_q_hits.tsv" | sort -u | wc -l)"

    python3 - all.faa "${gene}_q_hits.tsv" "${gene}_og_q_hits.tsv" "${gene}_A.faa" <<'PY'
import sys
allf, hf, ho, out = sys.argv[1:5]
pool, k = {}, None
for ln in open(allf):
    if ln.startswith('>'):
        k = ln[1:].strip().replace('|', '_'); pool[k] = []
    else: pool[k].append(ln.strip())
pool = {k: ''.join(v) for k, v in pool.items()}

def hits(p, minaln=150):
    out = []
    for ln in open(p):
        f = ln.rstrip('\n').split('\t')
        if int(f[3]) >= minaln:
            out.append(f[1].replace('|', '_'))
    return list(dict.fromkeys(out))

fam, og = hits(hf), hits(ho)
og_only = [s for s in og if s not in fam]
with open(out, 'w') as fh:
    for k in fam: fh.write(f">{k}\n{pool[k]}\n")
    for k in og_only: fh.write(f">OG_{k}\n{pool[k]}\n")
print(f"  tree A: {len(fam)} family candidates + {len(og_only)} outgroup")
PY

    # --- tree A, identity -------------------------------------------------
    famsa -t "${TARO_THREADS:-24}" "${gene}_A.faa" "${gene}_A.aln" 2>/dev/null
    trimal -in "${gene}_A.aln" -out "${gene}_A.trim.aln" -automated1 2>/dev/null
    local b a
    b=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${gene}_A.aln")
    a=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${gene}_A.trim.aln")
    echo "  tree A alignment: $b -> $a ($(( 100 * a / b ))% retained)"
    "$IQ" -s "${gene}_A.trim.aln" -m MFP -B 1000 -alrt 1000 \
          -T "${TARO_THREADS:-24}" -pre "${gene}_A_iq" -quiet -redo
    echo "  tree A model: $(grep -m1 'Best-fit model' "${gene}_A_iq.log" | sed 's/.*: //')"

    python3 - "${gene}_A_iq.treefile" "${gene}_A.faa" "${gene}_B.faa" "$gene" <<'PY'
import sys
from Bio import Phylo
tree, fa, out, gene = sys.argv[1:5]

t = Phylo.read(tree, "newick"); t.root_at_midpoint()
tips = t.get_terminals()
og = {l.name for l in tips if l.name and l.name.startswith("OG_")}

best = None
for nd in t.get_nonterminals():
    names = {l.name for l in nd.get_terminals()}
    if og <= names and (best is None or len(names) < len(best.get_terminals())):
        best = nd

if best is None or len(best.get_terminals()) == len(tips):
    print("  outgroup not monophyletic after midpoint rooting.")
    print("  The boundary cannot be placed from this tree; family membership")
    print("  is reported as untested and no copy-number claim is made.")
    keep = {l.name for l in tips if l.name not in og}
else:
    inside = {l.name for l in best.get_terminals()}
    intr = sorted(n for n in inside if n not in og)
    print(f"  outgroup clade: {len(inside)} tips, "
          f"{len(intr)} family candidates reclassified as outgroup")
    for n in intr:
        print(f"    -> outgroup: {n}")
    keep = {l.name for l in tips if l.name not in inside}

seqs, k = {}, None
for ln in open(fa):
    if ln.startswith('>'):
        k = ln[1:].strip(); seqs[k] = []
    else: seqs[k].append(ln.strip())
with open(out, 'w') as fh:
    for k in keep:
        fh.write(f">{k}\n{''.join(seqs[k])}\n")

by = {}
for n in keep:
    by.setdefault(n.split('_')[0], []).append(n.split('_', 1)[1])
print(f"  in the {gene} family: {len(keep)} sequences")
for s in sorted(by):
    print(f"    {s:<8} {len(by[s]):>2}  {', '.join(sorted(by[s]))}")
PY

    # --- tree B, resolution -----------------------------------------------
    famsa -t "${TARO_THREADS:-24}" "${gene}_B.faa" "${gene}_B.aln" 2>/dev/null
    trimal -in "${gene}_B.aln" -out "${gene}_B.trim.aln" -automated1 2>/dev/null
    b=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${gene}_B.aln")
    a=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${gene}_B.trim.aln")
    echo "  tree B alignment: $b -> $a ($(( 100 * a / b ))% retained)"
    [ $(( 100 * a / b )) -lt 40 ] && \
        echo "  NOTE: below 40% retention; reported with the tree"

    local og2; og2=$(grep '^>ambtri' "${gene}_B.trim.aln" | head -1 | sed 's/^>//')
    if [ -n "$og2" ]; then
        "$IQ" -s "${gene}_B.trim.aln" -m MFP -B 1000 -alrt 1000 -o "$og2" \
              -T "${TARO_THREADS:-24}" -pre "${gene}_B_iq" -quiet -redo
        echo "  tree B rooted on $og2"
    else
        "$IQ" -s "${gene}_B.trim.aln" -m MFP -B 1000 -alrt 1000 \
              -T "${TARO_THREADS:-24}" -pre "${gene}_B_iq" -quiet -redo
        echo "  tree B: no Amborella sequence, midpoint rooted"
    fi
    echo "  tree B model: $(grep -m1 'Best-fit model' "${gene}_B_iq.log" | sed 's/.*: //')"
}

# only copy-number families
awk -F'\t' 'NR>1 && $4=="copy-number" {print $2"\t"$3"\t"$5}' "$ANCH" \
| while IFS=$'\t' read -r gene locus outpar; do
    if [ -z "$outpar" ]; then
        echo
        echo "SKIP $gene: no outparalog declared in 02_."
        echo "  Declare one and record it in docs/METHODS_part1.md before"
        echo "  running, or the family boundary is untestable and the"
        echo "  copy-number claim cannot be made."
        continue
    fi
    build "$gene" "$locus" "$outpar"
done

echo
echo "trees in $W"
echo "next: python scripts/05_report.py"
