#!/usr/bin/env bash
# 21_subgroup_resolve.sh — two trees, one for identity and one for resolution.
#
# WHY TWO TREES
#   Squalene synthase is needed to decide which sequences are phytoene
#   synthase, because it is the nearest outparalog and has known identity.
#   But it is divergent enough that including it cost 63% of alignment columns
#   in 18_, against 45% for the PSY-only alignment. Those lost columns are
#   exactly the ones that would resolve relationships among close relatives.
#
#   So: tree A includes squalene synthase and settles identity. Tree B drops
#   it, keeps only the sequences tree A placed on the PSY side, and is rooted
#   on Amborella. Tree B has more sites and more monocot taxa, which is where
#   the resolution has to come from.
#
#   Running one tree for both purposes, as 18_ did, trades away the resolution
#   needed for the question actually being asked.
#
# SUBGROUP ANCHORS, from Lisboa et al. 2022 Table 1
#   The MSU and Phytozome locus identifiers encode the chromosome, so our
#   sequences can be labelled from genomic position without resolving an
#   identifier mapping.
#
#     Oryza sativa      LOC_Os12g43130   chr12  M1
#                       LOC_Os09g38320   chr9   M2
#                       LOC_Os06g51290   chr6   M2
#     Sorghum bicolor   Sobic.008G180800 chr8   M1
#                       Sobic.002G292600 chr2   M2
#
#   Two species with independent anchors is the point. If taro groups with M2
#   in both, the assignment does not rest on rice alone.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

IQ=$(command -v iqtree2 || command -v iqtree)
D="$TARO_REFS/proteomes"
W="$TARO_WORK/05_orthology/trees"
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
cd "$W"

# ---------------------------------------------------------------- search
echo "=========================================================="
echo "COLLECT  PSY and squalene synthase across all species"
echo "=========================================================="

getq () {  # locus -> query fasta appended
    local loc="$1" tag="$2"
    for A in $(awk -F'\t' -v l="$loc" '$4==l{print $1}' "$MAP"); do
        python3 - "$D/aratha.primary.faa" "$A" "$tag" <<'PY'
import sys
faa, acc, tag = sys.argv[1:4]
keep, buf = False, []
for ln in open(faa):
    if ln.startswith('>'):
        keep = ln.rstrip().endswith('|' + acc)
        if keep: buf.append(f">{tag}_aratha_{acc}\n")
    elif keep:
        buf.append(ln)
sys.stdout.writelines(buf)
PY
    done
}

: > q_psy.faa; getq AT5G17230 PSY  >> q_psy.faa
: > q_sqs.faa; getq AT4G34640 SQS  >> q_sqs.faa; getq AT4G34650 SQS >> q_sqs.faa
echo "  PSY queries: $(grep -c '^>' q_psy.faa)   SQS queries: $(grep -c '^>' q_sqs.faa)"

for f in psy sqs; do
    diamond blastp -q "q_${f}.faa" -d all -o "h_${f}.tsv" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 300 --evalue 1e-20 --quiet
done
echo "  PSY hits: $(cut -f2 h_psy.tsv | sort -u | wc -l)"
echo "  SQS hits: $(cut -f2 h_sqs.tsv | sort -u | wc -l)"

python3 - all.faa h_psy.tsv h_sqs.tsv A_psy_sqs.faa <<'PY'
import sys
allf, hp, hs, out = sys.argv[1:5]
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

psy, sqs = hits(hp), hits(hs)
sqs_only = [s for s in sqs if s not in psy]
# a sequence hitting both is kept once, on the PSY side, and tree A decides
with open(out, 'w') as fh:
    for k in psy:
        fh.write(f">{k}\n{pool[k]}\n")
    for k in sqs_only:
        fh.write(f">SQS_{k}\n{pool[k]}\n")
print(f"  tree A input: {len(psy)} PSY-family + {len(sqs_only)} outgroup "
      f"= {len(psy) + len(sqs_only)}")
PY

# ---------------------------------------------------------------- tree A
echo
echo "=========================================================="
echo "TREE A  identity: which sequences are phytoene synthase?"
echo "=========================================================="
famsa -t "${TARO_THREADS:-24}" A_psy_sqs.faa A_psy_sqs.aln 2>/dev/null
trimal -in A_psy_sqs.aln -out A_psy_sqs.trim.aln -automated1 2>/dev/null
B=$(awk '/^>/{next}{n+=length($0)}END{print n}' A_psy_sqs.aln)
A=$(awk '/^>/{next}{n+=length($0)}END{print n}' A_psy_sqs.trim.aln)
echo "  alignment: $B -> $A ($(( 100 * A / B ))% retained)"
"$IQ" -s A_psy_sqs.trim.aln -m MFP -B 1000 -alrt 1000 \
      -T "${TARO_THREADS:-24}" -pre A_iq -quiet -redo
echo "  model: $(grep -m1 'Best-fit model' A_iq.log | sed 's/.*: //')"

python3 - A_iq.treefile B_psy.faa A_psy_sqs.faa <<'PY'
import sys
from Bio import Phylo
tree, out, fa = sys.argv[1:4]

t = Phylo.read(tree, "newick"); t.root_at_midpoint()
tips = t.get_terminals()
sqs = {l.name for l in tips if l.name and l.name.startswith("SQS_")}

best = None
for nd in t.get_nonterminals():
    names = {l.name for l in nd.get_terminals()}
    if sqs <= names and (best is None or len(names) < len(best.get_terminals())):
        best = nd

if best is None or len(best.get_terminals()) == len(tips):
    print("  outgroup not monophyletic — tree A cannot define the boundary")
    print("  falling back to: PSY side = everything not prefixed SQS_")
    keep = {l.name for l in tips if l.name not in sqs}
else:
    inside = {l.name for l in best.get_terminals()}
    intr = sorted(n for n in inside if n not in sqs)
    print(f"  outgroup clade: {len(inside)} tips, "
          f"{len(intr)} query sequences inside it")
    for n in intr:
        print(f"    reclassified as squalene synthase: {n}")
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
print(f"\n  PSY sequences carried into tree B: {len(keep)}")
for s in sorted(by):
    print(f"    {s:<8} {len(by[s]):>2}  {', '.join(sorted(by[s]))}")
PY

# ---------------------------------------------------------------- tree B
echo
echo "=========================================================="
echo "TREE B  resolution: PSY only, Amborella outgroup"
echo "=========================================================="
famsa -t "${TARO_THREADS:-24}" B_psy.faa B_psy.aln 2>/dev/null
trimal -in B_psy.aln -out B_psy.trim.aln -automated1 2>/dev/null
B=$(awk '/^>/{next}{n+=length($0)}END{print n}' B_psy.aln)
A=$(awk '/^>/{next}{n+=length($0)}END{print n}' B_psy.trim.aln)
echo "  alignment: $B -> $A ($(( 100 * A / B ))% retained)"
echo "  compare with tree A above: dropping the outgroup buys back sites"

OG=$(grep '^>ambtri' B_psy.trim.aln | head -1 | sed 's/^>//')
if [ -n "$OG" ]; then
    echo "  outgroup: $OG"
    "$IQ" -s B_psy.trim.aln -m MFP -B 1000 -alrt 1000 -o "$OG" \
          -T "${TARO_THREADS:-24}" -pre B_iq -quiet -redo
else
    echo "  no Amborella sequence; midpoint rooting"
    "$IQ" -s B_psy.trim.aln -m MFP -B 1000 -alrt 1000 \
          -T "${TARO_THREADS:-24}" -pre B_iq -quiet -redo
fi
echo "  model: $(grep -m1 'Best-fit model' B_iq.log | sed 's/.*: //')"

# ---------------------------------------------------------------- anchors
echo
echo "=========================================================="
echo "ANCHORS  label rice and sorghum by chromosome"
echo "=========================================================="

chrom_of () {  # tag accession -> chromosome name
    local tag="$1" acc="$2"
    local g="$TARO_WORK/05_orthology/${tag}_gff/genomic.gff"
    [ -s "$g" ] || { echo "?"; return; }
    local seqid
    seqid=$(grep -m1 "protein_id=${acc}" "$g" 2>/dev/null | cut -f1)
    [ -z "$seqid" ] && { echo "?"; return; }
    awk -v s="$seqid" 'BEGIN{FS="\t"} $1==s && $3=="region" {
        match($9,/chromosome=[^;]+/)
        if (RSTART) print substr($9,RSTART+11,RLENGTH-11); exit}' "$g"
}

: > anchors.tsv
for tag in orysat sorbic; do
    [ -s "$TARO_WORK/05_orthology/${tag}_gff/genomic.gff" ] || {
        echo "  $tag: no GFF3 on disk, skipping"; continue; }
    for ACC in $(grep "^>${tag}_" B_psy.faa | sed "s/^>${tag}_//"); do
        C=$(chrom_of "$tag" "$ACC")
        SG="-"
        if [ "$tag" = orysat ]; then
            case "$C" in 12) SG=M1;; 9|6) SG=M2;; esac
        else
            case "$C" in 8) SG=M1;; 2) SG=M2;; esac
        fi
        printf "%s\t%s\t%s\t%s\n" "$tag" "$ACC" "chr${C}" "$SG" >> anchors.tsv
    done
done
printf "  %-8s %-20s %-8s %s\n" SPECIES PROTEIN CHR SUBGROUP
awk -F'\t' '{printf "  %-8s %-20s %-8s %s\n",$1,$2,$3,$4}' anchors.tsv

echo
echo "  tree B: $W/B_iq.treefile"
echo "  next:   python scripts/22_read_subgroup.py"
