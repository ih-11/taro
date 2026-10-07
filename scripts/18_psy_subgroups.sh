#!/usr/bin/env bash
# 18_psy_subgroups.sh
#
# Two corrections and one new question, all prompted by Lisboa et al. 2022
# (Genetics and Molecular Biology 45:e20210411), who analysed 351 PSY genes
# across 166 species.
#
# ---------------------------------------------------------------------------
# CORRECTION 1. The clade rule used in 16_ was unsound.
#
# That script took "the largest clade containing exactly one Arabidopsis
# sequence" to be the PSY clade, on the reasoning that Arabidopsis has one PSY.
# Lisboa et al. show the premise does not hold in the way the rule needs:
#
#   "The absence of subgroup E3 PSY in some species, such as those from
#    Brassicaceae, suggests that this paralog was lost in the ancestor of
#    the family."
#
# Arabidopsis lost subgroup E3. A rule anchored on Arabidopsis therefore
# cannot tell "outside the PSY family" apart from "inside a subgroup that
# Arabidopsis does not have".
#
# The rule also never evaluated the root, because Bio.Phylo.get_path returns
# nodes below the root rather than the root itself. The walk printed a root
# clade of 24 tips still holding one Arabidopsis sequence, which by the rule
# as written should have been selected.
#
# Replacing a flawed rule with a different flawed rule is not progress, so
# this script tests the boundary against an outgroup of known identity
# instead. Squalene synthase is the nearest outparalog of phytoene synthase
# and is present in every plant genome. If the five sequences 16_ excluded
# fall with squalene synthase, they are not PSY and the counts stand. If they
# fall inside PSY, the counts were too low.
#
# ---------------------------------------------------------------------------
# CORRECTION 2. The literature disagrees with our verified counts.
#
# Lisboa et al. report three PSY paralogs in cassava, rice, maize and sorghum.
# Our naive hit count gave 3 for each; our "verified" clade count gave 2.
# The naive count agreeing with the literature is a reason to doubt the
# verification, not to prefer it.
#
# ---------------------------------------------------------------------------
# NEW QUESTION. Which monocot subgroup does taro keep?
#
# Angiosperm PSY splits into subgroups E1, E2 and E3 in eudicots and M1 and M2
# in monocots, with the duplications predating the monocot-eudicot divergence.
# Monocots carry M1 and M2 only. A monocot with a single PSY has therefore lost
# one of the two, and which one it lost is a more informative statement than
# the count on its own.
#
# Lisboa et al. Table 1 assigns rice genes to subgroups. The MSU locus
# identifier encodes the chromosome, which is enough to label our rice
# sequences without resolving an identifier mapping:
#
#   LOC_Os12g43130   chr12   M1   shoot, leaf, blade, central vein, flag leaf
#   LOC_Os09g38320   chr9    M2   seedling, coleoptile, leaf, ROOT, radicle,
#                                 panicle, spikelet, floret, stigma, ovary;
#                                 also abiotic-stress inducible
#   LOC_Os06g51290   chr6    M2   seedling, leaf, panicle, floret, ovary, shoot
#
# If taro kept M1 it kept the copy that in rice is biased to photosynthetic
# tissue, and lost the one that reaches root and storage-adjacent tissue. For
# a corm crop that bears directly on whether the organ can accumulate
# carotenoid, and it is a question the collaborator's design does not raise.
#
# Taro is not among the 166 species Lisboa et al. screened. Among monocots
# their single-copy list holds only Zostera marina and Saccharum spontaneum,
# so a taro count would be an addition to it rather than a restatement.
# ---------------------------------------------------------------------------
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

IQ=$(command -v iqtree2 || command -v iqtree)
D="$TARO_REFS/proteomes"
W="$TARO_WORK/05_orthology/trees"
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
cd "$W"

echo "=========================================================="
echo "PART 1  squalene synthase outgroup"
echo "=========================================================="

# Arabidopsis squalene synthase, the nearest outparalog of PSY.
# AT4G34640 = SQS1, AT4G34650 = SQS2. Both are in the Arabidopsis proteome we
# already built, so no new download is needed.
: > sqs_q.faa
for LOC in AT4G34640 AT4G34650; do
    for A in $(awk -F'\t' -v l="$LOC" '$4==l{print $1}' "$MAP"); do
        python3 - "$D/aratha.primary.faa" "$A" "$LOC" >> sqs_q.faa <<'PY'
import sys
faa, acc, loc = sys.argv[1:4]
keep, buf = False, []
for ln in open(faa):
    if ln.startswith('>'):
        keep = ln.rstrip().endswith('|' + acc)
        if keep: buf.append(f">aratha_{acc}\n")
    elif keep:
        buf.append(ln)
sys.stdout.writelines(buf)
PY
    done
done
echo "  squalene synthase queries: $(grep -c '^>' sqs_q.faa)"
grep '^>' sqs_q.faa | sed 's/^/    /'

[ -s sqs_q.faa ] || { echo "  no SQS found in the Arabidopsis proteome"; exit 1; }

# collect squalene synthase homologs across all ten species
diamond blastp -q sqs_q.faa -d all -o sqs_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore \
    --max-target-seqs 100 --evalue 1e-20 --quiet

echo "  SQS homologs found: $(cut -f2 sqs_hits.tsv | sort -u | wc -l)"

echo
echo "=========================================================="
echo "PART 2  rebuild the PSY tree with the outgroup included"
echo "=========================================================="

python3 - all.faa psy.faa sqs_hits.tsv psy_og.faa <<'PY'
import sys
allf, psyf, hits, out = sys.argv[1:5]

def read(p, strip_pipe=False):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith('>'):
            k = ln[1:].strip()
            if strip_pipe: k = k.replace('|', '_')
            d[k] = []
        else:
            d[k].append(ln.strip())
    return {k: ''.join(v) for k, v in d.items()}

pool = read(allf, strip_pipe=True)
psy  = read(psyf)

# squalene synthase homologs, alignment at least 150 aa, same filter the PSY
# set used so the two halves of the tree are comparable
sqs = []
for ln in open(hits):
    f = ln.rstrip('\n').split('\t')
    if int(f[3]) >= 150:
        sqs.append(f[1].replace('|', '_'))
sqs = [s for s in dict.fromkeys(sqs) if s not in psy]

with open(out, 'w') as fh:
    for k, v in psy.items():
        fh.write(f">{k}\n{v}\n")
    for k in sqs:
        fh.write(f">SQS_{k}\n{pool[k]}\n")

print(f"  PSY-family sequences : {len(psy)}")
print(f"  SQS outgroup added   : {len(sqs)}")
print(f"  total                : {len(psy) + len(sqs)}")
PY

famsa -t "${TARO_THREADS:-24}" psy_og.faa psy_og.aln 2>/dev/null
trimal -in psy_og.aln -out psy_og.trim.aln -automated1 2>/dev/null
B=$(awk '/^>/{next}{n+=length($0)}END{print n}' psy_og.aln)
A=$(awk '/^>/{next}{n+=length($0)}END{print n}' psy_og.trim.aln)
echo "  alignment: $B -> $A residues ($(( 100 * A / B ))% retained after trimming)"

OG=$(grep '^>SQS_' psy_og.trim.aln | head -1 | sed 's/^>//')
echo "  outgroup : $OG"
"$IQ" -s psy_og.trim.aln -m MFP -B 1000 -alrt 1000 -o "$OG" \
      -T "${TARO_THREADS:-24}" -pre psy_og_iq -quiet -redo
echo "  model    : $(grep -m1 'Best-fit model' psy_og_iq.log | sed 's/.*: //')"

echo
echo "=========================================================="
echo "PART 3  where does the PSY/SQS boundary actually fall?"
echo "=========================================================="

python3 - psy_og_iq.treefile <<'PY'
import sys, re
from Bio import Phylo

t = Phylo.read(sys.argv[1], "newick")
t.ladderize()
tips = t.get_terminals()

# The five sequences 16_ excluded from the PSY clade. If they are squalene
# synthase they will sit with the SQS_ outgroup; if they are PSY they will not.
QUESTIONED = ["XP_020521936.1", "XP_021608855.1", "XP_015611707.1",
              "XP_006354229.1", "XP_020397011.1"]

def is_sqs(leaf):
    return bool(leaf.name) and leaf.name.startswith("SQS_")

sqs_tips = [l for l in tips if is_sqs(l)]
psy_tips = [l for l in tips if not is_sqs(l)]
print(f"  tips: {len(tips)}  ({len(sqs_tips)} outgroup, {len(psy_tips)} query)")

# the smallest clade holding every outgroup sequence defines the SQS side
best = None
for nd in t.get_nonterminals():
    names = {l.name for l in nd.get_terminals()}
    if all(l.name in names for l in sqs_tips):
        if best is None or len(nd.get_terminals()) < len(best.get_terminals()):
            best = nd

if best is None:
    print("  outgroup is not monophyletic; inspect the tree by hand")
    raise SystemExit

inside = {l.name for l in best.get_terminals()}
contam = sorted(n for n in inside if not n.startswith("SQS_"))
print(f"  smallest clade containing all outgroup sequences: "
      f"{len(inside)} tips")
print(f"  non-outgroup sequences inside it: {len(contam)}")
for n in contam:
    print(f"    {n}")

print()
print("  verdict for each sequence 16_ excluded:")
for q in QUESTIONED:
    hit = next((l for l in tips if l.name and q in l.name), None)
    if hit is None:
        print(f"    {q:<18} not in this tree")
        continue
    side = "SQUALENE SYNTHASE side" if hit.name in inside else "PSY side"
    print(f"    {q:<18} {hit.name.split('_')[0]:<8} {side}")

print()
print("  If these are on the squalene synthase side, the counts from 16_")
print("  stand and the literature disagreement needs another explanation.")
print("  If they are on the PSY side, they are PSY and the counts were low.")

# full composition of the PSY side, which is the corrected count
psy_side = [l for l in tips if l.name not in inside]
by_sp = {}
for l in psy_side:
    sp = l.name.split('_')[0]
    by_sp.setdefault(sp, []).append(l.name.split('_', 1)[1])

NAMES = {"colesc": "Colocasia esculenta (taro)",
         "manesc": "Manihot esculenta (cassava)",
         "aratha": "Arabidopsis thaliana", "orysat": "Oryza sativa",
         "zeamay": "Zea mays", "musacu": "Musa acuminata",
         "soltub": "Solanum tuberosum", "nelnuc": "Nelumbo nucifera",
         "ambtri": "Amborella trichopoda", "zosmar": "Zostera marina"}

print()
print("  corrected PSY count per species (everything off the outgroup side):")
print(f"  {'species':<34}{'n':>3}  members")
for s in ["colesc", "manesc", "aratha", "orysat", "zeamay", "musacu",
          "soltub", "nelnuc", "ambtri", "zosmar"]:
    g = by_sp.get(s, [])
    print(f"  {NAMES[s]:<34}{len(g):>3}  {', '.join(g) if g else '-'}")

print()
print("  Lisboa et al. 2022 report, for comparison:")
print("    cassava, rice, maize, sorghum   3 each")
print("    Arabidopsis, Zostera marina     1 each")
PY

echo
echo "=========================================================="
echo "PART 4  M1 or M2 — which subgroup does taro keep?"
echo "=========================================================="
echo
echo "  Lisboa et al. Table 1, rice subgroup assignments."
echo "  The MSU locus identifier encodes the chromosome, so our rice"
echo "  sequences can be labelled from their genomic position alone."
echo
echo "    LOC_Os12g43130  chr12  M1  shoot, leaf, blade, central vein"
echo "    LOC_Os09g38320  chr9   M2  seedling, coleoptile, leaf, ROOT,"
echo "                               radicle, panicle, ovary; stress inducible"
echo "    LOC_Os06g51290  chr6   M2  seedling, leaf, panicle, floret, ovary"
echo

RG="$TARO_WORK/05_orthology/orysat_gff"
if [ ! -s "$RG/genomic.gff" ]; then
    mkdir -p "$RG" && ( cd "$RG" \
        && datasets download genome accession GCF_001433935.1 --include gff3 \
             --filename r.zip --no-progressbar \
        && unzip -oq r.zip \
        && find . -name genomic.gff -exec mv {} ./genomic.gff \; \
        && rm -rf ncbi_dataset r.zip README.md md5sum.txt 2>/dev/null || true )
fi

echo "  rice PSY sequences in our tree, with chromosome:"
printf "  %-20s %-22s %s\n" PROTEIN SEQUENCE_ID SUBGROUP
for ACC in $(grep '^>orysat_' psy_og.faa | sed 's/^>orysat_//'); do
    LINE=$(grep -m1 "protein_id=${ACC}" "$RG/genomic.gff" 2>/dev/null || true)
    if [ -z "$LINE" ]; then
        printf "  %-20s %-22s %s\n" "$ACC" "-" "not found in GFF3"
        continue
    fi
    SEQID=$(echo "$LINE" | cut -f1)
    # NCBI rice sequence ids are NC_0890xx.x for chromosomes 1..12
    CHR=$(awk -v s="$SEQID" 'BEGIN{FS="\t"} $0 ~ "^##sequence-region" {next}
          $1==s && $3=="region" {match($9,/chromosome=[^;]+/);
          if (RSTART) print substr($9,RSTART+11,RLENGTH-11); exit}' "$RG/genomic.gff")
    case "$CHR" in
        12) SG="M1  (LOC_Os12g43130, shoot/leaf biased)" ;;
        9)  SG="M2  (LOC_Os09g38320, broad incl. root)"  ;;
        6)  SG="M2  (LOC_Os06g51290, broad)"             ;;
        *)  SG="unassigned (chr ${CHR:-?})"              ;;
    esac
    printf "  %-20s %-22s %s\n" "$ACC" "chr${CHR:-?}" "$SG"
done

echo
echo "  Taro's PSY subgroup is read from which labelled rice sequence it"
echo "  groups with in psy_og_iq.treefile. Support at that node decides"
echo "  whether the assignment can be stated or only noted."
echo
echo "  tree: $W/psy_og_iq.treefile"
