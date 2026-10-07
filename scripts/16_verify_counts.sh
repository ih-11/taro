#!/usr/bin/env bash
# 16_verify_counts.sh — two checks that must pass before the counts are quoted.
#
# Check A. Are Ces13250 and Ces13251 two independent CCD4 genes, or one locus
# split across two models? Their gene IDs are adjacent, which is the same
# signature that turned out to be a split model for the PSY-like locus on
# Superscaffold7. The three-copy CCD4 claim depends on the answer.
#
# Check B. Per-species PSY copy number is currently taken from homology search
# hit counts with a 150 aa filter, not from tree topology. Some hits may belong
# to the wider squalene and phytoene synthase superfamily rather than being
# true phytoene synthases. Counting from the tree means counting only sequences
# that fall inside the clade containing the Arabidopsis anchor.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

GFF="$TARO_REFS/proteomes/Colocasia_esculenta.Genome.V1.gff3"
FA="$TARO_REFS/proteomes/orthofinder_input/colesc.fa"
W="$TARO_WORK/05_orthology/trees"

echo "============================================================"
echo "CHECK A — are the three CCD4 genes independent loci?"
echo "============================================================"
echo
printf "%-10s %-16s %12s %12s %3s %8s %7s\n" GENE SCAFFOLD START END STR SIZE_KB AA
for g in Ces13251 Ces13250 Ces03723; do
    read -r chr s e st < <(awk -v g="$g" -F'\t' \
        '$3=="gene" && $9 ~ g {print $1, $4, $5, $7; exit}' "$GFF")
    aa=$(python - "$FA" "$g" <<'PY'
import sys
fa, g = sys.argv[1], sys.argv[2]
cur, n = None, 0
for ln in open(fa):
    if ln.startswith('>'): cur = ln[1:].strip()
    elif cur == f"colesc|{g}": n += len(ln.strip())
print(n)
PY
)
    printf "%-10s %-16s %12d %12d %3s %8.1f %7s\n" \
        "$g" "$chr" "$s" "$e" "$st" "$(echo "($e-$s)/1000" | bc -l)" "$aa"
done

echo
echo "everything annotated between Ces13250 and Ces13251:"
read -r C1 S1 E1 < <(awk -F'\t' '$3=="gene" && $9 ~ /Ces13250/ {print $1,$4,$5; exit}' "$GFF")
read -r C2 S2 E2 < <(awk -F'\t' '$3=="gene" && $9 ~ /Ces13251/ {print $1,$4,$5; exit}' "$GFF")
if [ "$C1" = "$C2" ]; then
    LO=$(( S1 < S2 ? S1 : S2 )); HI=$(( E1 > E2 ? E1 : E2 ))
    GAP=$(( S2 > E1 ? S2 - E1 : S1 - E2 ))
    awk -v c="$C1" -v lo="$LO" -v hi="$HI" -F'\t' \
        '$1==c && $3=="gene" && $4>=lo-10000 && $5<=hi+10000 {
            match($9,/ID=[^;]+/)
            printf "  %12d %12d %s  %s\n", $4, $5, $7, substr($9,RSTART+3,RLENGTH-3)}' "$GFF"
    echo
    echo "  gap between them: ${GAP} bp"
    echo
    echo "  interpretation:"
    echo "    same scaffold, same strand, small gap, nothing between,"
    echo "    protein lengths summing to roughly one -> split model"
    echo "    both near full length, or opposite strands -> tandem duplication"
else
    echo "  different scaffolds ($C1 vs $C2) -> independent loci"
fi

echo
echo
echo "============================================================"
echo "CHECK B — PSY copy number from tree topology"
echo "============================================================"

python - "$W/psy.tree" "$TARO_WORK/05_orthology/aratha_protein2gene.tsv" <<'PY'
from Bio import Phylo
import sys

tree_f, map_f = sys.argv[1], sys.argv[2]

# the Arabidopsis anchor: AT5G17230 is PSY and Arabidopsis has exactly one
anchor_accs = {ln.split('\t')[0] for ln in open(map_f)
               if ln.rstrip('\n').split('\t')[-1] == "AT5G17230"}

t = Phylo.read(tree_f, "newick")
t.root_at_midpoint()
leaves = t.get_terminals()

anchor = next((l for l in leaves
               if l.name and l.name.split('_',1)[1] in anchor_accs), None)
if anchor is None:
    sys.exit("Arabidopsis PSY anchor not in tree")

print(f"\nanchor: {anchor.name}  (AT5G17230, Arabidopsis PSY)")

# Walk outward from the anchor. At each ancestral node, list what the clade
# contains. The true PSY clade is the largest one that still holds exactly one
# Arabidopsis sequence, since Arabidopsis has a single PSY: once a clade picks
# up a second Arabidopsis gene it has expanded into the wider superfamily.
path = t.get_path(anchor)
print("\nwalking outward from the anchor:\n")
print(f"{'depth':>5} {'leaves':>7} {'aratha':>7}  composition")

best = None
for i in range(len(path)-1, -1, -1):
    clade = path[i]
    tips = clade.get_terminals()
    sp = {}
    for l in tips:
        s = l.name.split('_')[0]
        sp[s] = sp.get(s, 0) + 1
    n_ara = sp.get("aratha", 0)
    comp = " ".join(f"{k}:{v}" for k, v in sorted(sp.items()))
    print(f"{i:>5} {len(tips):>7} {n_ara:>7}  {comp}")
    if n_ara == 1:
        best = clade

# also check the root, in case the whole tree is one family
tips = t.get_terminals()
sp = {}
for l in tips:
    s = l.name.split('_')[0]; sp[s] = sp.get(s,0)+1
print(f"{'root':>5} {len(tips):>7} {sp.get('aratha',0):>7}  "
      + " ".join(f"{k}:{v}" for k,v in sorted(sp.items())))

if best is None:
    print("\nno clade with exactly one Arabidopsis sequence; "
          "the anchor's immediate neighbourhood needs manual inspection")
    sys.exit()

SPNAME = {"colesc":"Colocasia esculenta (taro)","manesc":"Manihot esculenta (cassava)",
          "aratha":"Arabidopsis thaliana","orysat":"Oryza sativa","zeamay":"Zea mays",
          "musacu":"Musa acuminata","soltub":"Solanum tuberosum",
          "nelnuc":"Nelumbo nucifera","ambtri":"Amborella trichopoda",
          "zosmar":"Zostera marina"}

tips = best.get_terminals()
by_sp = {}
for l in tips:
    s, g = l.name.split('_', 1)
    by_sp.setdefault(s, []).append(g)

print("\n" + "="*60)
print("PSY clade — largest clade containing exactly one Arabidopsis gene")
print("="*60)
print(f"\n{len(tips)} sequences\n")
print(f"{'species':<34}{'n':>3}  members")
for s in ["colesc","manesc","aratha","orysat","zeamay","musacu",
          "soltub","nelnuc","ambtri","zosmar"]:
    if s in by_sp:
        print(f"{SPNAME[s]:<34}{len(by_sp[s]):>3}  {', '.join(by_sp[s])}")
    else:
        print(f"{SPNAME[s]:<34}{0:>3}  —")

print("\nsequences OUTSIDE this clade (wider superfamily, not counted as PSY):")
inside = {l.name for l in tips}
out = {}
for l in leaves:
    if l.name not in inside:
        s, g = l.name.split('_', 1)
        out.setdefault(s, []).append(g)
if out:
    for s in sorted(out):
        print(f"  {SPNAME.get(s,s):<34}{len(out[s]):>3}  {', '.join(out[s])}")
else:
    print("  none — the whole tree is one clade")
PY
