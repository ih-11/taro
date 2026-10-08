#!/usr/bin/env bash
# 17_trees_proper.sh — phylogenies to the standard an evolutionary genomics
# paper would expect.
#
# The trees from 12_ were a screen: FAMSA alignment fed straight to FastTree
# with SH-like local support and midpoint rooting. That was enough to separate
# CCD4 from NCED and to place taro's PSY, but four things fall short of what
# would be reported in a gene family paper.
#
#   1. No alignment trimming. Gappy and poorly aligned columns, especially at
#      the termini, distort branch lengths and can inflate support.
#   2. SH-like local support is a fast approximation, not bootstrap.
#   3. No model selection. FastTree assumes JTT+CAT and tests nothing.
#   4. Midpoint rooting is a convenience. Amborella is in the tree and is the
#      conventional basal angiosperm outgroup.
#
# This script addresses all four: trimAl for trimming, IQ-TREE with
# ModelFinder for model choice, 1000 ultrafast bootstrap replicates plus
# SH-aLRT for support, and explicit outgroup rooting.
#
# Runtime is minutes, not hours, at 19 and 96 sequences.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

IQ=$(command -v iqtree2 || command -v iqtree)
[ -z "$IQ" ] && { echo "iqtree not found"; exit 1; }
echo "iqtree: $IQ"
echo "trimal: $(command -v trimal)"

W="$TARO_WORK/05_orthology/trees"
cd "$W"

run_tree () {
    local fam="$1" outgroup="$2"
    echo
    echo "=================================================="
    echo "$fam"
    echo "=================================================="

    [ -s "${fam}.aln" ] || { echo "  ${fam}.aln missing; run 12_ first"; return; }

    # --- trim ---------------------------------------------------------
    # -automated1 chooses between gappyout and strict based on alignment
    # characteristics. It is the usual choice for maximum likelihood input.
    trimal -in "${fam}.aln" -out "${fam}.trim.aln" -automated1 2>/dev/null

    local before after
    before=$(awk '/^>/{next}{n+=length($0)}END{print n}' "${fam}.aln")
    after=$(awk  '/^>/{next}{n+=length($0)}END{print n}' "${fam}.trim.aln")
    local nseq
    nseq=$(grep -c '^>' "${fam}.trim.aln")
    echo "  sequences            : $nseq"
    echo "  residues before trim : $before"
    echo "  residues after trim  : $after  ($(( 100 * after / before ))% retained)"

    # --- tree ---------------------------------------------------------
    # -m MFP          ModelFinder Plus: selects the substitution model
    # -B 1000         ultrafast bootstrap, 1000 replicates
    # -alrt 1000      SH-aLRT branch test, reported alongside UFBoot
    # -o              outgroup; rooting is stated, not assumed
    #
    # A node is conventionally regarded as well supported when
    # SH-aLRT >= 80 and UFBoot >= 95.
    local og
    og=$(grep '^>' "${fam}.trim.aln" | sed 's/^>//' | grep "^${outgroup}" | head -1)
    if [ -n "$og" ]; then
        echo "  outgroup             : $og"
        "$IQ" -s "${fam}.trim.aln" -m MFP -B 1000 -alrt 1000 \
              -o "$og" -T "${TARO_THREADS:-24}" \
              -pre "${fam}_iq" -quiet -redo
    else
        echo "  outgroup             : none found matching '$outgroup', midpoint rooting"
        "$IQ" -s "${fam}.trim.aln" -m MFP -B 1000 -alrt 1000 \
              -T "${TARO_THREADS:-24}" -pre "${fam}_iq" -quiet -redo
    fi

    echo
    echo "  model selected: $(grep -m1 'Best-fit model' "${fam}_iq.log" | sed 's/.*: //')"
    echo "  tree          : $W/${fam}_iq.treefile"
    echo "  report        : $W/${fam}_iq.iqtree"
}

run_tree psy ambtri
run_tree ccd ambtri

echo
echo "=================================================="
echo "support summary"
echo "=================================================="
python - "$W" <<'PY'
import re, sys
from pathlib import Path
W = Path(sys.argv[1])

for fam in ("psy", "ccd"):
    f = W / f"{fam}_iq.treefile"
    if not f.exists():
        print(f"\n{fam}: no tree"); continue
    # IQ-TREE writes internal labels as "SHaLRT/UFBoot"
    labels = re.findall(r"\)(\d+\.?\d*)/(\d+)", f.read_text())
    if not labels:
        print(f"\n{fam}: no support labels found"); continue
    strong = sum(1 for a, b in labels if float(a) >= 80 and int(b) >= 95)
    print(f"\n{fam}: {len(labels)} internal nodes")
    print(f"  well supported (SH-aLRT >= 80 and UFBoot >= 95): {strong}")
    print(f"  weakly supported                                : {len(labels)-strong}")
PY
