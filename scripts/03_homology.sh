#!/usr/bin/env bash
# 03_homology.sh — candidate recovery, with the sensitivity grid.
#
# What this step does and does not establish.
#
#   Reciprocal best hit identifies likely taro counterparts of Arabidopsis
#   pathway genes. That is a DETECTION result. It is a lower bound on family
#   size, because a paralog divergent enough that its best Arabidopsis hit is a
#   different gene never appears.
#
#   Copy number is established two steps later, by phylogeny (04_) and genomic
#   coordinates (05_). Nothing here is a copy-number claim.
#
# The grid sweeps alignment length and coverage to ask whether the CANDIDATE SET
# is stably recovered. An earlier draft of the method described the grid as
# testing copy number, which confused a detection parameter with an inference.
#
# None of these thresholds has a principled justification; they are
# conventional. The 150 aa filter is what excluded the two PSY fragments in the
# exploratory pass, which is the reason for sweeping rather than defending them.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
W="$TARO_WORK/03_homology"
mkdir -p "$W" && cd "$W"

[ -f colesc.dmnd ] || diamond makedb --in "$D/orthofinder_input/colesc.fa" -d colesc --quiet
[ -f aratha.dmnd ] || diamond makedb --in "$D/orthofinder_input/aratha.fa" -d aratha --quiet

python3 "$TARO_CODE/scripts/03_homology.py"
