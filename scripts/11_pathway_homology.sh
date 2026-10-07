#!/usr/bin/env bash
# 11_pathway_homology.sh
#
# Pathway gene inventory by direct homology search rather than orthogroup
# membership.
#
# Why. OrthoFinder reported zero taro phytoene synthase genes. Direct search
# found three PSY-like proteins, one of them a full-length reciprocal best hit
# to Arabidopsis AT5G17230 at 78% identity. All three were in
# Orthogroups_UnassignedGenes.tsv: MCL failed to cluster them at all.
#
# 6.6% of genes went unassigned in that run. Any pathway gene could have met
# the same fate, so orthogroup counts systematically undercount. Direct
# homology search does not have this failure mode.
#
# Method: for each Arabidopsis pathway protein, search the taro proteome,
# then search each taro hit back against Arabidopsis. Report reciprocal best
# hits separately from weaker or partial matches.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

D="$TARO_REFS/proteomes"
W="$TARO_WORK/05_orthology/pathway_homology"
MAP="$TARO_WORK/05_orthology/aratha_protein2gene.tsv"
mkdir -p "$W" && cd "$W"

[ -f colesc.dmnd ] || diamond makedb --in "$D/orthofinder_input/colesc.fa" -d colesc --quiet
[ -f aratha.dmnd ] || diamond makedb --in "$D/aratha.primary.faa" -d aratha --quiet

python "$TARO_CODE/scripts/11_pathway_homology.py"
