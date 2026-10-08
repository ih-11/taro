#!/usr/bin/env bash
# 04_families.sh — family trees for the copy-number targets.
#
# Thin wrapper. The logic is in 04_families.py, for the same reason 03_ moved
# there: the previous version carried its tree reading inside bash heredocs, and
# a schema change had to be applied to it by string replacement, which is how
# the raw.get(a["gene"]) crash happened.
#
# WHAT CHANGED FROM THE PREVIOUS VERSION
#
# 1. The database cache is validated, not merely tested for existence.
#
#    The old guard was 'if [ ! -f all.dmnd ]'. A zero-byte all.dmnd was sitting
#    in work/04_families/ next to a 201 MB all.faa built from the NINE-species
#    set, so the next run would have found the file present, skipped the
#    rebuild, and handed DIAMOND an empty index while printing a species count
#    taken from the variable rather than the database. 04_families.py writes a
#    sidecar recording the species list and sequence count and rebuilds when
#    anything disagrees.
#
# 2. Trees are built per FAMILY GROUP, not per target.
#
#    CCD1, CCD4 and the NCEDs share one pool of twelve taro genes, which
#    reciprocal best hit partitions 1 / 2 / 6. That partition is a result. A
#    tree containing only one target's candidates has assumed it, and cannot
#    distinguish a CCD4 paralog from a divergent NCED. So the CCD group gets one
#    tree. GGPPS and PSY are ungrouped and keep a tree each.
#
# 3. The root group is tested rather than trusted.
#
#    Arabidopsis has nine genes in the carotenoid cleavage dioxygenase family
#    and every one of them is inside the family being measured, so there is no
#    outparalog for it. The CCD tree is rooted on the strigolactone-branch
#    members CCD7 and CCD8, which is a sister clade inside the family and not an
#    outgroup outside it. If those two do not come out monophyletic, rooting on
#    the pair is the same error as rooting on a single SQS tip, which is what
#    made 18_psy_subgroups.sh report every count as zero. The script tests it
#    and falls back to midpoint rooting, reported as such, with the copy-number
#    claim weakened accordingly rather than quietly kept.
#
# 4. Which taro sequences enter a tree is declared, not tuned.
#
#    A 33 aa reciprocal hit exists in this data (Ces16580 against ZDS, 6%
#    coverage) and so do several 78 to 143 aa fragments. Putting a 113 aa piece
#    of a 437 aa protein into an alignment contributes mostly gaps and removes
#    the columns the close relatives need, which is the same trade-off that put
#    the divergent outparalog in its own tree. The floor is the most permissive
#    cell of the sensitivity grid already declared in 03_, alignment >= 100 aa
#    and query coverage >= 50%. Nothing is chosen after looking at the trees.
#
#    Sequences below that floor are not discarded from the analysis. They are
#    reported as unplaceable by phylogeny and handed to 05_, where genomic
#    coordinates decide what they are. The two PSY fragments are the case in
#    point: 26% and 34% coverage of non-overlapping stretches of the anchor,
#    which a tree cannot place but adjacent coordinates on one strand can.
set -euo pipefail
: "${TARO_CODE:?source envs/activate.sh first}"
: "${TARO_REFS:?}"
: "${TARO_WORK:?}"

for tool in diamond famsa trimal iqtree2; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        if [ "$tool" = "iqtree2" ] && command -v iqtree >/dev/null 2>&1; then
            continue
        fi
        echo "missing: $tool" >&2
        echo "the taro env should carry diamond, famsa, trimal and iqtree2" >&2
        exit 1
    fi
done

python3 "$TARO_CODE/scripts/04_families.py" "$@"
