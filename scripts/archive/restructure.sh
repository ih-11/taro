#!/usr/bin/env bash
# restructure.sh — archive the exploratory pass, leave a clean skeleton.
#
# Part 1's results survived external validation, but the repository records the
# search rather than the method: 27 scripts, three notebooks and fifteen stale
# figures, several superseded, with no way for a reader to tell which steps
# constitute the pipeline.
#
# Nothing is deleted that the logbook refers to. Superseded scripts move to
# scripts/archive/ with a note on why each was replaced, because the methods
# section has to be able to explain the pipeline's shape.
#
# Figures go, because we agreed figures are displayed in the notebook and not
# committed, and fifteen stale PNGs in git contradict that.
set -euo pipefail
cd ~/Code/taro

echo "== archive superseded scripts =="
mkdir -p scripts/archive
for f in 00_fetch_proteomes.sh 05_fetch_outgroups.sh 08_find_pathway_orthogroups.py \
         07_orthofinder.sh 09_pathway_by_geneid.sh 09_pathway_report.py \
         10_psy_direct_search.sh 12_gene_trees.sh 15_read_trees.py \
         16_verify_counts.sh 17_trees_proper.sh 18_inspect_support.py \
         18_psy_subgroups.sh 19_subgroup_call.py 21_subgroup_resolve.sh \
         22_read_subgroup.py 11_pathway_homology.py 11_pathway_homology.sh \
         13_psy_fragments.sh 14_add_cassava.sh 01_fetch_taro.sh \
         02_primary_isoform.sh 03_probe_outgroups.sh 04_probe_yin_set.sh \
         06_rebuild_primary.sh 20_add_monocots.sh; do
    [ -f "scripts/$f" ] && git mv "scripts/$f" "scripts/archive/$f" 2>/dev/null \
        || { [ -f "scripts/$f" ] && mv "scripts/$f" "scripts/archive/$f"; }
done
ls scripts/archive | wc -l | xargs echo "  archived:"

cat > scripts/archive/README.md <<'EOF'
# archive

The exploratory first pass. Kept because `../../LOGBOOK.md` refers to these by
name and because two of them contain errors worth recording.

The current pipeline is in `../`. These are not part of it and are not run.

## Why each was superseded

| Script | Replaced by | Reason |
|---|---|---|
| `00_fetch_proteomes.sh` | `01_inputs.sh` | first fetch attempt, probe and fetch split awkwardly |
| `01_fetch_taro.sh`, `02_primary_isoform.sh` | `01_inputs.sh` | merged |
| `03_probe_outgroups.sh`, `04_probe_yin_set.sh` | — | one-off probes; findings are in the logbook |
| `05_fetch_outgroups.sh` | `06_rebuild_primary.sh`, then `01_inputs.sh` | **error** — grouped isoforms by description text, merging unrelated genes sharing a description. Arabidopsis came out 14,297 against 27,562 expected |
| `06_rebuild_primary.sh`, `14_add_cassava.sh`, `20_add_monocots.sh` | `01_inputs.sh` | merged; same logic applied to all species at once |
| `07_orthofinder.sh` | — | **contributed nothing.** Placed zero taro genes in the PSY orthogroup; all three candidates were unassigned by MCL. Clustering is unreliable for targeted family questions on this dataset |
| `08_find_pathway_orthogroups.py` | `02_pathway_set.py` | **error** — parsed the NCBI gene JSON wrongly and resolved no symbols |
| `09_pathway_by_geneid.sh`, `09_pathway_report.py` | `02_pathway_set.py` | orthogroup-based; superseded with OrthoFinder |
| `10_psy_direct_search.sh`, `11_pathway_homology.*` | `03_homology.sh` | merged; thresholds now swept rather than fixed |
| `12_gene_trees.sh` | `04_families.sh` | FastTree with SH-like support, no trimming, no model selection |
| `13_psy_fragments.sh` | `05_report.py` | coordinate checks folded into the report |
| `15_read_trees.py`, `18_inspect_support.py` | `05_report.py` | merged |
| `16_verify_counts.sh` | `04_families.sh` | **error** — defined the PSY clade as the largest clade with exactly one Arabidopsis sequence. Arabidopsis lost subgroup E3, so the rule cannot distinguish "outside the family" from "inside a subgroup Arabidopsis lacks". It discarded rice LOC_Os09g38320, a published phytoene synthase. It also never evaluated the root |
| `17_trees_proper.sh`, `21_subgroup_resolve.sh` | `04_families.sh` | correct approach, generalised to every family |
| `18_psy_subgroups.sh` | — | **error** — rooted on a single outgroup tip, so no internal clade could contain all outgroup sequences; every count came out zero |
| `19_subgroup_call.py`, `22_read_subgroup.py` | — | subgroup question recorded as not determinable; see `METHODS_part1.md` |
EOF
echo "  wrote scripts/archive/README.md"

echo
echo "== notebooks =="
for n in part1_figures.ipynb part1_figures_v2.ipynb; do
    [ -f "notebooks/$n" ] && { git rm -q "notebooks/$n" 2>/dev/null || rm -f "notebooks/$n"; echo "  removed $n"; }
done

echo
echo "== figures =="
n=$(ls results/figures/*.png 2>/dev/null | wc -l)
git rm -q results/figures/*.png 2>/dev/null || rm -f results/figures/*.png
echo "  removed $n stale PNGs (figures are displayed, not committed)"
cat > results/figures/.gitkeep <<'EOF'
EOF
grep -q '^results/figures/\*' .gitignore 2>/dev/null || \
    printf '\n# figures are displayed in the notebook, not committed\nresults/figures/*\n!results/figures/.gitkeep\n' >> .gitignore

echo
echo "== tables =="
mkdir -p results/tables/archive
for t in count_verification.txt psy_subgroups.txt subgroup_call.txt \
         subgroup_trees.txt subgroup_verdict.txt tree_summary.txt \
         tree_methods.txt monocots_added.txt pathway_orthogroups.tsv; do
    [ -f "results/tables/$t" ] && { git mv "results/tables/$t" "results/tables/archive/$t" 2>/dev/null \
        || mv "results/tables/$t" "results/tables/archive/$t"; }
done
echo "  moved exploratory output to results/tables/archive/"

echo
echo "== clean pipeline skeleton =="
cat > scripts/README.md <<'EOF'
# scripts

The Part 1 pipeline. Five steps, run in order. `METHODS_part1.md` in `../docs/`
states the method they implement and the thresholds they use.

| Script | Does |
|---|---|
| `01_inputs.sh` | proteomes for 17 species, one protein per gene, counts verified against NCBI |
| `02_pathway_set.py` | Arabidopsis anchors, cross-checked against KEGG ko00906 and ko00900 |
| `03_homology.sh` | reciprocal best hits, with the threshold sensitivity grid |
| `04_families.sh` | two-tree design per family: identity then resolution |
| `05_report.py` | fixed decision rules applied; everything reported, including what fails them |

Figures are in `../notebooks/part1_landscape.ipynb` and are displayed rather
than committed.

The exploratory first pass is in `archive/`, with a note on why each script was
replaced. Two of them contain errors that are part of the record.
EOF
echo "  wrote scripts/README.md"

echo
echo "== result =="
tree -L 2 -I 'archive' . 2>/dev/null || ls -R .
echo
echo "next: the five pipeline scripts, then re-run end to end"
