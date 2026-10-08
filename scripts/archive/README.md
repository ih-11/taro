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
