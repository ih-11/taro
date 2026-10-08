# Archived tables

Output of superseded steps, kept because several were cited before the error
behind them was found. Each entry says what replaced it. The corrections
themselves are in `LOGBOOK.md`.

| file | from | replaced by | why |
|---|---|---|---|
| `pathway_orthogroups.tsv` | `08_find_pathway_orthogroups.py` | nothing | orthogroup clustering is not part of the method. It placed zero taro genes in the PSY orthogroup |
| `tree_summary.txt`, `tree_methods.txt` | `12_gene_trees.sh`, `15_read_trees.py` | `tree_assignments.tsv`, `family_assignment.tsv` | one tree per target rather than per family; support read without separating membership from placement |
| `count_verification.txt` | `16_verify_counts.sh` | `04b_read_trees.py` | its clade rule discarded a published rice PSY |
| `monocots_added.txt` | `20_add_monocots.sh` | `docs/proteome_manifest.tsv` | species list, superseded by the manifest with accessions and checksums |
| `psy_subgroups.txt` | `18_psy_subgroups.sh` | nothing | rooted on a single squalene synthase tip, so every count came out zero and the zeros were printed as a table |
| `subgroup_call.txt`, `subgroup_trees.txt`, `subgroup_verdict.txt` | `19_`, `21_`, `22_` | nothing | the PSY subgroup is recorded as **not determinable**. Two analyses reached opposite conclusions, M2 at 74.8/59 on ten species and M1 at 75.2/77 on seventeen, and neither cleared threshold. By the stopping rule no third was run |

## Also archived during the Part 1 rebuild

| file | replaced by | why |
|---|---|---|
| `kegg_only_loci.tsv` | `kegg_only_loci_raw.tsv` + `kegg_scope_review.tsv` | `02b_` overwrote `02_`'s output in place, so the raw KEGG difference could not be recovered. The two files now have one writer each |
| `pathway_homology.tsv` | `homology_hits_all.tsv`, `homology_candidates.tsv` | coverage computed as alignment length over query length, which counts gap columns and exceeded 100%; no query coordinates, so split models were untestable; sub-threshold hits discarded before writing |
| `scope_exclusions.tsv` | `kegg_scope_review.tsv` | an early scope list, before the KEGG cross-check. The committed scope carries a decision and a reason for all 69 loci |
