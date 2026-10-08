# Archived scripts

Nothing here runs as part of the method. Each was superseded, and the reason
matters more than the file, because in several cases the superseded script
produced a result that was reported before the error was found. Those are in
`LOGBOOK.md`.

**Do not run anything in this directory.** `patch_schema.sh` in particular would
corrupt the current `03_homology.py` and `04_families.sh`.

## Fetching and preparing references

| script | replaced by | why |
|---|---|---|
| `00_fetch_proteomes.sh` | nothing; one-off | reference proteomes are now in the store, with provenance in `docs/proteome_manifest.tsv` |
| `01_fetch_taro.sh` | nothing; one-off | taro assembly and annotation from figshare |
| `02_primary_isoform.sh` | `06_rebuild_primary.sh` | grouped isoforms by description text, giving Arabidopsis 14,297 proteins against the correct 27,562 |
| `03_probe_outgroups.sh`, `04_probe_yin_set.sh`, `05_fetch_outgroups.sh` | — | exploratory species selection |
| `06_rebuild_primary.sh` | — | rebuilt the primary-isoform sets from GFF3 `protein_id=` and `Dbxref=GeneID:`. Its output is still in use |

## Orthogroup clustering, abandoned

| script | why |
|---|---|
| `07_orthofinder.sh` | OrthoFinder placed **zero** taro genes in the phytoene synthase orthogroup; all three candidates landed in `Orthogroups_UnassignedGenes.tsv`, with 6.6% unassigned overall. That run used **nine** species, not the final seventeen, so the failure is recorded with that qualifier |
| `08_find_pathway_orthogroups.py` | parsed the NCBI gene JSON wrongly, assuming a `transcripts` array that is not there |
| `09_pathway_by_geneid.sh`, `09_pathway_report.py` | superseded with the clustering approach |

Summaries from that run are in `provenance/orthofinder_20261007/` in the store,
with a README recording the nine-species set.

## Targeted search and trees, first attempts

| script | replaced by | why |
|---|---|---|
| `10_psy_direct_search.sh`, `13_psy_fragments.sh` | `03_homology.py` | exploratory. Printed to the terminal and never wrote its tables, which is why the PSY figures quoted in early logbook entries had no file behind them |
| `11_pathway_homology.py`, `11_pathway_homology.sh` | `03_homology.py` | no reverse-search reporting, no query coordinates, coverage computed as alignment length over query length |
| `12_gene_trees.sh`, `17_trees_proper.sh` | `04_families.py` | one tree per target rather than per family, so a CCD4 tree could not tell a CCD4 paralog from a divergent NCED |
| `14_add_cassava.sh`, `20_add_monocots.sh` | — | species additions, folded into the reference set |
| `15_read_trees.py`, `18_inspect_support.py` | `04b_read_trees.py` | read support without separating membership from placement |
| `16_verify_counts.sh` | `04_families.py` | its clade rule discarded a **published** rice PSY. Arabidopsis lost PSY subgroup E3, so an Arabidopsis-anchored rule cannot tell "outside the family" from "inside a subgroup Arabidopsis lacks". It also used `Bio.Phylo.get_path()`, which never evaluates the root |
| `18_psy_subgroups.sh` | `19_subgroup_call.py`, then abandoned | rooted on a **single** squalene synthase tip, so no internal clade contained all eighteen and every count came out zero. The zeros were printed as a table |
| `19_subgroup_call.py`, `21_subgroup_resolve.sh`, `22_read_subgroup.py` | — | the subgroup question is recorded as not determinable in `docs/open_questions.md`. Two analyses reached opposite conclusions and neither cleared support, so by the stopping rule no third was run |

## Replaced during the Part 1 rebuild

| script | why |
|---|---|
| `03_homology_v1.py` | coverage as alignment length over query length, exceeding 100%; a reciprocal best hit could be suppressed from the printed table; no query coordinates; reciprocity tested against the anchor alone |
| `04_families_v1.sh` | tree B removed only the *Arabidopsis* root sequences, leaving every other species' outgroup member in place; the root group was read from the first row of a group; one coverage floor served both reference tips and taro candidates |
| `05_report_v1.py` | "no annotated gene between them" used as evidence for a split model, in a genome with 82 kb mean gene spacing; exon counts keyed on mRNA identifiers and never folded to genes; a cumulative-length rule that classed a 105 Mb chromosome as unplaceable |
| `patch_schema.sh` | a one-off that renamed columns inside `03_` and `04_` by string replacement. It missed one occurrence, which crashed `03_` after it had printed both of its tables. **Running it now would corrupt the current scripts** |
| `restructure.sh` | one-off repository move. Died partway when `git rm` of the figures removed the directory holding `.gitkeep` |
