# scripts

Run in order. Each reads `$TARO_*` paths set by `envs/activate.sh`.

Some scripts are superseded but kept, because `LOGBOOK.md` refers to them and
the errors they contain are part of the record.

## Current

| Script | Purpose |
|---|---|
| `01_fetch_taro.sh` | taro proteome and GFF3 from Figshare |
| `02_primary_isoform.sh` | 34,340 transcripts to 28,253 genes, longest isoform per gene |
| `03_probe_outgroups.sh` | test which Araceae proteomes NCBI holds (none, as it turned out) |
| `04_probe_yin_set.sh` | test the Yin et al. 2021 substitute species set |
| `06_rebuild_primary.sh` | outgroup proteomes, one protein per gene, via GFF3 gene IDs |
| `07_orthofinder.sh` | orthogroups across all nine species |
| `09_pathway_by_geneid.sh` | pathway orthogroups, joined on Arabidopsis locus |
| `10_psy_direct_search.sh` | PSY by direct homology, after OrthoFinder reported zero |
| `11_pathway_homology.sh` | full pathway inventory by reciprocal best hit |

## Superseded

| Script | Why kept |
|---|---|
| `00_fetch_proteomes.sh` | first fetch attempt, replaced by 01 and 06 |
| `05_fetch_outgroups.sh` | grouped isoforms by description text, merging unrelated genes |
| `08_find_pathway_orthogroups.py` | parsed the NCBI gene JSON wrongly, resolved nothing |

## Next

| Script | Purpose |
|---|---|
| `12_psy_tree.sh` | phytoene synthase gene tree across all nine species |
| `13_ccd_tree.sh` | separate true CCD4 from NCED in the three-copy result |
