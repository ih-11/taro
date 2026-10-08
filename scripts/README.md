# scripts

The Part 1 pipeline. Five steps, run in order. `../docs/METHODS_part1.md`
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
