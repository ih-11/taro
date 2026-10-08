# Pipeline

Numbered in run order. Every script prints what it did and writes its full
output to `results/tables/log_<step>.txt`.

Two conventions hold throughout, both learned the hard way and both recorded in
`LOGBOOK.md`.

**A step that gets nothing exits non-zero.** `02_pathway_set.py` once queried
the wrong KEGG identifier, received an empty response, and printed "only in our
list: 28" as though that were a finding. An empty result is a failure, not a
datum.

**A cache is validated, never assumed present.** `04_families.sh` used to test
`[ ! -f all.dmnd ]`, and a zero-byte database sat next to a 201 MB FASTA built
from an earlier species set. Existence is not validity.

## Inputs

| script | does | writes |
|---|---|---|
| `01_inputs.sh` | verifies the seventeen reference proteomes are present and correct, and fetches only what is missing | the store, `log_01_inputs.txt` |

It verifies before fetching, because everything was downloaded during the
exploratory pass and re-downloading 454,114 proteins to prove the pipeline runs
would take an hour and change nothing.

It also carries the rule the whole inventory depends on: **one protein per
gene, the longest isoform, with the gene identity taken from the GFF3 CDS
attributes** (`protein_id=` paired with `Dbxref=GeneID:`) and never from the
protein header. An earlier attempt grouped isoforms by description text, which
merged unrelated genes sharing a description — "F-box protein",
"uncharacterized protein" — and produced 14,297 Arabidopsis genes against the
27,562 expected. Every count downstream would have been wrong by a factor that
varied per species.

## Scope, before any taro data

| script | does | writes |
|---|---|---|
| `02_pathway_set.py` | hand-curated anchors, cross-checked against KEGG `ath00906` and `ath00900` | `pathway_anchors.tsv`, `kegg_only_loci_raw.tsv` |
| `02b_annotate_kegg.py` | gives every KEGG-only locus its name and a bucket | `kegg_scope_review.tsv` |
| `02c_propose_scope.py` | proposes a decision and a reason for each; leaves the genuinely unclear as REVIEW | `kegg_scope_review.tsv` |
| `02d_finalise_scope.py` | resolves the remaining REVIEW rows; GGR and GPS1 are judgement, the rest mechanical | `kegg_scope_review.tsv` |
| `02e_group_families.py` | de-duplicates targets, records which share a gene family | `pathway_anchors.tsv` |

`kegg_scope_review.tsv` is committed at this point. After `03_` has seen taro
results, a scope change cannot be told apart from a result-driven one.

## Evidence

| script | does | writes |
|---|---|---|
| `03_preflight.sh` | checks the scripts downstream read the anchors schema as it now stands, and that every input exists | — |
| `03_homology.sh` → `03_homology.py` | reciprocal best hit both directions, the nine-cell threshold grid, split-model candidates from query coordinates | `homology_hits_all.tsv`, `homology_candidates.tsv`, `homology_sensitivity.tsv` |
| `04_families.sh` → `04_families.py` | one tree per family group: tree A with the root group, tree B without | `work/04_families/` in the store, `jobs.json` |
| `04b_read_trees.py` | reads the trees. Builds nothing | `tree_assignments.tsv`, `family_assignment.tsv` |
| `05_report.py` | genomic coordinates, locus independence, split and over-merged models | `part1_inventory.tsv` |
| `06_readme_block.py` | regenerates the README findings section from the master table | `README.md` between its markers |

`03_` writes **every** forward hit at e < 1e-5 whether or not it passes any
threshold, and never suppresses a reciprocal best hit from the printed table.
The two PSY fragments fail both default thresholds, and an earlier print filter
hid them while the file held them all along.

`04_` and `04b_` are separate because a tree takes half an hour to infer and
should not be rebuilt every time the reading changes, and because a reading
step that cannot reach back and alter how the tree was made is harder to fool.

## Maintenance

| script | does |
|---|---|
| `cleanup_store.sh` | plans by default, deletes only with `--apply`. Rescues the audit trail into `provenance/` and verifies it before removing anything |
| `extract_provenance.sh` | recovers proteome accessions from NCBI GFF3 headers. Reads and copies only |
| `polish_repo.sh` | archives superseded files by `git mv`. Plans by default, `--apply` does it |

## What a step does not establish

Each script's docstring states this, and the statements compound:

- `03_` gives **detection**, a lower bound on family size
- `04b_` gives **family assignment**, not a count
- `05_` gives the **locus**, which is what turns an assignment into a copy

A copy number appears in the master table only where an assignment clearing
SH-aLRT ≥ 80 and UFBoot ≥ 95 meets an independent locus on one of the fourteen
anchored sequences.
