# Data sources

Replaces the ten-species version, which was wrong in most of its particulars
after the reference set grew to seventeen.

## Reference proteomes

Seventeen species, one protein per gene. The authoritative record, with
accession, assembly build, annotation source, protein count and a
sequence-content checksum per proteome, is:

**`docs/proteome_manifest.tsv`**

That file was reconstructed rather than kept. `manifest.tsv` in the store was
intended to be the provenance record and held four rows, for *Spirodela
polyrhiza*, *Pistia stratiotes*, *Amorphophallus konjac* and *Zantedeschia
elliottiana*, all marked "manual" — a list of genomes that were wanted and
never fetched. It contained no row for any of the seventeen proteomes the
inventory is built on, so not one had a recorded accession, assembly version or
annotation release.

Accessions were recovered from the `#!genome-build-accession` and
`#!annotation-source` headers of the nine NCBI GFF3 files before those files
were deleted, cross-checked against the accessions written into the archived
download scripts, and written out with a per-proteome checksum of the sorted
sequence set so content drift becomes detectable. The headers themselves are
kept in `provenance/sources/*_gff_header.txt`.

### Why these seventeen

| role | species | tag |
|---|---|---|
| target | *Colocasia esculenta* | colesc |
| anchor, all gene identities come from here | *Arabidopsis thaliana* | aratha |
| most distant angiosperm, roots the family trees | *Amborella trichopoda* | ambtri |
| Alismatales, taro's own order | *Zostera marina* | zosmar |
| Araceae-adjacent monocots | *Ananas comosus*, *Asparagus officinalis*, *Elaeis guineensis*, *Musa acuminata*, *Phoenix dactylifera* | anacom, aspoff, elagui, musacu, phodac |
| Poales, where PSY copy number is best described | *Oryza sativa*, *Zea mays*, *Sorghum bicolor*, *Setaria italica*, *Brachypodium distachyon* | orysat, zeamay, sorbic, setita, bradis |
| eudicot outgroups | *Manihot esculenta*, *Solanum tuberosum*, *Nelumbo nucifera* | manesc, soltub, nelnuc |

Cassava (*Manihot esculenta*) was added after an earlier clade rule discarded a
published rice PSY, to give the eudicot side more than Arabidopsis alone.

The pipeline reads `refs/proteomes/orthofinder_input/<tag>.fa`, whose headers
are `<tag>|<identifier>`. The directory keeps its name from the abandoned
OrthoFinder step; nothing in the current method uses OrthoFinder.

## Taro genome and annotation

**`Colocasia_esculenta.Genome.V1`** (Sun et al.), from figshare. The API
response recording which figshare article it came from is kept at
`provenance/sources/taro_figshare_files.json`.

| | |
|---|---|
| genes | 28,253 |
| assembly carrying genes | 2,318 Mb |
| mean spacing | 82 kb per gene |
| anchored sequences | 14, named `Superscaffold1` to `Superscaffold14`, 2,261 Mb, 26,986 genes |
| unanchored | 65, named `unanchor*`, 57 Mb, 1,267 genes |
| heterozygosity, as published | 1.83% |

Taro is 2n = 28, so fourteen chromosomes, and the fourteen `Superscaffold`
sequences match that. Anchored status is taken from the naming rather than from
a length threshold: a cumulative-length rule tried first cut between
Superscaffold13 at 113 Mb and Superscaffold12 at 105 Mb, treating a 105 Mb
sequence carrying 1,963 genes as unplaceable.

The 82 kb mean spacing is load-bearing for the method. It is why genomic
distance is not used as evidence for a split gene model: two genes tens of
kilobases apart with nothing annotated between them is the most ordinary
arrangement in this annotation.

The 1.83% heterozygosity is why a gene on an unanchored sequence is reported
but contributes no copy. A haplotig of a gene already counted cannot be
distinguished from a paralog by coordinates alone.

### A second taro annotation exists and disagrees

**GCA_009445465.1** reports 56,238 genes against this annotation's 28,253, a
roughly twofold difference. Nothing in Part 1 uses it. The discrepancy is
recorded because every count in this inventory is a count *in one annotation*,
and because comparing the two is a legitimate Part 2 question in its own right.

## Pathway definition

Hand-curated anchors cross-checked against KEGG `ath00906` (carotenoid
biosynthesis) and `ath00900` (terpenoid backbone biosynthesis), both queried
with the `ath` organism prefix. An earlier attempt queried `ko00906`, received
nothing, and printed the empty result as though it were data.

The frozen scope is `results/tables/kegg_scope_review.tsv`: every KEGG locus
absent from the hand list, with a decision and the reason for it. 54 excluded,
11 recorded as members of a family already represented by a target, 4 excluded
from scope but retained as phylogenetic references. The target table that
results is `results/tables/pathway_anchors.tsv`, 28 targets.

Both files were committed before the taro search ran. After `03_` has seen taro
results a scope change cannot be distinguished from a result-driven one.

## Part 2 read data

Not yet downloaded. `data/longread`, `data/shortread` and `data/derived` are
empty.

| accession | platform | cultivar | role |
|---|---|---|---|
| PRJNA1073178 | ONT | Lipu Taro No.1 | primary experiment |
| SRR34972528 | PacBio | Bun Long | secondary replication, confound declared |

The PacBio reads are from Bun Long, the same cultivar the `asm2026` assembly
was built from, so cultivar and assembly quality are confounded in that
comparison. The ONT data is therefore the primary experiment and the PacBio
comparison is a replication with the confound stated rather than controlled.
See `docs/ANALYSIS_PLAN_part2.md`.

## Tools

In `envs/taro.yml`. The versions that produced Part 1:

| tool | role |
|---|---|
| DIAMOND 2.2.8 | all homology searching |
| FAMSA | protein alignment |
| trimAl | column trimming, `-gt 0.80` |
| IQ-TREE 2 | tree inference, `-m MFP -B 1000 -alrt 1000` |
| Biopython | reading trees |

## What was deleted, and what was kept

`CLEANUP.md` in the store records the removal of 113,955 files, almost all of
one abandoned OrthoFinder run, taking the store from 3.9 GB to about 350 MB.
The summary tables from that run are kept in
`provenance/orthofinder_20261007/` with a README recording that the run used
**nine** species, not the final seventeen, so its clustering failure was
observed on nine and whether it would also fail on seventeen is untested.
