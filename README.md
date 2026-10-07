# taro

Carotenoid biosynthesis pathway in taro (*Colocasia esculenta*): gene
inventory, copy number, and the reliability of the published annotation.

This repository supports a collaboration on biofortification of Indonesian
taro landraces, where the plan is to raise beta-carotene content using the
MePSY1 promoter and gene from cassava. Before a construct is designed, it is
worth knowing what taro's own carotenoid genes look like, and whether the
annotation describing them can be trusted.

## The question

Taro has three published genome assemblies of very different quality, and the
same transcript data mapped against different assemblies can produce different
gene models. So:

> How much of taro's transcript level annotation reflects biology, and how much
> is an artifact of which reference assembly the reads were mapped to?

There is already a concrete reason to ask. Two published taro annotations
disagree about the number of protein coding genes by roughly a factor of two,
56,238 against 28,253.

## Workflow

```mermaid
flowchart TD
    A1["Taro proteome<br/>28,253 genes"] --> A3["OrthoFinder"]
    A2["9 outgroup proteomes<br/>incl. cassava"] --> A3
    A3 --> A4["Pathway gene set"]
    A4 --> A5["Gene trees<br/>PSY and CCD/NCED"]
    A5 --> A6["Copy number<br/>per species"]

    B1["ONT corm reads"] --> B3["minimap2<br/>3 assemblies"]
    B2["PacBio reads"] --> B3
    B3 --> B4["Isoform models"]
    B4 --> B5["Reference<br/>stability"]

    B4 --> C1["Empirical TSS"]
    C1 --> C2["Upstream regions"]
    C2 --> C3["FIMO motif scan"]

    A4 -.-> C3
```

**Part 1** inventories the pathway and counts copies. **Part 2** asks whether
transcript models survive a change of reference assembly. **Part 3** locates
where taro's own carotenogenic genes start transcribing, and whether the
cassava promoter has a plausible home there.

Parts 1 and 3 need only the genome, so they survive even if the public long
read data proves too shallow for Part 2. If Part 2 comes back empty, that is
itself a usable result: a quantitative statement of the sequencing depth the
question would require.

## Why Part 1 came first

It could not fail. Orthology on a proteome with 96.9% BUSCO completeness
produces an answer regardless of how the rest goes, and it needs one proteome
plus outgroups rather than seven gigabytes of assemblies.

It also answers the question that most affects construct design: how many
phytoene synthase genes taro carries. If taro had several, then "add a cassava
PSY" would mean something different from what the plan assumes, because
endogenous paralogs may already differ in where and when they are expressed.

## Findings so far

Part 1 is complete except for one verification step. Everything below comes
from public data.

**Taro has one full-length phytoene synthase.** Ces24605 on Superscaffold13,
428 aa against the 437 aa Arabidopsis reference, 78% identity, reciprocal best
hit. A second PSY-like locus on Superscaffold7 is split across two gene models,
Ces12496 and Ces12497, lying 259 bp apart on the same strand with nothing
annotated between them. Combined they recover only 61% of a full length
protein. Whether this is a real second copy or a degenerate locus cannot be
decided from the published annotation, which is itself an instance of the
problem this project set out to document.

**Taro has three CCD4 carotenoid cleavage dioxygenases.** Ces13251, Ces13250
and Ces03723, confirmed by gene tree against Arabidopsis anchors for both CCD4
and NCED3, with every gene roughly twice as far from the wrong anchor as from
the right one. CCD4 degrades carotenoids, and in potato and peach CCD4
variation is a principal determinant of flesh colour.

**Taro has both ORANGE chaperones.** Ces06664 and Ces26954, which control
phytoene synthase stability post-translationally. This matters because OR acts
on PSY protein, not PSY transcript, so transcript level engineering does not
address it.

Taken together: one synthesis gene against three degradation genes. Raising
flux through a single PSY may not produce the expected accumulation if the
degradation step is amplified. This bears directly on construct design and was
not anticipated when the project began.

**Still to verify.** Per-species PSY copy number is currently taken from
homology search hit counts, not from tree topology. Some hits may belong to the
wider squalene and phytoene synthase superfamily. The numbers should not be
quoted until confirmed.

Full reasoning, including two methodological errors made and corrected along
the way, is in `LOGBOOK.md`.

## Repository layout

Code and data are kept apart. This repository holds what is small, text based
and worth version controlling. Large files live on a separate drive and are
never committed, because all of them are either downloadable from public
archives or reproducible by running the scripts here.

```
~/Code/taro/                           this repository
├── envs/
│   ├── taro.yml                       conda environment specification
│   └── activate.sh                    sets paths, TMPDIR, thread count
├── scripts/                           numbered pipeline steps
├── notebooks/                         exploration and figures
├── docs/                              data sources, analysis plan, questions
├── results/                           tables and figures
├── LOGBOOK.md                         what was done, when, and why
└── store ->  /mnt/f/RA/Downstream/Project5_TARO

/mnt/f/RA/Downstream/Project5_TARO/    large files, not tracked
├── refs/        assemblies, annotations, proteomes
├── data/        downloaded sequencing reads
├── work/        intermediates, one directory per step
└── logs/        run logs
```

`store` is a gitignored symlink. It lets you reach the data side from inside
the repository without the data entering version control.

## Setup

```bash
conda env create -f envs/taro.yml      # once
source envs/activate.sh                # at the start of every session
```

`activate.sh` sets the paths the scripts rely on:

| Variable | Points to |
|---|---|
| `TARO_CODE` | this repository |
| `TARO_REFS` | reference sequences |
| `TARO_DATA` | downloaded reads |
| `TARO_WORK` | analysis intermediates |
| `TARO_LOGS` | run logs |

It also pins `TMPDIR` to a scratch directory inside the repository, because
`/tmp` on this machine is tmpfs, meaning RAM backed. A large sort writing there
would consume memory rather than disk and fail in a way that is awkward to
diagnose.

## Scripts

Run in order. See `scripts/README.md` for which are superseded and why.

| Script | Does |
|---|---|
| `01_fetch_taro.sh` | taro proteome and GFF3 from Figshare |
| `02_primary_isoform.sh` | 34,340 transcripts to 28,253 genes |
| `03_probe_outgroups.sh` | test which Araceae proteomes NCBI holds |
| `04_probe_yin_set.sh` | test the substitute species set |
| `06_rebuild_primary.sh` | outgroup proteomes, one protein per gene |
| `07_orthofinder.sh` | orthogroups across nine species |
| `09_pathway_by_geneid.sh` | pathway orthogroups via Arabidopsis locus |
| `10_psy_direct_search.sh` | PSY by direct homology |
| `11_pathway_homology.sh` | pathway inventory by reciprocal best hit |
| `12_gene_trees.sh` | PSY and CCD/NCED family trees |
| `13_psy_fragments.sh` | coordinates of the PSY-like loci |
| `14_add_cassava.sh` | add the transgene source species |
| `15_read_trees.py` | copy number and clade membership from trees |

## Species set

Nine outgroups plus taro. Eight come from Yin et al. 2021, who used them for
taro gene family clustering, which anchors the comparison to a published taro
analysis. Cassava was added separately because it is the source of the
transgene.

| Species | Role |
|---|---|
| *Arabidopsis thaliana* | functional reference for PSY and OR |
| *Manihot esculenta* | cassava, source of the proposed transgene |
| *Oryza sativa*, *Zea mays* | monocot references |
| *Musa acuminata* | monocot, starch storage |
| *Solanum tuberosum* | tuber forming eudicot, carotenoid literature |
| *Nelumbo nucifera* | basal eudicot |
| *Amborella trichopoda* | basal angiosperm, tree root |
| *Zostera marina* | the only other Alismatid available |

**Limitation.** No member of the Araceae has a protein set at NCBI, so taro's
closest sequenced relative, *Spirodela polyrhiza* at roughly 73 million years
divergence, is absent. Duplication timing therefore cannot be dated more
precisely than "after the Alismatales split". This matters specifically because
*Zostera*, the only other Alismatid here, also carries a single PSY, so the low
count may be an Alismatales characteristic rather than something particular to
taro. *Spirodela* would settle that.

## Data and permissions

Every input is public. Nothing in this pilot uses Indonesian plant material or
Indonesian derived data, so no material transfer agreement and no BRIN foreign
research permit is engaged at this stage. That changes as soon as landrace
material enters the project, and should be raised before rather than after.

Accession records are mirrored in a companion repository, `seqfetch`, under its
`reference/`, `nanopore/` and `pacbio/` modules.

## Reading

`docs/data_sources.md` documents the three taro assemblies, the two long read
datasets, and the reasoning behind the species choice.
`docs/open_questions.md` holds the points to raise with the collaborator.
`LOGBOOK.md` is the running record.