# Data sources

All public. No Indonesian material or Indonesian-origin data in this pilot,
so no MTA and no BRIN foreign research permit is engaged at this stage.
That changes the moment landrace material enters — flag before Stage 1.

## Assemblies — the three-way comparison

| Tag | Study | Accession | Platform | Contig N50 | BUSCO | Note |
|---|---|---|---|---|---|---|
| `asm2019` | Bellinger et al. 2020, G3 10:2763 | GCA_009445465.1 | — | — | — | American Samoa consortium; **the one Zou et al. mapped to** |
| `asm2021` | Yin et al. 2021, Mol Ecol Resour 21:68 | CNP0001082 (CNGBdb) | PacBio+ONT+Illumina+Hi-C | 400 kb | 85.7% | 'Longxiangyu'; het 0.45% |
| `asm2026` | Sun et al. 2026, Sci Data 13:802 | JBQGWC000000000 | PacBio Revio HiFi + Hi-C | 18.02 Mb | 96.9% | 'Bun Long'; het 1.83%; LAI 16.32 — **current best** |

asm2026 annotation: Figshare doi 10.6084/m9.figshare.29917034
(`Colocasia_esculenta.Genome.V1.gff3` / `.cds` / `.pep`)

CNGBdb can be slow from outside China. If asm2021 annotation is unobtainable,
build models de novo on all three and demote the "vs published annotation"
panel to secondary.

## Long-read transcriptome

| Tag | Study | Accession | Platform | Tissue | Caution |
|---|---|---|---|---|---|
| `ont_corm` | Zou et al. 2025, Biology 14:173 | PRJNA1073178 | ONT PromethION cDNA | corm, 3 stages x 3 reps | **only ~2.75 GB clean across all 9 libraries**; N50 1191–1395 bp |
| `pb_mixed` | Sun et al. 2026 | SRR34972528 | PacBio Revio full-length | — | generated but annotation used short reads — unexploited |

## Short-read (context / junction support)

| Tag | Accession | Note |
|---|---|---|
| `rna_cncb` | CRA028667 | Sun et al. RNA-seq |
| `rna_pigment` | PRJNA639211 | He et al. corm pigmentation, 4 stages — ~4M reads/library, underpowered |
| `rna_ssr` | PRJNA387094 | Wang et al. 2017 de novo |

## Araceae proteomes — orthology outgroups

Same set Sun et al. used for homology-based prediction, so the comparison is
anchored to the published annotation's own assumptions:

- *Amorphophallus konjac* — Gao et al. 2022
- *Zantedeschia elliottiana* — Wang et al. 2023, Sci Data 10:605
- *Pistia stratiotes* — Qian et al. 2022, Mol Ecol Resour 22:2732
- *Spirodela polyrhiza* — An et al. 2019, PNAS 116:18893 (closest relative, ~73.23 MYA)

Plus reference proteins: *Manihot esculenta* MePSY1 (Manes.02G081700),
MePSY2 (Manes.01G124200), AtOR / AtOR-like, IbOr, CmOr.

## Target gene set

```
MEP          DXS DXR MCT CMK MDS HDS HDR
backbone     GGPPS  PSY  PDS Z-ISO ZDS CRTISO
cyclisation  LCYB LCYE
xanthophyll  CYP97A CYP97C BCH ZEP VDE NSX
degradation  CCD1  CCD4  NCED
sink         OR  OR-like  fibrillin  PES
```

**Positive control set** — highly expressed in corm, guarantees a result even
if carotenoid genes are too low to resolve:
`CeAGPL1-4  CeAGPS1-2  CeSBE1-3  CeSS1-3  CeGBSS1`
