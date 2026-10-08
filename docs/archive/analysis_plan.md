# Analysis plan

Metrics are fixed before results are seen. Do not revise afterwards.

## 01_qc
Read counts, length distributions, per-library yield, both long-read sets.
**Decision gate:** if `ont_corm` yields < ~100k full-length reads total,
Part 2 becomes Spine B and the depth-requirement analysis is the deliverable.

## 02_align
minimap2 vs each of asm2019 / asm2021 / asm2026 — six runs.
```
ONT     -ax splice -uf -k14
PacBio  -ax splice:hq -uf
```
Record wall time, peak RSS, primary/secondary/supplementary counts.
Chromosome naming differs across assemblies — normalise at ingest or joins
will silently break.

## 03_isoform
Model sets per (reads x assembly) via IsoQuant (or StringTie `-L`).
One GFF per combination.

## 04_compare — the core result
Coordinates are not comparable across assemblies. Two independent routes:

- **R1 sequence** — `gffread -w` per model set, then all-vs-all / clustering.
  Asks: did I recover the same molecule?
- **R2 projection** — Liftoff onto asm2026, compare splice-junction sets.

Agreement between R1 and R2 is the confidence check; they fail differently.

Metrics, fixed in advance:
- isoforms per locus
- splice-junction concordance (Jaccard)
- 5' end position shift (bp)
- ORF intact / truncated / absent
- locus split or fused between assemblies

## 05_orthology — Part 1
OrthoFinder: taro proteome + 4 Araceae + reference PSY/OR set.
Deliverable: copy number per pathway step, with the PSY paralogy question
stated explicitly (cassava has PSY1/PSY2; saffron has four).

## 06_promoter — Part 3
Empirical TSS from long-read 5' ends — reuse
`Longread_pipeline/Annotation/annotation_for_longread_mod_v1.4.0_hs.py`,
which already emits TSS / CPS / Pair tables.
Extract upstream regions from asm2026. FIMO (MEME .sif in
`/mnt/f/RA/Containers/MEME`) against JASPAR plants plus the Carvita25-specific
motifs from the dissertation: AT-1, ACGTTBOX, SGBFGMGMAUX28.
Background model from taro intergenic sequence.

**State in the paper:** ONT/Iso-Seq 5' ends suffer degradation; motif
prediction has a high false-positive rate. This is hypothesis generation
grounded in empirical TSS — not functional validation. Still a step up from
a guessed upstream window.

## Figures
- **F1** mapping + recovery stats across three assemblies
- **F2** UpSet — isoform concordance (all three / two / one)
- **F3** browser tracks — same reads, three assemblies, one locus *(the memorable one)*
- **F4** pathway table, per-locus verdicts
