# Part 1 — stated method

Part 1 is a pathway inventory: which carotenoid and MEP biosynthesis genes
taro carries, and in what copy number. It motivates Part 2 and does not stand
as a study on its own.

**This document is written after the analysis, not before, and says so.** The
first pass was exploratory: five methods were tried for one question, two were
wrong, and the thresholds were chosen in the moment. The results survived
external validation, but the record was a history of attempts rather than a
method.

What follows is the method as it should have been stated at the outset. The
pipeline implementing it is re-run end to end, so the reported numbers are
produced by this method and not by the path that found them. Anything this
method cannot reproduce is dropped.

The exploratory history is in `LOGBOOK.md` and the superseded scripts are in
`scripts/archive/`. Neither is deleted, because a methods section eventually
has to explain why the pipeline looks the way it does.

---

## The governing rule

> **A copy-number claim requires a gene tree with an outgroup of known
> identity. A presence claim requires only a reciprocal best hit.**

This distinction was absent from the first pass, which is why three findings
ended up at three evidentiary standards without that being stated.

- *Taro carries one PSY* is a copy-number claim. It depends on where the family
  boundary falls, which is what a tree with a known outparalog establishes and
  what a similarity threshold cannot.
- *Taro carries ORANGE and ORANGE-like* is a presence claim. A reciprocal best
  hit at high identity and coverage settles it. A tree would add nothing.

Applying the rule decides which families get the full treatment, rather than
that being decided by which ones turned out to be interesting.

---

## Inputs

Seventeen proteomes: taro plus sixteen others. Eight follow Yin et al. 2021,
who used that set for taro gene-family clustering; cassava was added because it
is the source of the proposed transgene; seven monocots were added to widen the
sampling within the monocots.

One protein per gene throughout, the longest isoform, with the gene identifier
taken from the GFF3 CDS attributes rather than from the protein header. Every
species' gene count is checked against the count NCBI reports, and a mismatch
is an error rather than a note.

**Limitation, stated once and carried everywhere:** no member of the Araceae
has a protein set at NCBI. Taro's closest sequenced relative in this set is
*Zostera marina*, which diverged far earlier than *Spirodela polyrhiza* would
have. Nothing in this analysis dates a duplication or a loss, and no claim is
made that requires doing so.

---

## Defining the pathway gene set

Arabidopsis anchors, one locus identifier per pathway step, cross-checked
against KEGG ko00906 (carotenoid biosynthesis) and ko00900 (terpenoid backbone)
rather than assembled from memory.

The first pass assembled this list by hand and never verified it. Since a
published rice phytoene synthase was missed once already, a gene set that was
never checked for completeness is a plausible place for the same failure.

Anchors are resolved to protein accessions through the Arabidopsis GFF3, not by
gene symbol, because the NCBI gene query returns a gene identifier and no
protein accession, and because symbol lookups failed on this dataset.

---

## Homology

Reciprocal best hit. Forward: each Arabidopsis anchor against the taro
proteome. Reverse: each taro hit against the Arabidopsis proteome. A taro gene
is a one-to-one ortholog when it is a reciprocal best hit covering at least 70%
of the query.

Orthogroup clustering is **not** used. OrthoFinder was run in the first pass
and placed zero taro genes in the phytoene synthase orthogroup, because MCL
failed to cluster them at all — all three taro candidates appeared in
`Orthogroups_UnassignedGenes.tsv`, and 6.6% of all genes went unassigned in
that run. Clustering is not reliable for targeted gene-family questions on this
dataset, and including it would be keeping a step that contributed nothing.

### Thresholds, and a sensitivity analysis instead of a justification

| Parameter | Value |
|---|---|
| DIAMOND e-value | 1e-5 (inventory), 1e-20 (family trees) |
| minimum alignment length for tree inclusion | 150 aa |
| minimum query coverage for a one-to-one call | 70% |

None of these has a principled justification; they are conventional. The first
pass set them in the moment and never tested them, and the 150 aa filter is
what excluded the two PSY fragments.

So the pipeline reports the inventory across a grid — alignment length 100,
150, 200 aa and coverage 50%, 70%, 80% — and any gene whose copy number changes
across that grid is flagged in the output. A count that moves with an arbitrary
parameter belongs in the paper as such.

---

## Gene trees, for copy-number claims only

Families with more than one taro ortholog, or where the boundary is in
question, get the following. Every family gets the same treatment; none is
exempted for being less interesting.

**Two trees, because one cannot do both jobs.**

*Tree A, identity.* The family plus a declared outparalog of known identity.
Establishes where the family boundary falls.

*Tree B, resolution.* The sequences tree A placed inside the family, without
the outparalog, rooted on *Amborella trichopoda*. The outparalog is divergent
enough that including it cost 50% of alignment columns against 45% without it,
and those columns are what resolve close relationships.

| Family | Anchor | Outparalog for tree A |
|---|---|---|
| PSY | AT5G17230 | squalene synthase, AT4G34640 / AT4G34650 |
| CCD4 | AT4G19170 | NCED3, AT3G14440 — same superfamily, mutually anchoring |
| NCED | AT3G14440 | CCD4, AT4G19170 |
| DXS | AT4G15560 | declared at run time, recorded before use |
| HDR | AT4G34350 | declared at run time, recorded before use |

**Pipeline, identical for every family:** FAMSA alignment, trimAl `-automated1`,
IQ-TREE with ModelFinder, 1000 ultrafast bootstrap replicates and 1000 SH-aLRT
replicates. Model is selected per tree by BIC and reported per tree, because
different sequence sets select different models.

Alignment retention after trimming is reported for every tree. A tree built on
less than 40% of its columns is reported with that figure attached.

### Support

A node is well supported at **SH-aLRT ≥ 80 and UFBoot ≥ 95**. A claim resting
on a node below that is reported as an observation and not stated as a result.

The first pass used FastTree's SH-like local support at ≥0.70, which is an
approximation rather than a test. Under proper bootstrap, PSY support fell from
17 of 21 nodes to 9 of 21. That gap is why the standard is fixed here.

### Stopping rule

If two analyses of the same question reach opposite conclusions and neither
clears the support threshold, the question is recorded as not determinable and
is not attempted a third time.

This is not hypothetical. The PSY subgroup assignment gave M2 on ten species
and M1 on seventeen, at 74.8/59 and 75.2/77. More data, opposite answer, both
weak. A third attempt would have been method-shopping.

---

## External validation

Copy numbers are compared against Lisboa et al. 2022, who analysed 351 PSY
genes across 166 species, for every species where they publish a count.

Agreement is reported. **Disagreement is reported as disagreement and not
resolved by changing the method**, which is what happened in the first pass:
the `16_` clade rule was adopted, found to contradict the literature, and the
contradiction was initially attributed to the literature rather than to the
rule. The rule was wrong — it discarded rice LOC_Os09g38320, a published and
functionally characterised phytoene synthase.

---

## What Part 1 reports

| Claim | Kind | Evidence required |
|---|---|---|
| taro carries one PSY | copy number | tree A boundary + external validation |
| taro carries three CCD4 | copy number | tree A boundary + genomic coordinates |
| CCD4 pair is tandem, not a split model | structure | coordinates + tree node support |
| taro carries ORANGE and ORANGE-like | presence | reciprocal best hit |
| pathway inventory | copy number | RBH + sensitivity grid |
| taro PSY subgroup | — | **not determinable; recorded as such** |

---

## What Part 1 does not claim

- the timing of any duplication or loss
- that single-copy PSY is ancestral to the Alismatales — eleven monocots
  sampled and only the two Alismatales at one copy is an **observation about
  counts**, and the tree cannot carry it to an inference
- anything about expression or enzyme activity; copy number is not flux, and
  the asymmetry between one synthesis gene and three degradation genes is a
  hypothesis about carotenoid turnover rather than a measurement of it
