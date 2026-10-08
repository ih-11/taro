# Part 1 — stated method

Part 1 is a pathway inventory: which carotenoid and MEP biosynthesis genes taro
carries, and in what copy number. It motivates Part 2 and does not stand as a
study on its own.

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
`scripts/archive/`. Neither is deleted, because a methods section eventually has
to explain why the pipeline looks the way it does.

---

## The governing rule

> **Copy number = family assignment (phylogeny) + independent genomic locus
> (coordinates). Both limbs, or it is not a copy-number claim.**
>
> **Presence = reciprocal homology with adequate coverage.**
>
> **A detection count is not a copy number.** Reciprocal best hit returns a
> lower bound on family size and is reported as *candidates detected*.

An earlier draft of this document required only a tree for a copy-number claim.
That was wrong, and wrong in a way this project had already demonstrated.

A tree establishes what a sequence **is**. It cannot establish that two proteins
correspond to two genomic copies rather than to one gene split across two gene
models, to two haplotypes of one locus, or to a duplicated assembly region. The
`Ces12496`/`Ces12497` locus is exactly that case: two proteins, both genuinely
PSY-like, one gene.

Conversely, coordinates alone establish that two loci are independent but say
nothing about whether either belongs to the family in question.

Neither limb is sufficient. Both together are.

### Criteria for "independent locus"

Non-overlapping coordinates are not enough. Taro is clonally propagated and
heterozygous — 1.83% in Bun Long, 0.45% in Longxiangyu — so allelic sequence
retained as separate contigs is a live possibility rather than a formality.

A locus counts as independent when:

- it sits on a different anchored chromosome from its putative paralog, **or**
- it sits on the same chromosome with non-overlapping coordinates **and**
  distinct flanking gene content, **and**
- it is not on an unanchored scaffold. Sun et al. anchored 96.86% of the
  assembly to 14 chromosomes, so unanchored sequence is where allelic
  duplication would concentrate.

Flanking gene content is read from the GFF3 and reported, not assumed.

### Why "presence" is a different kind of claim

*Taro carries ORANGE and ORANGE-like* asks whether a convincing taro counterpart
of an Arabidopsis gene exists. A reciprocal best hit at high identity and
coverage answers it. A tree would add nothing, and a coordinate check would
answer a question nobody asked.

Separating the two kinds is what keeps the evidence standard uniform rather than
set by which findings turned out to be interesting.

---

## Inputs

Seventeen proteomes: taro plus sixteen others. Eight follow Yin et al. 2021, who
used that set for taro gene-family clustering; cassava was added because it is
the source of the proposed transgene; seven monocots were added to widen
sampling within the monocots.

One protein per gene throughout, the longest isoform, with the gene identifier
taken from the GFF3 CDS attributes rather than the protein header. Every
species' gene count is checked against the count NCBI reports, and a mismatch is
an error rather than a note.

**Limitation, stated once and carried everywhere:** no member of the Araceae has
a protein set at NCBI. Taro's closest sequenced relative here is *Zostera
marina*, which diverged far earlier than *Spirodela polyrhiza* would have.
Nothing in this analysis dates a duplication or a loss, and no claim is made
that requires doing so.

---

## Step 1 — the pathway gene set

Arabidopsis anchors, one locus identifier per pathway step, cross-checked
against KEGG ko00906 (carotenoid biosynthesis) and ko00900 (terpenoid backbone)
rather than assembled from memory.

The first pass assembled this list by hand and never verified it. Since a
published rice phytoene synthase was missed once already, an unverified gene set
is a plausible place for the same failure.

Anchors resolve to protein accessions through the Arabidopsis GFF3, not by gene
symbol, because the NCBI gene query returns a gene identifier and no protein
accession, and because symbol lookups failed on this dataset.

---

## Step 2 — candidate recovery, and what the sensitivity grid tests

```
Arabidopsis anchor
  → permissive homology search
  → candidate family sequences
  → sensitivity analysis: is the candidate set stably recovered?
  → phylogenetic classification        (identity)
  → genomic locus validation           (independence)
  → copy number
```

The grid tests **candidate recovery**, not copy number. An earlier draft said
"any gene whose copy number changes across the grid is flagged", which confused
a detection parameter with an inference. Copy number is established two steps
later, by phylogeny and coordinates.

| Parameter | Values swept |
|---|---|
| minimum alignment length | 100, 150, 200 aa |
| minimum query coverage | 50%, 70%, 80% |
| DIAMOND e-value | 1e-5 (inventory), 1e-20 (family search) |

None of these has a principled justification; they are conventional. The first
pass set them in the moment and never tested them, and the 150 aa filter is what
excluded the two PSY fragments. A candidate set that is not stable across the
grid is reported as unstable.

Orthogroup clustering is **not** used. OrthoFinder was run in the first pass and
placed zero taro genes in the phytoene synthase orthogroup, because MCL failed
to cluster them at all — all three candidates appeared in
`Orthogroups_UnassignedGenes.tsv`, and 6.6% of genes went unassigned in that
run. Clustering was unsuitable for targeted recovery *on this dataset*; that is
an observation about this analysis, not a general property of the method.
Keeping a step that contributed no evidence would make the pipeline look
thorough while misleading about where the conclusions came from.

---

## Step 3 — gene trees, for family assignment

Families where a copy-number claim is made get the following. Every such family
gets the same treatment; none is exempted for being less interesting.

**Two trees, because one cannot do both jobs.**

*Tree A, identity.* The family plus a declared outparalog of known identity.
Establishes where the family boundary falls.

*Tree B, resolution.* The sequences tree A placed inside the family, without the
outparalog, rooted on *Amborella trichopoda*. The outparalog is divergent enough
that including it cost 50% of alignment columns against 45% without it, and
those columns are what resolve close relationships.

| Family | Anchor | Outparalog for tree A |
|---|---|---|
| PSY | AT5G17230 | squalene synthase, AT4G34640 / AT4G34650 |
| CCD4 | AT4G19170 | **CCD7 (AT2G44990) and CCD8 (AT4G32810)** |
| NCED | AT3G14440 | **CCD7 and CCD8** |
| DXS | AT4G15560 | declared at run time, recorded before use |
| HDR | AT4G34350 | declared at run time, recorded before use |

The CCD4 and NCED outgroup changed from an earlier draft, which used each as the
other's outgroup. That is close to circular: it uses NCED to define the CCD4
boundary while using CCD4 to define NCED's, and works only if the two are
reciprocally monophyletic, which is the thing being tested. CCD7 and CCD8 are
carotenoid cleavage oxygenases outside both subfamilies and serve as a genuine
outgroup.

**Pipeline, identical for every family:** FAMSA alignment, trimAl
`-automated1`, IQ-TREE with ModelFinder, 1000 ultrafast bootstrap replicates and
1000 SH-aLRT replicates. Model selected per tree by BIC and reported per tree,
because different sequence sets select different models.

Alignment retention after trimming is reported for every tree. A tree built on
fewer than 40% of its columns is reported with that figure attached.

### Support

A node is well supported at **SH-aLRT ≥ 80 and UFBoot ≥ 95**. A claim resting on
a node below that is reported as an observation, not stated as a result.

The first pass used FastTree's SH-like local support at ≥0.70, an approximation
rather than a test. Under proper bootstrap, PSY support fell from 17 of 21 nodes
to 9 of 21. That gap is why the standard is fixed here.

### Stopping rule

If two analyses of the same question reach opposite conclusions and neither
clears the support threshold, the question is recorded as not determinable and
is not attempted a third time.

Not hypothetical. The PSY subgroup assignment gave M2 on ten species and M1 on
seventeen, at 74.8/59 and 75.2/77. More data, opposite answer, both weak. A third
attempt would have been method-shopping.

---

## Step 4 — the master table

The deliverable is one machine-readable table, not a figure. Part 2 modifies
exactly these columns, so this is the join between the parts.

`results/tables/part1_inventory.tsv`

| Column | Meaning |
|---|---|
| `family` | pathway step |
| `at_locus` | Arabidopsis anchor |
| `taro_gene` | taro gene identifier |
| `assignment` | family assignment from tree A, or `RBH-only` |
| `protein_aa` | protein length |
| `completeness` | intact / fragmented / unresolved |
| `scaffold`, `start`, `end`, `strand` | genomic coordinates |
| `independent_locus` | yes / no / unresolved, by the criteria above |
| `tree_support` | SH-aLRT/UFBoot at the assigning node, or `n/a` |
| `evidence` | `copy-number verified` / `candidate (RBH)` / `presence` |
| `status` | copy 1, copy 2, unresolved locus, split model, … |

The `evidence` column is what stops a detection count being read as a copy
number. Any figure drawn from this table encodes it; the pathway landscape shows
evidence level alongside the count rather than a bare number for every enzyme,
because a figure that displays 28 numbers identically is making 28
copy-number claims regardless of what the caption says.

---

## External validation

Copy numbers are compared against Lisboa et al. 2022, who analysed 351 PSY genes
across 166 species, for every species where they publish a count.

Agreement is reported. **Disagreement is reported as disagreement and not
resolved by changing the method**, which is what happened in the first pass: the
`16_` clade rule was adopted, found to contradict the literature, and the
contradiction was initially attributed to the literature rather than to the
rule. The rule was wrong — it discarded rice LOC_Os09g38320, a published and
functionally characterised phytoene synthase.

---

## What Part 1 reports

| Claim | Kind | Evidence required |
|---|---|---|
| taro carries one intact PSY locus | copy number | tree A boundary + coordinates + external validation |
| a second PSY-like locus is unresolved | structure | coordinates + protein completeness |
| taro carries three CCD4 | copy number | tree A boundary + independent-locus criteria |
| the CCD4 pair is tandem, not a split model | structure | coordinates + flanking genes + node support |
| taro carries ORANGE and ORANGE-like | presence | reciprocal best hit |
| pathway inventory | **detection** | RBH + candidate sensitivity grid |
| taro PSY subgroup | — | **not determinable; recorded as such** |

---

## What Part 1 does not claim

- the timing of any duplication or loss
- that single-copy PSY is ancestral to the Alismatales — eleven monocots sampled
  and only the two Alismatales at one copy is an **observation about counts**,
  and the tree cannot carry it to an inference
- that RBH candidate counts are family sizes; they are lower bounds
- anything about expression or enzyme activity. Copy number is not flux, and the
  asymmetry between one synthesis gene and three degradation genes is a
  hypothesis about carotenoid turnover rather than a measurement of it

---

## Revision note

This document was revised after external methodological review. Three changes:
the governing rule now requires a genomic-locus limb as well as a phylogenetic
one; the pathway inventory is labelled a detection claim rather than a
copy-number claim; and the sensitivity grid is stated as testing candidate
recovery rather than copy number. A fourth change was made on our own reading:
the CCD4 and NCED outgroup moved from each other to CCD7 and CCD8.
