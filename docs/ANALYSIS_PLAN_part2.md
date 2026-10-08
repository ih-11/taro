# Part 2, pre-registered

Written before any Part 2 data is downloaded. Not to predict what the data will
contain, but to stop the inferential decisions moving once it is known which
answer they produce.

The decisions being frozen are: which loci were nominated, what counts as
informative evidence, what counts as sufficient coverage, what counts as
resolution, what happens when evidence is insufficient, and that ONT is primary
with PacBio as replication and no platform claim.

Supersedes three drafts. The second froze a read-count gate; the third replaced
it after the dataset was looked up. Section 4 keeps both, because a
pre-registration that edits out its own mistakes is not one.

---

## 1. What Part 1 handed over, frozen before any read is seen

Part 1 identified these loci from protein and genome evidence alone, and the
inventory was committed before the transcriptome was downloaded. **That
ordering is the design.** Transcript evidence is now an independent test of
predictions already on record, rather than a search through the transcriptome
for examples that look like validation afterwards.

### Candidate split models: one gene across two annotation records

| locus | models | anchor coverage | genomic span |
|---|---|---|---|
| PSY-like, Superscaffold7 | `Ces12496` + `Ces12497` | 126-238 and 281-427, zero overlap, collinear | 3,883 bp |
| DXS-like, Superscaffold2 | `Ces02306` + `Ces02307` | 482-577 and 585-686 | 11,451 bp |
| DXS-like, Superscaffold2 | `Ces02308` + `Ces02309` | 64-206 and 482-716 | 6,680 bp |

### Candidate over-merged model: one record spanning more than one gene

`Ces14428`, the sole taro GGPPS candidate: 2,232 residues against a family
reference median of 360, 13,255 bp and 11 exon records on Superscaffold8. The
prenyltransferase domain is a GGPS family member at SH-aLRT 100 / UFBoot 100;
1,906 residues, 85% of the protein, match nothing in Arabidopsis.

### Carried as observations

Six NCED genes assigned by both trees with neither clearing support. Two
near-identical HDR genes 138 Mb apart. Four near-complete DXS genes, three in a
170 kb tandem array.

---

## 2. The question

**Primary.** How does long-read sequencing depth affect the ability to
adjudicate uncertain taro gene models?

**Secondary.** Are those conclusions reproduced in an independent PacBio
dataset from another cultivar?

An earlier draft asked whether the answer "depends on the platform" while also
stating that no platform ranking is possible, because ONT is Lipu Taro No. 1
and PacBio is Bun Long. Those two statements contradict each other. Any
ONT-to-PacBio difference is platform plus genotype plus library preparation
plus sampling, and the question has been rewritten so it does not ask for
something the design cannot deliver.

---

## 3. The data, as it actually is

Checked against the publication and the BioProject record rather than assumed.

| | PRJNA1073178 | SRR34972528 |
|---|---|---|
| platform | **Oxford Nanopore PromethION, R9.4** | PacBio |
| library | full-length cDNA, SQK-PCS109, **PCR-barcoded** | — |
| cultivar | **Lipu Taro No. 1** | Bun Long |
| tissue | **corm only**, 30 / 60 / 90 days | — |
| samples | 9 (3 stages × 3 replicates) | — |
| volume | 66 Gbases in the BioProject record; the paper reports 2.75 GB clean, read N50 **1,191–1,395 bp** | — |

Whether "2.75 GB clean" is per sample or in total is not in the main text; the
supplementary tables settle it. Both readings are workable.

**The PacBio reads are from the same cultivar the `asm2026` assembly was built
from.** Cultivar and assembly quality are confounded in any PacBio-to-reference
comparison. Platforms are never pooled: a pooled assembly makes platform,
cultivar and depth inseparable.

---

## 4. Feasibility, in two stages

### The gate that was wrong, kept on the record

An earlier draft proposed proceeding at ≥500,000 full-length reads and pivoting
below 100,000. Two faults. It never said per sample or pooled, so its verdict
on this dataset depends on an ambiguity in its own units. And it gated on total
reads when **the experimental unit is usable coverage at the loci being
adjudicated**: 800,000 reads with none from PSY resolves nothing, and 150,000
resolves a locus that happens to be well expressed. Total read count is QC, not
biology.

### Stage 0, before downloading anything

The corm expansion paper's supplementary tables carry per-gene abundance across
the three stages. Look up every locus in section 1 there. This costs nothing
and is the most informative check available, because **a gene not transcribed
in corm cannot have its model resolved by corm RNA at any depth**, and the
targets are carotenoid pathway genes in a starch storage organ. Record the
result before fetching a read.

### Stage 1, dataset-level technical feasibility — permissive

**≥100,000 classified full-length reads: proceed with the depth-series
experiment. Below that: do not interpret transcriptome-wide depth saturation;
restrict to simulation and descriptive locus evidence.**

No higher threshold is set. The depth series is itself the experiment that
determines how much is enough, so deciding beforehand that some number makes
reconstruction valid pre-empts the thing being measured.

Read count, N50 and read-length distribution, mapping rate and full-length
classification rate are reported as **QC**, never as a biological pass or fail.

### Stage 2, locus-level eligibility — the one that matters

**A locus is evaluable at a given depth only when at least five independent
informative reads cover the region required to distinguish the competing
models.** Three states, and the third is not optional:

| state | meaning |
|---|---|
| **resolved** | informative reads discriminate between the competing models |
| **unresolved despite coverage** | enough informative reads, and they do not discriminate |
| **not evaluable** | fewer than five informative reads |

Collapsing the last two would report absence of evidence as evidence of
absence, which is the specific error section 5 exists to prevent.

**Five is a judgement, not a biological law**, and is stated as one. Unlike a
total read count it is tied to the inferential unit, and it stops a gene-model
claim resting on one lucky molecule.

**"Independent" is operational, not rhetorical.** Reads count as independent
when their alignment termini differ, after duplicate collapse. The library is
PCR-amplified, so five copies of one molecule are one observation.

---

## 5. What counts as evidence for a split model

An earlier draft said: *"A transcript spanning both models proves one gene. A
transcript ending inside the gap proves two."* **The second half is wrong** and
is corrected here.

A read can terminate for reasons that have nothing to do with gene structure:
RNA degradation, incomplete reverse transcription, sequencing truncation,
transcript processing, or simply a molecule shorter than the span. In a library
whose N50 is about 1.3 kb that is the common case, not the exception. Absence
of a bridging read is not positive evidence for two genes.

| conclusion | what it requires |
|---|---|
| **one transcriptional unit** | multiple independent reads bridge the two records with coherent splice structure |
| **separate transcriptional units** | independent transcript populations repeatedly terminate and initiate around the respective models, **with adequate local coverage**, and no bridging molecules |
| **insufficient evidence** | coverage exists and does not distinguish them |

**A bridging call additionally requires at least two independent bridging reads
with distinct termini.** SQK-PCS109 is a PCR library and PCR chimeras produce
precisely the artifact that would falsely bridge two adjacent models. A single
bridging read is the one observation this library can manufacture, so it is not
sufficient on its own.

---

## 6. What read length allows, before any coverage question

The library N50 is 1,191–1,395 bp. Coverage cannot rescue a question whose
answer needs a molecule longer than the library contains.

| question | span a read must cover | against N50 ≈ 1.3 kb |
|---|---|---|
| `Ces12496` + `Ces12497` one transcript? | ~0.8 kb CDS plus UTRs | within reach |
| `Ces02306` + `Ces02307` one transcript? | ~0.6 kb CDS | within reach |
| `Ces02308` + `Ces02309` one transcript? | ~1.1 kb CDS | at the edge |
| `Ces14428` 5′ region and GGPS domain one molecule? | **~6.7 kb CDS** | **far out in the tail** |

**The Ces14428 question is most likely not answerable with this dataset**, on
read length alone and independently of depth. Recording that now, before
downloading, is worth more than discovering it afterwards.

---

## 7. Design

```
depth series:  10, 25, 50, 75, 100% of classified full-length reads
               three fixed-seed subsamples per level, seeds recorded

per locus:     informative reads, then one of
               resolved / unresolved despite coverage / not evaluable

reported:      a resolution curve per locus, descriptive
```

**Stable resolution is declared now: a locus is depth-resolved at the lowest
depth where the same interpretation is obtained in 3/3 subsampling replicates
and at every greater sampled depth.** The reported quantity is that depth. It
answers the primary question directly, where a significance test between two
depths does not.

**A control locus set** carries the curve, because if the carotenoid genes are
quiet in corm a curve fitted on four loci estimates nothing. Twenty loci drawn
from the corm transcriptome before any Part 2 analysis, stratified across the
expression range so the low end is represented, each with a multi-exon model of
comparable structure. The curve is described on the controls; the Part 1 loci
are what it is applied to.

**Reads are not molecules.** Subsampling a PCR library subsamples duplicates
with originals, so a naive curve looks better than the same number of
independent molecules would. Duplicates are collapsed where the barcode
structure permits; where it does not, the curve is labelled as being in reads,
and no absolute molecule requirement is read off it.

---

## 8. Statistics

**No significance test between depths.** An earlier draft proposed McNemar for
a locus resolved or not at two depths. Two reasons it is dropped, and the
second applies at any sample size:

- four hand-selected annotation anomalies are nowhere near the regime where it
  adds inference;
- **nested subsamples are not paired observations.** The larger subsample
  contains the smaller one's reads. The dependence is directional and
  structural, and no amount of extra loci repairs it.

| quantity | treatment |
|---|---|
| locus resolution across depths | descriptive curve, plus replicate consistency |
| a proportion near 1 | **Wilson interval**, since the normal approximation fails at the boundary |
| coverage against a parity line, paired by locus | **Wilcoxon signed-rank** |
| copy number between species | **none**. Tips on a tree are not independent samples |

Replicate subsamples at one depth are technical replicates of one library.

---

## 9. What Part 2 will not establish

- **Anything about a gene not transcribed in corm.** Stage 0 measures this
  before any read is downloaded.
- **The Ces14428 question, most likely**, for the read-length reason in
  section 6.
- **A platform ranking.** ONT is Lipu Taro No. 1 and PacBio is Bun Long.
- **Separate genes from absent bridging reads alone**, per section 5.
- **An absolute molecule requirement**, per section 7.
- **A corrected annotation.** Part 2 can show a model is wrong and show the
  evidence; re-annotation is downstream.
- **Expression level**, or anything about carotenoid content.

---

## 10. Stopping rules, carried forward from Part 1

- Two analyses of the same question reaching opposite conclusions with neither
  clearing threshold: **not determinable**. No third attempt.
- A threshold declared before the data is not moved after it.
- A step that returns nothing exits non-zero rather than printing an empty
  result as a finding.
- A quantity that cannot be checked before the data arrives is not used as a
  gate on whether the data arrives.

---

## Sources

- Biology 14(2):173 — corm expansion study; platform, cultivar, stages, N50
- NCBI BioProject PRJNA1073178 — 66 Gbases, 9 experiments
- Pardo-Palacios et al., Nat Methods 2024 — LRGASP, annotation-free metrics
- Long-read depth benchmarking in iNeurons, bioRxiv 2026
