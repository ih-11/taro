# Part 2 — pre-registered analysis plan

**Status:** draft for agreement. Once committed, changes go in the Deviations
log at the end with a date and a reason, and are never made silently.

**Written before any read is downloaded or any alignment is run.** That is the
point. Part 1 chose its methods after seeing each result, which is why five
methods were tried for one question and why two of them were wrong. This
document exists so that cannot happen again.

---

## 1. Question

Taro has three published genome assemblies of very different quality. Does the
same long-read transcript data produce the same transcript models against each?

This is a methodological question. The carotenoid pathway is the motivating
case, not the subject.

## 2. Hypotheses, and what each predicts

**H1.** Transcript structures reconstructed from identical reads differ
systematically with the reference assembly, and the differences track assembly
contiguity rather than being random.

*Predicts:* cross-assembly concordance of intron chains is materially below the
ceiling set by simulation and by the aligner control; and pairwise concordance
is highest between the two best assemblies.

**H0.** Differences between assemblies are no larger than those produced by
read sampling and aligner choice alone.

*Predicts:* cross-assembly concordance falls inside the band set by the
controls, with no ordering by assembly quality.

**Both outcomes are reportable and both will be reported.** H1 says
"re-annotate before you design"; H0 says "existing taro annotations are more
robust than their provenance suggests, and can be reused." The second is less
exciting and more useful to the field. Committing to publish either is what
makes this a test rather than a demonstration.

### Prior expectation

Zhang et al. (RNA 2019) mapped one rice RNA-seq dataset to three rice genomes.
Gene-level calls were fairly robust: 76% of syntenic orthologs gave the same
differential-expression call across all three. Splicing was not: only 24% of
differentially-used loci were found by both of two genomes, and one reference
yielded 3.4-fold more such loci than another.

Transcript-level analysis is therefore expected to be far more
reference-sensitive than gene-level. If taro behaves like rice, the effect
should be large rather than marginal. Stating this in advance means a small
effect is a real result rather than a disappointment to be explained away.

---

## 3. Data, fixed now

### Assemblies

| Tag | Accession | Platform | Contig N50 | BUSCO | Cultivar |
|---|---|---|---|---|---|
| `asm2019` | GCA_009445465.1 | — | — | — | — |
| `asm2021` | CNP0001082 (CNGBdb) | PacBio + ONT + Illumina + Hi-C | 400 kb | 85.7% | Longxiangyu |
| `asm2026` | JBQGWC000000000 | PacBio Revio HiFi + Hi-C | 18.02 Mb | 96.9% | Bun Long |

### Reads

| Tag | Accession | Platform | Cultivar | Note |
|---|---|---|---|---|
| `ont` | PRJNA1073178 | ONT PromethION cDNA | Lipu Taro No.1 | 9 libraries, ~2.75 GB clean total |
| `pb` | SRR34972528 | PacBio Revio | Bun Long | **same cultivar as asm2026** |

### A confound that must be declared, not discovered later

The PacBio reads come from Bun Long, which is the cultivar `asm2026` was built
from. Those reads will map better to `asm2026` partly because it is the same
plant, not because the assembly is better. Cultivar and assembly quality are
confounded in that comparison, and no analysis can separate them.

The ONT reads are from Lipu Taro No.1, which is none of the three assemblies.
They are therefore equally foreign to all three.

**Decision: `ont` is the primary experiment. `pb` is a secondary replication
with the confound stated wherever it appears.** Results from `pb` that agree
with `ont` strengthen the conclusion; results that disagree are attributed to
the confound unless shown otherwise, and that asymmetry is declared here rather
than chosen afterwards.

Taro is clonally propagated and highly heterozygous (1.83% in Bun Long, 0.45%
in Longxiangyu). Cross-cultivar mapping is therefore a real divergence effect,
analogous to the SNP-density effect Zhang et al. measured between rice
subspecies, and it will be quantified the same way rather than assumed away.

---

## 4. Feasibility gate, evaluated before anything else

~2.75 GB of clean ONT data across nine libraries is thin for isoform work.

**Gate:** count full-length reads after pychopper across all `ont` libraries.

- **≥ 500,000 full-length reads** — proceed with the plan as written.
- **100,000 to 500,000** — proceed, but the primary outcome is restricted to
  loci above a per-locus read threshold fixed at the same time, and the
  restriction is reported.
- **< 100,000** — the plan is abandoned and replaced by the depth-requirement
  study: simulate at a range of depths, report the depth at which transcript
  structures become stably recoverable. That is a complete and useful result,
  not a failure, and committing to it now prevents the data being tortured.

This gate is evaluated once, before any mapping, and the outcome recorded.

---

## 5. Methods, fixed now

### Alignment

```
ONT     minimap2 -ax splice -uf -k14
PacBio  minimap2 -ax splice:hq -uf
```

**Aligner control.** Zhang et al. ran three aligners to show genome choice
dominated aligner choice. Without that control we cannot attribute an effect to
the reference rather than to minimap2 behaving differently on assemblies of
different contiguity.

A second long-read spliced aligner (uLTRA or deSALT; chosen at setup, recorded
here before use) is run on the same reads against the same three assemblies.
The aligner is a nuisance variable, not a treatment.

### Transcript models

IsoQuant, reference-guided, **with the assembly's own annotation withheld**.
Models are built from the reads and the genome sequence only.

This is deliberate. SQANTI3's structural categories — FSM, ISM, NIC, NNC — are
defined relative to a reference annotation. Classifying the same transcript
under three different annotations of three different qualities would measure
annotation quality while appearing to measure transcript stability. Withholding
the annotations removes that confound at the cost of losing the category
vocabulary, which is reintroduced only in Section 6 under a single common
annotation.

StringTie2 `-L` is run as a secondary method, treated as a nuisance variable in
the same way as the aligner.

### Common frame

Coordinates are not comparable across assemblies. Zhang et al. solved this by
restricting analysis to syntenic orthologs present in all three genomes.

`asm2026`'s annotation is lifted to `asm2021` and `asm2019` with Liftoff, and
the three are compared with LiftoffTools (variants, synteny, clusters). The
analysis set is the loci successfully lifted to all three with ≥90% coverage.

LiftoffTools' clusters module additionally reports copy-number change per gene
between assemblies, which answers directly whether the CCD4 and PSY loci from
Part 1 are annotated consistently across the three taro assemblies. That is the
Part 1 question asked across references, and it is how the two parts connect.

---

## 6. Outcomes

### Primary

**Intron-chain concordance.** For each locus in the common frame, the set of
transcript models produced under each assembly, compared as ordered intron
chains.

*Match definition, fixed now:* two transcripts match if every internal splice
junction is identical to the base pair. Terminal exon boundaries are **not**
required to match, because long-read 5′ and 3′ ends are unreliable — ONT and
Iso-Seq both suffer 5′ degradation. Loosening this later would be a deviation.

Reported as the proportion of loci whose model set is identical across all
three assemblies, with a **Wilson** confidence interval. Wilson rather than
normal approximation because the proportion may sit near 1, where the normal
approximation misbehaves.

### Secondary, annotation-free

Following LRGASP Challenge 3, which faced the same absence of trustworthy
ground truth:

- transcripts per locus
- transcript model length distribution
- ORF completeness and coding potential

These compare cleanly across assemblies because none of them references an
annotation.

### Tertiary, within the common annotation only

SQANTI3 structural categories, computed against the lifted `asm2026`
annotation for all three — one annotation, three assemblies. This is the only
form in which the categories are interpretable here, and the alluvial diagram
of category change is drawn only in this form.

### Statistical treatment

- The comparison is **paired** — the same locus under two assemblies — so
  category change is tested with **McNemar**, not chi-square.
- Where an effect size is reported for a paired proportion, it is reported with
  a confidence interval and not only a p-value.
- No test is run on anything with fewer than 20 informative loci; the number is
  reported instead.

---

## 7. The null, and why it is simulation rather than resampling

An earlier draft proposed bootstrap-resampling the reads against one assembly
to get a baseline churn rate. That gives a variance estimate but no ground
truth.

LRGASP uses simulation instead, with Trans-NanoSim for ONT and IsoSeqSim for
PacBio, because simulated reads come from a known transcript set. That yields
sensitivity and precision directly — LRGASP defines sensitivity as true
positives over known transcripts, precision as true positives over mapped
transcripts — rather than only a statement that something changed.

**Design.** Simulate ONT reads from the lifted `asm2026` transcript set at the
observed depth and error profile. Map and reconstruct against all three
assemblies with the identical pipeline. Because the source transcripts are
known, this gives:

1. a **ceiling** — the concordance achievable when the transcripts really are
   identical, which is below 100% because of mapping and reconstruction error
2. per-assembly **sensitivity and precision** against known truth
3. the attribution: how much of the real-data discordance is pipeline noise and
   how much is the reference

Without this, "N% of models changed" has no denominator and a referee asks
"compared to what?" with no answer available.

Simulation is run **before** the real comparison is interpreted, so the ceiling
is not chosen after seeing the result it is used to judge.

---

## 8. Decision rules

Fixed now, applied as written.

**H1 is supported** if, in the `ont` primary analysis, cross-assembly
intron-chain concordance is below the simulated ceiling by more than the
aligner-control band, **and** pairwise concordance is ordered by assembly
contiguity.

**H0 is supported** if cross-assembly concordance falls within the
aligner-control band, with no ordering by assembly quality.

**Inconclusive** if the simulated ceiling is itself so low that the comparison
cannot discriminate. In that case the reported result is the ceiling, framed as
a statement about what current tools can resolve on a 2.3 Gb, 85%-repeat
genome.

**Partial** outcomes — effect present but unordered, or ordered but inside the
control band — are reported as such and not rounded toward either hypothesis.

---

## 9. Confirmatory versus exploratory

**Confirmatory**, covered by this plan: everything in Sections 6 to 8.

**Exploratory**, labelled as such wherever it appears and never presented as
hypothesis-tested: the carotenoid pathway loci specifically, the PSY
`Ces12496`/`Ces12497` locus as a case study, and any per-gene observation.

The PSY locus is the figure most likely to persuade a reader, and it is a
single locus chosen because Part 1 made it interesting. It is an illustration
of the general result, never evidence for it.

---

## 10. What would invalidate this analysis

Stated now so they are not rationalised later.

- `asm2021` unobtainable from CNGBdb. The study drops to two assemblies, which
  weakens the contiguity ordering to a single contrast. Recorded, not hidden.
- Liftoff transfers too few loci to all three to leave a usable common frame.
- The simulated ceiling is indistinguishable from the real concordance.
- The two read sets disagree in direction, which given the cultivar confound
  would mean neither can be interpreted alone.

---

## 11. What this plan does not cover

Part 3, the promoter analysis, is out of scope and will need its own plan. One
correction carried forward: Part 3 was previously described as needing only the
genome. It does not. Empirical transcription start sites require transcript
evidence, and long reads give observed 5′ ends that bound the start region
rather than pinpointing it. Part 3 is comparative cis-regulatory analysis using
long-read-supported 5′ ends, and will be described that way.

---

## 12. Deviations log

Every departure from the above, with date and reason. An empty log is a claim
that nothing changed, so it is left empty only if that is true.

| Date | Section | Change | Reason |
|---|---|---|---|
| | | | |

---

## References

- Zhang et al. 2019, *RNA* 25:669 — choice of reference genome affects
  differential expression and alternative splicing
- Pardo-Palacios et al. 2024, *Nat Methods* — LRGASP consortium
- Pardo-Palacios et al. 2023 — SQANTI3
- Shumate & Salzberg — Liftoff and LiftoffTools
