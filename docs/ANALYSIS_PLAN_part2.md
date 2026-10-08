# Part 2, pre-registered

Written before any Part 2 data is downloaded. The point of writing it now is
that the decisions below cannot then be made by whoever is holding the result.

Supersedes the version drafted before Part 1 finished, which could not name the
loci Part 2 has to resolve because they had not been found yet.

---

## 1. What Part 1 handed over

Part 1 did not only count genes. It found that this annotation is wrong in two
opposite ways at specific, named loci, and **long reads are the instrument that
settles both**. That is the link between a taro carotenoid question and a
long-read methods question, and it is what makes Part 2 worth doing rather than
a generic assembler comparison.

### One gene across two annotation records

| locus | models | evidence |
|---|---|---|
| PSY-like, Superscaffold7 | `Ces12496` + `Ces12497` | cover anchor residues 126-238 and 281-427, zero overlap, collinear, 259 bp apart on one strand with nothing annotated between |
| DXS-like, Superscaffold2 | `Ces02306` + `Ces02307` | 482-577 and 585-686, 54 bp apart |
| DXS-like, Superscaffold2 | `Ces02308` + `Ces02309` | 64-206 and 482-716, 1,541 bp apart |

**A transcript spanning both models proves one gene. A transcript ending inside
the gap proves two.** No amount of further genomic analysis decides this.

### One record spanning more than one gene

`Ces14428`, the sole taro GGPPS candidate: 2,232 residues against a family
reference median of 360, 13,255 bp and 11 exon records on Superscaffold8. Its
prenyltransferase domain is a GGPS family member at SH-aLRT 100 / UFBoot 100,
and **1,906 residues, 85% of the protein, match nothing in Arabidopsis**. Every
one of the fourteen Arabidopsis prenyltransferases that does match lands on
residues 1907-2232, so the matched part is one domain recognised by one
superfamily, not several genes merged.

Whether the 5′ 1,906 residues are transcribed in the same molecule as the
domain is a question only transcript evidence answers.

### Carried as observations, not claims

Six NCED genes assigned by both trees with neither clearing support (UFBoot 59
and 43 against 95). Two near-identical HDR genes 138 Mb apart. Four
near-complete DXS genes, three of them in a 170 kb tandem array. Read depth at
each locus bears on all three.

---

## 2. The question

**How much long-read transcriptome data is needed to adjudicate a gene model,
and does the answer depend on the platform?**

Operationally: assemble the taro transcriptome from long reads, ask what each
of the loci above looks like in the assembly, and measure how that answer
changes with depth.

This is a benchmarking question with a biological deliverable, which is the
shape the lab works in and the shape that serves the collaboration.

---

## 3. Data, and the confound declared up front

| accession | platform | cultivar | role |
|---|---|---|---|
| PRJNA1073178 | ONT | Lipu Taro No.1 | **primary** |
| SRR34972528 | PacBio | Bun Long | secondary, replication |

**The PacBio reads are from Bun Long, the same cultivar the `asm2026` assembly
was built from.** Cultivar and assembly quality are therefore confounded in any
PacBio-to-reference comparison, and nothing in the design removes that. The ONT
data is the primary experiment for that reason, and the PacBio comparison is a
replication with the confound stated rather than controlled.

Pooling the two platforms is not done. A pooled assembly would make platform,
cultivar and depth inseparable.

---

## 4. The feasibility gate — OPEN, NEEDS A NUMBER

> **This is the one decision still outstanding, and it has to be made before
> PRJNA1073178 is downloaded.**

My proposal, which is judgement and not a standard:

| full-length ONT reads | action |
|---|---|
| ≥ 500,000 | proceed with the design below |
| 100,000 – 500,000 | proceed, with the depth series truncated and said so |
| < 100,000 | pivot: the study becomes "what depth would have been required", using simulation alone |

Set the number now, for the same reason the claim kinds were fixed before the
Part 1 search ran. A gate chosen after seeing the read count is not a gate.

---

## 5. Metrics

### Not SQANTI3 structural categories across assemblies

FSM, ISM, NIC and NNC are defined **relative to a reference annotation**.
Comparing them across assemblies with different annotations measures annotation
quality, not transcript stability. Part 1 is itself the demonstration: the two
taro annotations differ about twofold in gene count, so the same transcript
would land in different categories against each.

SQANTI3 is still run, against one fixed annotation, to describe a single
assembly. It is not used to compare two.

### Annotation-free, following LRGASP Challenge 3

For the no-ground-truth case: transcripts per locus, model length distribution,
coding potential. These are properties of the assembly, not of an annotation it
is scored against.

### Simulation as the null

Trans-NanoSim for ONT, IsoSeqSim for PacBio. A bootstrap resample of reads
gives a variance estimate and no ground truth; a simulation gives actual ground
truth, and therefore **sensitivity (TP / known)** and **precision (TP / mapped)**
rather than a spread.

Simulation parameters are set from the real data's own error and length
profiles, and are recorded before the real assembly is scored.

### Cross-assembly annotation comparison

LiftoffTools, all three modules: `variants`, `synteny`, `clusters`.

---

## 6. Design

```
depth series: subsample the ONT reads to 10, 25, 50, 75, 100% of full-length
              reads, three replicate draws per level, fixed seeds recorded

per subsample: assemble, map to the reference, measure
               - the annotation-free metrics above
               - the state of each named locus from section 1
               - sensitivity and precision against the simulated null

per locus:     a transcript spanning both models of a split pair, or not
               a transcript covering Ces14428's 5' region and its domain, or not
```

The reported curve is **locus resolution against depth**, which is the
deliverable: a statement of how much data is needed before a gene model
question can be answered, grounded in gene models that genuinely needed
answering.

---

## 7. Statistics

| comparison | test | why |
|---|---|---|
| a locus resolved or not, same reads, two depths | **McNemar** | paired binary, and the pairs are the same locus |
| a proportion near 1 | **Wilson interval** | the normal approximation fails at the boundary, and several of these will be near 1 |
| coverage against a parity line | **Wilcoxon signed-rank** | paired, no distributional assumption |
| copy number between species | **none** | tips on a tree are not independent samples. Any such test is phylogenetic pseudoreplication |

Replicate subsamples at one depth are technical replicates of one library.
They bound sampling variation and are not independent biological observations,
and no test treats them as such.

---

## 8. What Part 2 will not establish

- **Expression level.** This is assembly and model resolution. Quantification
  across tissues or cultivars is a different experiment.
- **A platform ranking.** The cultivar confound means a PacBio-to-ONT
  difference cannot be attributed to platform.
- **A corrected annotation.** Part 2 can say a model is wrong and show the
  transcript evidence. Producing a re-annotation is downstream of that.
- **Anything about carotenoid content.** No metabolite data is involved, and
  the pathway genes here are a target set, not a phenotype.

---

## 9. Stopping rules, carried forward from Part 1

- Two analyses of the same question reaching opposite conclusions with neither
  clearing threshold: recorded as **not determinable**. No third attempt.
- A threshold declared before the data is not moved after it.
- A step that returns nothing exits non-zero rather than printing an empty
  result as a finding.
