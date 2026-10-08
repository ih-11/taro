# Open questions

Questions this project has raised and not answered. Each says what would answer
it. A question here is not a defect; it is a claim that was not made.

## From Part 1

### 1. What is Ces14428?

Taro's only reciprocal best hit to GGPPS11 is a **2,232-residue** gene model in
a family whose reference median across sixteen species is 360. It spans
13,255 bp on Superscaffold8 with 11 exon records. Its prenyltransferase domain
is a GGPS family member at SH-aLRT 100 / UFBoot 100.

**1,906 of its 2,232 residues, 85% of the protein, match nothing in Arabidopsis
at e < 1e-5.** Every one of the fourteen Arabidopsis prenyltransferases that
does match lands on residues 1907-2232, so the matched part is one domain
recognised by one superfamily, not several genes merged.

So taro's GGPPS copy number is withheld. Either this is a lineage-specific
protein with a GGPS domain at its C-terminus, or an annotation that ran one
gene into adjacent sequence.

**What would answer it:** the second taro annotation (GCA_009445465.1) at the
same locus; RNA-seq or long-read coverage across the model, which would show
whether the 5' 1,900 residues are transcribed as one molecule with the domain;
or a protein-domain scan of the unmatched region against a database wider than
Arabidopsis.

### 2. Is the taro NCED family size six?

Both trees group the same six taro genes with the five Arabidopsis NCEDs, and
neither holds the grouping under resampling: UFBoot 59 in tree A and 43 in
tree B, against a threshold of 95 declared before the data. SH-aLRT clears in
tree A at 94.4 and collapses in tree B at 30.5.

Four of the six are on anchored sequences and two are on `unanchor*`, so even a
clean tree would give four copies and two unresolved.

**What would answer it:** more sequence per column. The CCD alignment retains
504 columns across 136 sequences, and the NCED subclade is the part of that
tree with the least signal. A denser monocot sample, or a codon-level alignment
of the NCED clade alone, would be the way in.

### 3. Is the DXS array real or an annotation artifact?

Ces01301, Ces01303 and Ces01305, at 87%, 73% and 91% anchor coverage, sit
inside 170 kb on Superscaffold1 between 173.10 and 173.24 Mb with one to three
annotated genes between each. Ces23591 at 93% sits elsewhere. Separately, four
fragments cluster inside 18 kb on Superscaffold2 and resolve into two
split-model loci.

DXS is a detection target with no declared root group, so no copy number
attaches to it and none is claimed. But a four-member DXS family with a tandem
array would matter for the first committed step of the MEP pathway.

**What would answer it:** the same questions as Ces14428, applied to the array.
Whether the three near-complete models are each transcribed, and whether the
second annotation agrees they are three genes.

### 4. Are the two HDR genes a duplication or a haplotig pair?

Ces00461 and Ces01140 are 78.2% and 78.0% identical to the anchor over 445 and
441 aa, cover the same part of it, and sit 138 Mb apart on Superscaffold1 with
678 genes between them on opposite strands. Both are on an anchored sequence
and they are 138 Mb apart, which argues against haplotigs. HDR is a detection
target, so this is recorded rather than counted.

**What would answer it:** read depth at both loci. A haplotig pair carries half
the coverage of a true two-copy locus.

### 5. Does orthogroup clustering fail on seventeen species?

The OrthoFinder run that justified abandoning orthogroup clustering used
**nine** species. It placed zero taro genes in the phytoene synthase
orthogroup, with 6.6% of genes unassigned overall. Whether the same failure
occurs on the final seventeen-species set is untested.

This matters only if a reviewer asks why clustering was not used, and the honest
answer has to carry the qualifier.

**What would answer it:** re-running `scripts/archive/07_orthofinder.sh` on the
current seventeen proteomes. The result cannot be compared to the old one,
which is why it was not simply regenerated.

### 6. What is the provenance of `aratha_protein2gene.tsv`?

Every step joins on it, and it maps 48,265 Arabidopsis proteins to 27,562
genes. It was derived from an Arabidopsis GFF3 that is gone, and the
`md5sum.txt` that was supposed to record which release is zero bytes. Its own
checksum is now recorded (`72ef769e3a69038d2826322d62dc2b40`) but that pins the
file, not its origin.

**What would answer it:** re-deriving it from a named Araport or TAIR release
and comparing. If the mapping reproduces, the release is identified. If not,
every locus identifier in Part 1 rests on an unrecorded annotation version,
which would need stating.

## Withdrawn claims

Recorded here so they are not quietly reintroduced.

**PSY subgroup assignment.** Two analyses of whether taro's PSY sits in monocot
subgroup M1 or M2 reached opposite conclusions and neither cleared support:
M2 on a 10-species tree at 74.8/59, M1 on 17 species at 75.2/77. Taro's parent
node subtends 22 tips, so it does not group tightly with anything. Recorded as
not determinable; by the stopping rule, no third attempt was made.

**Duplication timing.** Withdrawn with the subgroup assignment, which it
depended on.

**"Single-copy PSY is ancestral to the Alismatales."** Taro and *Zostera* each
carry one PSY while eleven monocots across five orders carry two to four. That
is an observation about counts on a species list, not a phylogenetic inference,
and two tips do not establish an ancestral state.

## For Part 2

### 7. The feasibility gate

I proposed proceeding at **≥ 500,000 full-length ONT reads** and pivoting to a
depth-requirement study below 100,000. That is my judgement, not a standard,
and it has not been agreed. It needs a number before the data is downloaded,
for the same reason the claim kinds were fixed before the search.

### 8. Annotation-relative metrics across assemblies

SQANTI3's structural categories (FSM, ISM, NIC, NNC) are defined relative to a
reference annotation. Comparing them across assemblies with different
annotations measures annotation quality, not transcript stability. Part 2 needs
the annotation-free alternative, which is what LRGASP Challenge 3 was designed
around: transcripts per locus, model length, coding potential.

### 9. Simulation as the null

A bootstrap resample of reads gives a variance estimate and no ground truth.
Trans-NanoSim or IsoSeqSim give actual ground truth and therefore sensitivity
(TP/known) and precision (TP/mapped). The plan calls for simulation; the
parameters have not been set.
