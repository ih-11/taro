# Part 1 methods: the carotenoid pathway gene inventory of *Colocasia esculenta*

This is a complete rewrite. The method changed materially during the analysis,
and an earlier version of this file described a design that several steps did
not implement. Where a decision was made after seeing data, it says so.

---

## 1. The question

How many copies of each carotenoid pathway gene does the taro genome carry, and
which pathway steps are present at all?

The question matters for a biofortification target because the dose of a
pathway enzyme depends on how many loci encode it, and because a cleavage
enzyme present in extra copies is a candidate sink for the product a
biofortification programme is trying to accumulate.

## 2. The governing rule

> **Copy number = family assignment (phylogeny) + independent genomic locus
> (coordinates). Presence = reciprocal homology with adequate coverage.
> A detection count is not a copy number.**

Each clause does separate work.

**Family assignment** is a statement about a clade. Reciprocal best hit says
only that a gene's closest relative in a reference proteome lies inside a
declared family, which is a property of a pairwise search. A tree is required
to say the gene falls inside the clade the family's members define.

**Independent genomic locus** is required because two sequences assigned to one
family can still be one gene counted twice. Two mechanisms, both present in
this data: a gene broken across two annotation records, and, in an assembly
reported at 1.83% heterozygosity, the same gene appearing twice as two
haplotigs.

**A detection count is not a copy number** because reciprocal best hit is a
lower bound. A paralog divergent enough that its closest relative in the
reference lies outside the declared family never appears in the result.

## 3. Claim kinds, fixed before any search

Every target was assigned one of three claim kinds before taro was searched,
and the assignment was frozen in a committed file. The purpose is to stop the
evidence standard from being set by which findings turn out to be interesting.

| claim kind | targets | what is established | how |
|---|---|---|---|
| copy-number | 4 | a number of loci | homology, then phylogeny, then coordinates |
| detection | 21 | a lower bound on family size | homology only |
| presence | 3 | the gene exists | reciprocal homology with coverage |

A copy-number target additionally requires a declared group to root its tree
on. Two targets originally marked copy-number, DXS and HDR, had no such group
and were moved to detection before the search rather than being left in a state
where the machine-readable method and the narrative disagreed.

## 4. Scope

The target list began as a hand-curated set of 28 anchors and was cross-checked
against KEGG `ath00906` (carotenoid biosynthesis) and `ath00900` (terpenoid
backbone biosynthesis). KEGG returned 93 Arabidopsis loci across the two
pathways; 23 were already in the hand list, 70 were not, and 5 hand-list
entries were outside both KEGG pathways.

Each of the 70 was given a proposed decision and a reason, and the resulting
file was committed before the taro search ran. The boundary applied:

> **In**: carotenoid pathway enzymes, the plastidial MEP chain supplying them,
> and the accumulation regulators (ORANGE, ORANGE-like, fibrillin).
> **Out**: everything else, including pathways that compete for the same
> precursor.

Dispositions: 54 excluded, 11 recorded as members of a family already
represented by a target, 4 excluded from scope but retained as phylogenetic
references.

Two decisions were judgement rather than mechanism and are recorded as such.
**GGR** (geranylgeranyl reductase, AT4G38460) reduces GGPP to phytyl-PP and is
the main competing sink for the pool PSY draws on. It is excluded, because GGPP
also feeds gibberellins and protein geranylgeranylation, and a scope admitting
one competing sink admits all of them; it belongs in the Part 2 flux
discussion. **GPS1** (AT2G34630) is excluded on a simpler ground: it makes C10
GPP, not the C20 GGPP that PSY condenses.

Two scope changes were made after the taro search and both are additive
bookkeeping rather than changes to which loci are in scope.

**IDI was added** during the KEGG cross-check, before the search, because
skipping IPP/DMAPP isomerisation leaves a biochemical gap between HDR output
and prenyl-diphosphate synthesis.

**BCH was de-duplicated after the search.** The hand list carried both BCH1
(AT4G25700) and BCH2 (AT5G52570), the two Arabidopsis paralogs of one family,
so the target table held two targets for one family. The consequence was a
phantom result: Ces06731 is recovered by both, its reverse best hit is BCH1, so
the second target recorded a stable zero across all nine threshold
combinations, which would have read as a carotenoid hydroxylase family absent
from taro. Targets went from 29 to 28. The set of Arabidopsis loci in scope is
identical before and after.

## 5. Reference data

Seventeen proteomes, one protein per gene. Arabidopsis as the anchor species;
Amborella as the most distant angiosperm; cassava, potato and lotus as
eudicots; and twelve monocots spanning Alismatales, Zingiberales, Arecales,
Asparagales, Poales and Bromeliales.

The provenance of these proteomes had to be reconstructed. The file intended to
record it, `manifest.tsv`, contained four rows for species not in the project
(*Spirodela polyrhiza*, *Pistia stratiotes*, *Amorphophallus konjac*,
*Zantedeschia elliottiana*), all marked "manual": a list of genomes that were
wanted and never fetched. No row existed for any of the seventeen proteomes the
inventory is built on. Accessions were recovered from the NCBI GFF3 headers
before those files were deleted, and are now in `docs/proteome_manifest.tsv`
together with a sequence checksum per proteome so content drift is detectable.

Taro itself: `Colocasia_esculenta.Genome.V1`, 28,253 genes across 2,318 Mb.
Fourteen sequences named `Superscaffold1` to `Superscaffold14` carry 2,261 Mb
and 26,986 genes; 65 named `unanchor*` carry 57 Mb and 1,267 genes. Taro is
2n = 28, so fourteen chromosomes, and the naming matches. Mean spacing is
**82 kb per gene**, which matters in section 8.

A second taro annotation exists (GCA_009445465.1, 56,238 genes) and differs
about twofold in gene count. Nothing here uses it; the discrepancy is noted
because any count in this inventory is a count in one annotation.

### Why orthogroup clustering is not used

OrthoFinder was run and abandoned. It placed **zero** taro genes in the
phytoene synthase orthogroup: all three candidates appeared in
`Orthogroups_UnassignedGenes.tsv`, with 6.6% of genes unassigned overall. That
run used nine species, not the final seventeen, so whether clustering would
also fail on the full set is untested, and the rescued summary tables record the
species set explicitly. The method therefore uses targeted reciprocal best hit
with a declared root group, which asks a narrower question and answers it.

## 6. Candidate recovery

DIAMOND `blastp`, e < 1e-5, forward from each Arabidopsis anchor into the taro
proteome and reverse from every taro hit into the Arabidopsis proteome. Two
levels of reciprocity are recorded:

- `rbh_anchor` — the reverse best hit is the anchor locus
- `rbh_family` — the reverse best hit is the anchor **or** a declared member of
  the same family

The grid counts `rbh_family`, because the target is the family. Without it, a
taro NCED whose closest Arabidopsis relative is NCED5 rather than the NCED3
anchor is silently dropped. On this data the stricter test would have given the
same answer everywhere: all six taro NCEDs and all four near-complete taro DXS
genes name the anchor itself.

Coverage is DIAMOND's `qcovhsp`, the percentage of query residues inside the
aligned block. An earlier version used alignment length divided by query
length, which counts gap columns and therefore exceeded 100% wherever the
alignment contained insertions: LCYB 102%, LCYE 103%, ORANGE 107%, fibrillin
106%. The correction moved one call, Ces13250 against CCD4, from 70% to 69% and
so across the default cutoff. The tree resolved it (section 7).

**Every** forward hit at e < 1e-5 is written out, whether or not it passes any
threshold, and a reciprocal best hit is never suppressed from the printed
table. This is not a stylistic choice. The two PSY fragments Ces12496 and
Ces12497 fail both defaults, so an earlier print filter hid them and the output
appeared to show PSY recovering one gene while the file contained all three.

### Sensitivity

Alignment length is swept across 100, 150, 200 aa and query coverage across 50,
70, 80 percent, giving nine combinations. None of these thresholds has a
principled justification; they are conventional, which is the reason for
sweeping them rather than defending them. A target whose count moves across the
grid is reported as unstable.

Two properties of the grid itself became visible and are stated rather than
corrected after the fact:

- **The 200 aa alignment floor exceeds the length of some pathway proteins.**
  MDS and NSX align over 166 and 175 aa against proteins of roughly 230 and
  220 aa, so a complete alignment cannot meet it. For a target set spanning 220
  to 700 residues, absolute alignment length is the wrong parameter and
  coverage is the right one.
- **The 80 percent coverage tier sits above where genuine ortholog pairs at
  this divergence land.** Six single-copy targets (MCT, MDS, IDI, Z-ISO, VDE,
  NSX) have exactly one taro ortholog each, at 78, 72, 80, 79, 79 and 79
  percent coverage. All six are stable at 50 and 70 percent and vanish at 80.

Strip those two artifacts and two targets are genuinely threshold-dependent:
DXS at three or four, and CCD4 at two or three.

### Split gene models

Query start and end are recorded for every hit, because coverage alone cannot
distinguish one gene broken across two records from two separate partial genes:
26% and 34% sum to 60% whether the pieces cover different parts of the anchor
or the same part twice. Pairs of reciprocal hits whose aligned query ranges
barely overlap while together covering much more of the anchor than either
alone are reported as split-model candidates on sequence evidence, to be
confirmed or refused on coordinates in section 8.

## 7. Family assignment by phylogeny

### Trees are built per family, not per target

CCD1, CCD4 and the NCEDs are all carotenoid cleavage dioxygenases and recover
from one pool of twelve taro genes, which reciprocal best hit partitions one /
three / six. **That partition is a result.** A tree containing only one
target's candidates has assumed it and cannot distinguish a CCD4 paralog from a
divergent NCED. The three therefore share one tree. PSY and GGPPS are in no
group and keep a tree each.

PDS/ZDS/CRTISO, LCYB/LCYE and CYP97A/CYP97C are also same-family sets. All are
detection targets, so no tree is built, but the grouping is recorded because
their candidate counts are not independent of one another.

### The root group, and the absence of an outparalog for the CCDs

| job | rooted on | relationship to the family |
|---|---|---|
| PSY | SQS1, SQS2 | squalene synthase: different substrate, product and pathway |
| GGPPS | FPS1, FPS2 | farnesyl diphosphate synthase: different product |
| CCD | CCD7, CCD8 | **inside the family**, the strigolactone branch |

Arabidopsis has nine genes in the carotenoid cleavage dioxygenase family —
CCD1, CCD4, CCD7, CCD8 and NCED2/3/5/6/9 — and every one is inside the family
being measured. **There is no Arabidopsis outparalog for this family.** That is
a fact about the family, not a flaw in the design, and the CCD tree is not
rooted the way the other two are. It is rooted on the strigolactone-branch
members, a sister clade of established identity, and the CCD1, CCD4 and NCED
groups are delimited within that rooted tree by their own Arabidopsis members.

In one respect the CCD boundary is better constrained than PSY's: the
CCD4-versus-NCED partition is anchored by six Arabidopsis genes, one CCD4 and
five NCEDs, where the PSY family boundary rests on a single split. An
Arabidopsis-anchored boundary can mislead — Arabidopsis lost PSY subgroup E3,
and a clade rule built on that assumption discarded a published rice PSY in an
earlier attempt — and six anchors across a partition is a sturdier constraint
than one.

Whether CCD7 and CCD8 are monophyletic was treated as an open question, since
CCD8 is often the most divergent lineage in the family and rooting on a
non-monophyletic pair would repeat an error made earlier in this project. They
are monophyletic (section 9).

### Two trees per job

**Tree A** is the family plus the root group and asks one question: is the root
group monophyletic? If yes, the family boundary can be placed. If no, it cannot
and nothing further is claimed from it.

**Tree B** is the family alone and exists for resolution, because a divergent
root group costs alignment columns the close relatives need.

Tree B **cannot be rooted**. In a gene family tree the genes of any one species
do not form a clade — each Amborella CCD groups with its own orthologs, not
with the other Amborella CCDs — and tree B is by construction the family with
its outgroup removed. It is therefore read by splits only, which require no
root.

Every recovered sequence is assigned to whichever query it hits with the highest
bitscore, so a sequence recovered by a root-group query belongs to the root set
and a sequence recovered by a family query belongs to the family set. An earlier
version searched both query sets together and built tree B by removing only the
*Arabidopsis* root sequences, leaving every other species' squalene synthase in
place: PSY tree A had 72 tips and tree B had 70, where a family-only tree should
have lost the whole eighteen-tip SQS clade.

### Which sequences enter a tree

Two floors, for two roles.

- **Reference tips** are scaffolding and must be near-complete in both
  directions: qcov ≥ 70% **and** scov ≥ 70%. Half a protein supplies half a
  column set and gaps for the rest. Seven GGPPS-family queries at 50% coverage
  against 454,114 sequences pulled in cis-prenyltransferases and polyprenyl
  synthases of very different lengths, producing a 3,441-column alignment for a
  family of roughly 350-residue proteins.
- **Taro candidates** are the question, not the scaffolding, and keep the
  permissive cell of the sensitivity grid: alignment ≥ 100 aa and qcov ≥ 50%.

A taro sequence below that floor is not discarded. It is reported as
unplaceable by phylogeny and decided on coordinates instead. The two PSY
fragments are the case in point: a 113-residue piece of a 437-residue protein
contributes mostly gaps to an alignment, but adjacent coordinates on one strand
are decisive.

### Alignment and inference

FAMSA, then trimAl `-gt 0.80`, then IQ-TREE with `-m MFP`, 1000 ultrafast
bootstrap replicates and 1000 SH-aLRT replicates.

Trimming uses a fixed gap threshold rather than `-automated1`. `-automated1`
selects between gappyout and strict from each alignment's own statistics, so
two alignments can be trimmed by different methods and their retention
percentages do not measure the same quantity. The `-automated1` column count is
still computed and reported beside the fixed one so the change is auditable
against earlier figures.

**Retention percentage is a poor diagnostic and is not used as one.** It
measures how gappy the untrimmed alignment was. What matters is the absolute
number of retained columns against the length of the proteins in the alignment:
504 columns for a family of ~600-residue proteins, 319 for ~350, 386 for ~430.

### The three tests

1. **Membership.** Is the root group one side of a single edge in tree A, and
   which side is each taro tip on? Tested as a **split**, not a rooted clade: a
   set is monophyletic in an unrooted tree if and only if some edge separates
   it, so a clade whose terminals equal the set or its complement counts. This
   sidesteps rooting entirely. For a job with one Arabidopsis reference this is
   the only question a tree can answer, and it answers it well.

2. **Assignment**, in tree A rooted on its root group, which is valid because
   test 1 establishes the root group is monophyletic. For each taro tip and
   each target, take the common ancestor of the tip and that target's
   references and keep it only if it excludes every other target's references.
   The smallest surviving one is the assignment. **If none survives the tip is
   recorded as not determinable**, with no fallback to a least-bad option.

3. **Confirmation**, in tree B, by the smallest side of any edge containing a
   target's references plus its assigned taro tips and excluding every other
   target's references. That clade is the orthogroup and the other species
   belong in it; an exact-split test would always fail, because Arabidopsis plus
   taro is a strict subset of a seventeen-species clade.

IQ-TREE writes an unrooted tree and `common_ancestor` is a rooted operation. An
earlier version of the reading step applied it to the newick as parsed, with
whatever root its parenthesisation implied, and the five Arabidopsis NCEDs
reported a common ancestor spanning the entire tree. That is the same class of
error as `Bio.Phylo.get_path()` never evaluating the root, which discarded a
published rice PSY in an earlier attempt, and as rooting on a single tip.

### Support

**SH-aLRT ≥ 80 AND UFBoot ≥ 95, both required.** A node meeting one and not the
other is reported as below threshold. FastTree SH-like values used in
exploratory work are an approximation, not a test, and nothing here uses them.

### Stopping rule

If the two trees disagree about a tip and neither supporting node clears
support, the assignment is recorded as not determinable and no third analysis is
run.

## 8. Genomic locus

Two family members occupy independent loci if they are on different anchored
sequences, or on the same anchored sequence with non-overlapping coordinates
and at least one annotated gene between them. A member on an `unanchor*`
sequence is reported but contributes no copy, because a haplotig of a gene
already counted cannot be distinguished from a paralog by coordinates alone.

Anchored is read from the assembly's own naming, not from a length threshold. A
cumulative-length rule tried first put the cut between Superscaffold13 at
113 Mb and Superscaffold12 at 105 Mb, treating a 105 Mb sequence carrying 1,963
genes as unplaceable.

### Split models need sequence evidence and genome evidence

**Genomic distance is not a criterion.** With 28,253 genes in 2,318 Mb the mean
spacing is 82 kb, so two genes tens of kilobases apart with nothing annotated
between them is the most ordinary arrangement in this annotation rather than a
signal. An earlier version used "no annotated gene between" alone and merged
Ces13250 and Ces13251, which are 75,073 bp apart, into one locus.

A split model requires all four of:

- aligned query ranges overlapping by at most 25% of the shorter — the pieces
  must cover **different** parts of the anchor
- query order matching genomic order given the strand — the pieces must be
  collinear
- the same strand
- no annotated gene between

The combined genomic span is reported against the family's own gene spans as
context rather than as a criterion.

The query-overlap condition does most of the work. Ces13250 and Ces13251 cover
anchor residues 181-592 and 70-592: the same part of the same protein twice, so
they cannot be two pieces of one gene whatever their spacing. The collinearity
condition refused Ces02307 and Ces02308, which are 44 bp apart on one strand
with nothing between and zero query overlap, but cover 585-686 then 64-206
going downstream.

### Over-merged models

The mirror error is one record spanning more than one gene. A candidate whose
protein exceeds twice the median length of its family's reference proteins —
taken from tree B's input, which holds the other sixteen species and no taro —
is flagged, its genomic span and exon count reported, and its copy number
withheld rather than reported as verified.

A flagged protein is then searched against Arabidopsis keeping all HSPs, and
the HSPs are **merged into query regions** before judgement. Several
Arabidopsis genes matching one region is a domain matched by a superfamily.
Different genes matching different regions is a merge. The quantity that
matters is how much of the protein matches nothing at all.

## 9. What the method establishes, and what it does not

Established:

- **the root group is monophyletic in all three jobs**, so the family boundary
  can be placed in each: CCD7/CCD8 at SH-aLRT 97.5 / UFBoot 100, FPS1/FPS2 at
  100/100, SQS1/SQS2 at 100/100. The concern that CCD7 and CCD8 might not form
  a clade was tested and is unfounded.
- **family membership** for every taro candidate, at the support of its job's
  split.
- **three CCD4 loci**, recovered as a 21-tip clade in both trees at UFBoot 99,
  including the gene whose reciprocal hit turned on one percentage point of
  coverage.
- **one full-length PSY locus** plus **one further PSY-like locus** carrying a
  split model, confirmed independently by non-overlapping collinear anchor
  coverage and by adjacent coordinates on one strand with nothing between.
- **presence** of both ORANGE chaperones and of fibrillin.

Not established, and recorded as such:

- **the NCED family size as a copy number.** Both trees group the same six
  taro genes with the Arabidopsis NCEDs and neither holds the grouping under
  resampling: UFBoot 59 in tree A and 43 in tree B against a declared threshold
  of 95. Reported as a detection result.
- **the GGPPS copy number.** The single candidate is a 2,232-residue model in a
  family whose reference median is 360, and 85% of it matches nothing in
  Arabidopsis. The prenyltransferase domain is a family member at 100/100; the
  model around it is unresolved.
- **PSY subgroup assignment** (M1 versus M2), withdrawn in earlier work as not
  determinable and not revisited here.
- **any copy number for a detection target.** Twenty-one targets carry lower
  bounds only.

## 10. Reproducing

```bash
source envs/activate.sh
python scripts/02_pathway_set.py          # anchors, KEGG cross-check
python scripts/02b_annotate_kegg.py       # name the KEGG-only loci
python scripts/02c_propose_scope.py       # propose a decision for each
python scripts/02d_finalise_scope.py      # resolve the remainder, freeze
python scripts/02e_group_families.py      # de-duplicate, record family groups
bash   scripts/03_preflight.sh            # schema and input check
bash   scripts/03_homology.sh             # candidates, grid, split candidates
bash   scripts/04_families.sh             # trees
python scripts/04b_read_trees.py          # read them
python scripts/05_report.py               # coordinates, master table
```

Outputs are in `results/tables/`. `part1_inventory.tsv` is the master table:
one row per target-and-gene, carrying the homology evidence, the phylogenetic
assignment with its support, the genomic locus, and one `evidence` column
stating what the row establishes.

Every correction made during the analysis is recorded in `LOGBOOK.md`, with the
superseded scripts in `scripts/archive/` and the superseded tables in
`results/tables/archive/`.
