# Logbook

A running record of what was done, what it produced, and why decisions were
made. Commit messages in this repository are deliberately terse. This file
carries the reasoning instead, so that someone picking the project up later
(including a future version of me) can reconstruct not just the steps but the
thinking behind them.

Entries are newest last.

---

## 2026-10-07  Project setup

### Environment

Created a dedicated conda environment rather than extending the existing
general purpose one. The existing environment carries TensorFlow, PyTorch and
the single cell analysis stack, and adding bioconda genomics tools to it
invites a dependency solver conflict that costs an afternoon to untangle.

Verified the machine before starting: 28 threads, 88 GiB RAM, 23 GiB swap.
The heaviest single step in the whole project is building a minimap2 splice
index for a 2.3 Gb genome, which needs roughly 15 to 20 GiB. There is ample
headroom.

Two storage observations shaped the layout:

1. `/tmp` is tmpfs, meaning it is backed by RAM, not disk. Writing a large
   sort there consumes memory. `TMPDIR` is therefore pinned to a directory
   inside the repository.
2. `/mnt/f` is a Windows drive reached through the WSL v9fs layer, which is
   noticeably slower than native ext4 for the many small reads and writes that
   alignment produces. The home directory is ext4 with 869 GB free.

The final arrangement puts code in `~/Code/taro` (ext4, version controlled)
and large files in `/mnt/f/RA/Downstream/Project5_TARO` (the user's existing
project convention), connected by a gitignored symlink.

### Taro proteome

Downloaded from Figshare (doi 10.6084/m9.figshare.29917034), the data deposit
accompanying Sun et al. 2026, Scientific Data 13:802.

An initial attempt with `aria2c` failed with HTTP 403. Figshare's download
endpoint redirects to S3 and rejects that client's user agent. `curl -L`
handles the redirect correctly and was used instead.

**Result:** 34,340 protein sequences.

**Problem:** the paper reports 28,253 protein coding genes. The file therefore
contains alternative transcripts, at a ratio of about 1.22 transcripts per
gene.

**Why this matters.** OrthoFinder assumes one sequence per gene. If isoforms
are left in, each gene can scatter across its own orthogroup, which inflates
apparent gene family size. Since the entire purpose of this analysis is to
count how many phytoene synthase paralogs taro has, feeding in isoforms would
corrupt the one number we care about, and would do so in the direction that
produces a falsely interesting result.

**Fix.** Headers follow the pattern `>rna-Ces#####.N`, where `Ces#####` is the
gene and `.N` is the isoform index. Collapsed to the longest isoform per gene
and renamed to the gene identifier.

**Verification:** 28,253 sequences out, matching the published gene count
exactly. 3,775 genes carry more than one isoform; the maximum is 12.

That isoform figure is worth remembering. It comes from a short read derived
annotation, so it is a lower bound on the isoform diversity that long reads
would be expected to reveal.

### Outgroup species: an unplanned substitution

The original plan was to use the four Araceae species that Sun et al.
themselves used for homology based gene prediction: *Amorphophallus konjac*,
*Zantedeschia elliottiana*, *Pistia stratiotes*, and *Spirodela polyrhiza*.
Using the same species would have anchored our orthology to the assumptions
of the published annotation.

Probing NCBI showed this is not possible:

| Species | At NCBI | Protein set |
|---|---|---|
| *Spirodela polyrhiza* | 9 assemblies | none |
| *Amorphophallus konjac* | 1 assembly | none |
| *Pistia stratiotes* | not found under this name | none |
| *Zantedeschia elliottiana* | not found under this name | none |

All four are assembly only. Retrieving proteomes would mean chasing Chinese
repositories (NGDC/CNCB) and journal supplements, which would consume most of
the available time for an input that is helpful but not essential.

**Substitution.** Yin et al. 2021 (Molecular Ecology Resources 21:68) performed
gene family clustering for taro against a ten species set and published every
accession. Eight of those nine non taro species have annotated proteomes at
NCBI. Reusing their set is defensible precisely because it is the set a
published taro paper used for this exact kind of analysis.

Species used: *Arabidopsis thaliana*, *Oryza sativa*, *Zea mays*,
*Musa acuminata*, *Solanum tuberosum*, *Nelumbo nucifera*,
*Amborella trichopoda*, *Zostera marina*.

This gives three monocots, a tuber forming eudicot with an established
carotenoid literature (potato), a basal eudicot, a basal angiosperm outgroup,
and Arabidopsis as the functional reference for PSY and OR.

**Limitation to state in the methods.** *Spirodela polyrhiza* is taro's closest
sequenced relative, having diverged roughly 73 million years ago, and it is
absent from our set. Orthogroup boundaries near taro are therefore less
sharply resolved than they could be. This does not affect the paralog count
itself, which depends on taro's own sequences, but it does limit what can be
said about when any duplication occurred.

### Finding: two taro annotations disagree twofold

While probing NCBI, the 2019 taro assembly `GCA_009445465.1` (Bellinger et al.
2020, G3 10:2763) was found to report **56,238 protein coding genes**.
Sun et al. 2026 report **28,253**, which we confirmed by direct count.

Same species. Roughly a factor of two apart.

This is worth highlighting because it is exactly the instability the project
set out to investigate, and it surfaced from a metadata query before any
analysis had been run. It belongs in the introduction of whatever is written.

It is also directly relevant to Part 2, because `GCA_009445465.1` is the
assembly that Zou et al. 2025 mapped their ONT corm reads to. Re-mapping those
same reads to the much better 2026 assembly is the core experiment of Part 2.

### Outgroup proteomes: a second isoform problem

First attempt (`05_fetch_outgroups.sh`) grouped isoforms by their description
text, because RefSeq protein FASTA headers do not carry a gene identifier.

This failed, and failed quietly. Arabidopsis came out at 14,297 genes against
an expected 27,562, roughly half. The cause is that thousands of genuinely
distinct genes share a description string: "F-box protein", "disease
resistance protein", "uncharacterized protein". Grouping on description merged
them.

**Why this was dangerous.** The error collapses paralogs, which is the precise
opposite of the isoform problem but corrupts the same number. It was caught
only because the Arabidopsis gene count is a well known figure and the output
looked obviously wrong.

**Fix** (`06_rebuild_primary.sh`). Download the GFF3 alongside the protein
FASTA and build the protein to gene mapping from the CDS feature attributes,
where `protein_id=` and `Dbxref=GeneID:` appear together. Group on the real
gene identifier, keep the longest isoform.

**Verification:** every species now matches its NCBI reported gene count.

| Species | Proteins in | Genes out | Expected |
|---|---|---|---|
| *Arabidopsis thaliana* | 48,265 | 27,562 | 27,562 |
| *Oryza sativa* | 42,580 | 28,740 | 28,738 |
| *Zea mays* | 51,764 | 33,461 | 33,461 |
| *Musa acuminata* | 47,707 | 30,737 | 30,737 |
| *Solanum tuberosum* | 37,967 | 28,404 | 28,404 |
| *Nelumbo nucifera* | 38,191 | 24,073 | 24,073 |
| *Amborella trichopoda* | 31,494 | 17,106 | 17,106 |
| *Zostera marina* | 20,648 | 20,436 | 20,436 |

Rice is two over because two proteins have CDS entries with no GeneID in the
GFF3. They are retained as singletons, which is harmless.

### A note on error checking

Two errors were made and caught during this phase, both in isoform handling,
and both would have silently distorted the paralog count. The pattern is worth
naming: when a preprocessing step changes the number of sequences, check the
new number against an externally known value before proceeding.

This motivates a sanity check to apply once OrthoFinder finishes. Arabidopsis
has exactly one phytoene synthase gene, AT5G17230, and this is not in dispute
anywhere in the literature. If the PSY orthogroup contains two Arabidopsis
genes, something upstream is wrong and the taro count cannot be trusted
either. The check costs nothing and catches the class of error that has
already occurred twice.

### Next

OrthoFinder across all nine proteomes, then extraction of the carotenoid and
MEP pathway orthogroups, then the pathway table.

## 2026-10-07 evening  Part 1 results

### OrthoFinder completed

28 minutes on 24 threads. 223,125 genes (93.4%) assigned to 18,726
orthogroups; 8,597 orthogroups contain all nine species, 1,412 of those
entirely single copy. STRIDE rooted the species tree on *Amborella* with 2,155
duplications supporting against 8 contradicting, which is the expected
topology and a good sign the species set behaves sensibly.

### OrthoFinder reported zero taro phytoene synthase genes

This cannot be true. Taro accumulates carotenoids, so it has a PSY.

Direct homology search (`10_`) found three PSY-like taro proteins, and all
three appear in `Orthogroups_UnassignedGenes.tsv`. MCL failed to cluster them
at all rather than placing them wrongly. *Zostera marina* was also absent from
that orthogroup, which fits the same explanation.

6.6% of all genes went unassigned in this run. Any pathway gene could have met
the same fate, so **orthogroup membership systematically undercounts**. This
is a methodological finding worth stating: for targeted gene family questions,
clustering output should be checked against direct homology search rather than
trusted on its own.

### Pathway inventory by reciprocal best hit

Reciprocal best hit was used instead. Note that this does detect duplications,
contrary to a concern raised earlier: multiple taro paralogs can share the same
best Arabidopsis match, and DXS, HDR, CCD4 and NCED3 all returned more than one.

**Phytoene synthase: one full-length ortholog.** Ces24605, 78.0% identity,
83% query coverage, e = 1.4e-202, reciprocal best hit to AT5G17230 with no
other Arabidopsis hit above threshold.

Ces12497 and Ces12496 cover only 34% and 26% of the query. Adjacent gene IDs,
and their coverage sums to roughly one protein. Most likely a single gene
split across two models by the annotation. To be resolved by gene tree before
the copy number is stated publicly.

**Carotenoid cleavage dioxygenase 4: three copies.** Ces13251, Ces13250,
Ces03723, all reciprocal best hits at 86 to 90% coverage.

This matters more than the PSY result for the collaboration. CCD4 degrades
carotenoids, and in potato and peach CCD4 variation is a principal determinant
of flesh colour. Taro carrying one phytoene synthase against three CCD4
enzymes suggests that increasing flux through PSY alone may not produce the
expected accumulation, because the degradation step is amplified.

**Caveat.** Ces11733 and Ces27429 appear under both CCD4 and NCED3. CCD and
NCED belong to the same superfamily and cross-hit readily. The three-copy CCD4
claim requires a gene tree before it is stated to anyone.

**ORANGE: present as a clean pair.** Ces06664 (OR) and Ces26954 (OR-like),
each the reciprocal best hit of its own Arabidopsis counterpart. Taro has the
chaperone machinery that controls PSY protein stability.

**Also expanded:** DXS four copies, HDR two. DXS is the entry point to the MEP
pathway and its expansion is common in plants.

### Next

Gene trees, not more best hits. Two questions need topology to answer:
whether Ces12496/12497 fall inside the PSY clade or outside it, and whether
the three CCD4 candidates are genuinely CCD4 rather than NCED.

## 2026-10-07 late evening  Gene trees

### Why trees were needed

Reciprocal best hit gave an inventory but left two things unresolved, both of
which depend on topology rather than pairwise similarity.

First, taro showed one full-length phytoene synthase (Ces24605, 428 aa) plus
two short matches, Ces12497 and Ces12496, covering only 34% and 26% of the
Arabidopsis query. Second, three taro genes were called CCD4 orthologs, but
two of them also matched NCED3. CCD and NCED belong to the same superfamily
and cross-hit readily, so that result could not be trusted as it stood.

### The PSY fragments are one broken gene model, not a paralog

Genome coordinates settle it without needing a tree:

| Gene | Scaffold | Start | End | Strand | Protein |
|---|---|---|---|---|---|
| Ces24605 | Superscaffold13 | 96,419,008 | 96,426,240 | − | 428 aa |
| Ces12496 | Superscaffold7 | 79,028,474 | 79,030,838 | + | 112 aa |
| Ces12497 | Superscaffold7 | 79,031,097 | 79,032,356 | + | 155 aa |

Ces12496 and Ces12497 sit on the same scaffold and strand, separated by 259 bp,
with nothing else annotated between them. Their alignments to Arabidopsis PSY
cover non-overlapping regions. This is a single locus split across two gene
models by the annotation.

Their combined length is 267 aa against a 437 aa reference, so they recover
only 61% of a full phytoene synthase even when summed. Either the model is
missing a further segment, or the locus is degenerate. Distinguishing those two
possibilities requires the long-read data, which is Part 2.

Both sequences also fell below the 150 aa threshold used for tree building, so
they are too fragmentary to place phylogenetically. That is itself the finding.

This makes a clean illustration of the project's premise. The published
annotation gives one intact copy of the pathway's rate-limiting enzyme and one
locus that cannot be interpreted without better data.

### Cassava added to the proteome set

Cassava was not in the Yin et al. species set, because that set was assembled
for taro gene-family clustering rather than carotenoid engineering. But cassava
is the source of the transgene this collaboration proposes to move, and it
carries the PSY paralogs that motivate the entire copy-number question.
Without it the PSY tree had no anchor to the gene being transferred.

Added GCF_001659605.2: 49,290 proteins collapsing to 29,735 genes, matching the
NCBI count exactly.

The OrthoFinder run remains a nine-species result while the gene trees are
ten-species. This is deliberate. Orthogroups were a screening step and the
trees are the evidence; rerunning OrthoFinder for one additional species would
cost another 28 minutes and change none of the conclusions.

### CCD4 and NCED separate cleanly

Patristic distance from each taro gene to the two Arabidopsis anchors,
AT4G19170 (CCD4) and AT3G14440 (NCED3):

| Taro gene | to CCD4 | to NCED3 | Assignment |
|---|---|---|---|
| Ces13251 | 0.966 | 2.133 | CCD4 |
| Ces13250 | 1.194 | 2.362 | CCD4 |
| Ces03723 | 1.309 | 2.477 | CCD4 |
| Ces11733 | 2.056 | 0.808 | NCED |
| Ces27429 | 2.059 | 0.810 | NCED |
| Ces09230 | 2.071 | 0.823 | NCED |
| Ces28200 | 2.141 | 0.893 | NCED |
| Ces16312 | 2.142 | 0.894 | NCED |

Every gene is roughly twice as far from the wrong anchor as from the correct
one. There is no ambiguity.

**The earlier caveat is resolved.** Taro has three CCD4 genes: Ces13251,
Ces13250 and Ces03723. The cross-hits that prompted the caveat were genuine
NCED genes, which is expected within a superfamily, and not evidence against
the CCD4 count.

One point remains open: Ces13250 and Ces13251 have adjacent gene identifiers,
so they may be a tandem duplication or, as with the PSY fragments, another
split model. Worth checking coordinates.

### Taro has one PSY where most plants have three

Sequences recovered per species in the PSY family search:

| Species | PSY-family proteins |
|---|---|
| *Colocasia esculenta* (taro) | 1 |
| *Zostera marina* | 1 |
| *Nelumbo nucifera* | 2 |
| *Manihot esculenta* (cassava) | 3 |
| *Arabidopsis thaliana* | 1 |
| *Oryza sativa* | 3 |
| *Zea mays* | 3 |
| *Musa acuminata* | 3 |
| *Amborella trichopoda* | 3 |
| *Solanum tuberosum* | 4 |

**These counts are from a homology search with a 150 aa filter, not from tree
topology, and have not yet been verified against the tree.** Some may belong
to the broader squalene and phytoene synthase superfamily rather than being
true phytoene synthases. The count must be confirmed before being quoted.

Taro's nearest neighbours by patristic distance from Ces24605:

| Distance | Protein | Species |
|---|---|---|
| 0.396 | XP_009386425.1 | *Musa acuminata* |
| 0.432 | XP_010270744.2 | *Nelumbo nucifera* |
| 0.440 | XP_021592587.1 | *Manihot esculenta* |
| 0.443 | XP_015164579.1 | *Solanum tuberosum* |

The newick shows taro grouping with *Musa* inside a clade that also contains
*Zostera*, which is the expected monocot topology and a sign the tree is
behaving sensibly.

Note that *Zostera marina*, the only other member of the Alismatales in this
set, also has a single copy. The low count may therefore be a characteristic
of the Alismatales rather than something specific to taro. This is precisely
the question *Spirodela polyrhiza* would settle, which turns the earlier
limitation from a generic caveat into a concrete gap.

### What this means for the collaboration

Taro carries a single phytoene synthase, while cassava, the source of the
proposed transgene, carries more. Taro also carries three CCD4 carotenoid
cleavage dioxygenases, the enzymes that degrade the compound the project aims
to accumulate. In potato and peach, CCD4 variation is a principal determinant
of flesh colour.

Raising flux through a single PSY into a tissue carrying a triplicated
degradation step may not produce the accumulation the plan assumes. This is a
concrete, testable concern derived entirely from public data, and it bears
directly on construct design.

### Next

1. Verify PSY copy number per species from tree topology rather than hit
   counts. This is the claim most likely to be wrong and most likely to be
   repeated.
2. Figure: two panels, the PSY tree with taro and cassava highlighted and the
   CCD/NCED tree with clades coloured.
3. Check whether Ces13250 and Ces13251 are a tandem duplication or a split
   model, using the same coordinate approach that resolved the PSY fragments.

## 2026-10-07 night  External critique and corrections

An independent read of the repository raised several points. Most were
correct and are recorded here with the resulting changes.

### Accepted: the CCD4 interpretation overreached

The README stated that taro's degradation step is "amplified" relative to
synthesis. Copy number is not flux. Three CCD4 paralogs could include genes
that are silent in corm tissue, expressed in other organs, or catalytically
dead. Nothing in this analysis measures expression or activity.

Corrected claim: taro's genome encodes three CCD4 paralogs and one intact
phytoene synthase. Whether that asymmetry affects carotenoid accumulation in
corm is a hypothesis worth testing, not a result. It remains worth raising
with the collaborator, because it identifies a variable their design does not
currently account for.

### Accepted: Ces13250 and Ces13251 were advertised before being checked

The logbook flagged these adjacent gene IDs as possibly a tandem duplication
or another split model, and the README then reported three independent copies
anyway. Coordinate check pending, same method that resolved the PSY fragments.

### Accepted: Part 3 does not need only the genome

The README claimed Parts 1 and 3 survive without long-read data, while Part 3
begins with empirical TSS. A transcription start site cannot be derived from
genome sequence. This was a plain inconsistency.

Further, ONT cDNA and PacBio Iso-Seq both suffer 5' degradation, so neither
establishes true TSS in the way CAGE or 5'-RACE does. What they provide is
observed 5' ends, which bound the start region without pinpointing it.

This is still an improvement on the dissertation's approach, which analysed a
guessed upstream window, but it must be described accurately. Part 3 is
comparative cis-regulatory analysis using long-read-supported 5' ends, not
promoter characterisation.

### Accepted: motif presence is not regulatory compatibility

Finding cassava promoter motifs in taro upstream regions is sequence evidence
only. It does not show that the promoter would drive expression in taro, and
certainly not in corm storage parenchyma. The tobacco data in the dissertation
showed vascular expression, not storage tissue. Part 3 is framed accordingly.

### Partially accepted: ONT and PacBio should be analysed separately

Agreed, but the stated reason understates the problem. The two datasets come
from different cultivars: Zou et al. used 'Lipu Taro No.1' and Sun et al. used
'Bun Long'. Differences between them therefore confound sequencing technology
with genotype, and neither factor can be isolated.

They should be treated as two independent attempts at the same question.
Agreement between them is the evidence; a pooled analysis would obscure
exactly what makes them informative.

### Accepted: Part 2 must hold everything else constant

Same reads, same minimap2 preset, same isoform caller, same parameters, with
only the reference assembly changing. This will be written into the script as
an explicit constraint rather than left as a convention, because it is the
single assumption the entire experiment rests on.

### Standing: the missing Araceae relative

Already recorded. The critique confirms it as the largest evolutionary
limitation, and the *Zostera* observation sharpens rather than resolves it.
