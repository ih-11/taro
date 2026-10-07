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
