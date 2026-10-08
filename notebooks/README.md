# Notebooks

## `part1_landscape.ipynb`

Four figures from `results/tables/part1_inventory.tsv` and the two tables beside
it. Nothing is computed here that is not already in those files, and **nothing
is written**: every figure is displayed, none is saved. The analysis lives in
`scripts/`, the method in `docs/METHODS_part1.md`.

```bash
source envs/activate.sh
jupyter lab notebooks/part1_landscape.ipynb
```

### What the figures encode, and why it is not gene counts

A bar chart of genes per pathway step would be the obvious figure and would
undo the method, which exists to separate three kinds of statement a count
treats alike: a copy number, a detection lower bound, and a presence call. So
the encoded variable is the `evidence` column, and claim kind is printed beside
each target rather than folded into the same colour.

| figure | question |
|---|---|
| 1 | per target, what evidence did each gene reach, and do gene models differ from countable loci |
| 2 | for the copy-number families, what is each node's support, and is the node testing membership or sub-assignment |
| 3 | where on the Arabidopsis anchor does each model sit — the split-model evidence, drawn |
| 4 | does the candidate set survive the nine threshold combinations |

Figure 3 carries the most weight. Two bars tiling the PSY anchor end to end
beside three bars lying across the same CCD4 residues is the whole
split-annotation-versus-tandem-duplication argument, and no amount of prose
about genomic distance does the same work.

### Colour

One sequential blue hue for the ordered evidence ladder, two categorical hues
for the two support statistics. Both were validated rather than eyeballed, with
the dataviz skill's checker:

```
"#86b6ef,#5598e7,#2a78d6,#1c5cab,#0d366b"  --ordinal     all checks pass
"#2a78d6,#eb6834"                 --pairs all   CVD ΔE 24.7, normal 33.6
```

**Five tiers is not a preference, it is what passed.** A six-step ramp fails the
adjacent-lightness check anywhere in this hue's usable band, so eleven distinct
`evidence` strings collapse to five rungs. The full mapping prints in the
notebook before any figure: a figure that hides which categories were merged is
one you cannot check.

The ladder is ordered **within** a claim kind, not across. `present` for a
presence target is the complete claim that target ever made and sits on the
same rung as `candidate (RBH)`, because both rest on reciprocal homology with
adequate coverage.
