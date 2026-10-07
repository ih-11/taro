# Open questions

## For the collaborator meeting

1. **P1PSY1 tissue specificity.** In tobacco T1, GUS appeared in leaf/stem
   vascular bundles and was *not* observed in roots. The taro target is corm
   storage parenchyma. A vascular-inducible eudicot promoter driving a monocot
   storage organ is a mismatch worth confronting before construct design.
   (Source: their own dissertation data — raise it as close reading.)

2. **Which step breaks?** mRNA → ribosome → protein stability → sink capacity
   are four separable steps. OR acts on PSY *protein stability* (holdase
   chaperone) — ribo-seq does **not** measure that. Be precise about which
   assay covers which step; overclaiming here gets caught.

3. **Baseline gap.** Taro β-carotene ~0.01–0.10 mg/100 g vs Carvita25 cassava
   ~2.2–2.3 mg/100 g — a 20–200× gap. But β-carotene had the highest GCV
   (52.78%) and GAM (98.61%) of any trait in taro germplasm surveys, and
   Champagne et al. achieved >4× by breeding alone. Has anyone screened the
   Indonesian landraces first?

4. **Does taro have one PSY or several?** Changes what "add a cassava PSY"
   actually means. Answerable in Part 1, before the visit.

5. **Authorship and framing** — agreed before analysis starts, not after.

## Internal unknowns

- Does `annotation_for_longread_mod_*.py` scale from 110 Mb (Chlamydomonas)
  to 2.3 Gb at 85% repeat? Test memory assumptions early.
- Chromosome naming differs across the three assemblies.
- minimap2 secondary-alignment defaults matter far more in a repeat-dominated
  genome than in anything this pipeline has been run on.
- `/tmp` is tmpfs — TMPDIR is pinned in `envs/activate.sh`. Don't bypass it.

## Findings from the fetch phase (2026-10-07)

**Two taro annotations differ twofold in gene count.**
`GCA_009445465.1` (Bellinger et al. 2020, the assembly Zou et al. mapped to)
reports 56,238 protein-coding genes. Sun et al. 2026 report 28,253 — the
proteome used here, confirmed by direct count after isoform collapse.

Same species, ~2x apart. Direct support for the project's premise that taro's
annotation is unstable, obtained before any analysis was run. Worth stating in
the introduction.

**No Araceae relative has a protein set at NCBI.**
Probed *Spirodela polyrhiza*, *Pistia stratiotes*, *Amorphophallus konjac*,
*Zantedeschia elliottiana*: assemblies only, no annotation. The four aroids
Sun et al. used for homology-based prediction are therefore unavailable to us
as proteomes.

Substituted the species set from Yin et al. 2021, who used it for taro
gene-family clustering — defensible, and anchored to a published taro
analysis. *Spirodela*, the closest relative at ~73 MYA, is absent.
**State as a limitation.**

**Taro isoform content.** 34,340 transcripts across 28,253 genes;
3,775 genes multi-isoform, max 12. From a short-read-derived annotation —
a lower bound on what long reads would find.
