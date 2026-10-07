# scripts

Numbered to match `work/` subdirectories. Each writes to its own `work/NN_*`
and logs to `logs/`.

```
00_fetch_refs.sh        assemblies + annotations + proteomes
01_fetch_reads.sh       SRA -> data/longread
02_qc.sh                read stats, the Part 2 decision gate
03_align.sh             minimap2 x 2 read sets x 3 assemblies
04_isoform.sh           IsoQuant / StringTie per combination
05_compare.sh           R1 sequence + R2 projection
06_orthology.sh         OrthoFinder vs Araceae
07_promoter.sh          TSS extraction + FIMO
```

Every script begins with:

```bash
source "$(dirname "$0")/../envs/activate.sh"
```
