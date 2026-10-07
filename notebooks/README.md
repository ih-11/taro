# notebooks

Exploration and figure generation. Pipeline steps belong in `scripts/`.

```
01_read_stats.ipynb        QC + the depth decision gate
02_isoform_concordance.ipynb   UpSet, Jaccard  -> F2
03_locus_tracks.ipynb      browser-style tracks -> F3
04_pathway_table.ipynb     orthology summary    -> F4
```

Register the kernel once:

```bash
conda activate taro
python -m ipykernel install --user --name taro --display-name "Python (taro)"
```
