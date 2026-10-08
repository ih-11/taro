#!/usr/bin/env bash
# cleanup_store.sh — remove what the current pipeline does not use.
#
# The store holds 54,365 files. 54,058 of them are OrthoFinder per-orthogroup
# output: one sequence file and one alignment per orthogroup, across 18,726
# orthogroups. OrthoFinder is dropped from the method — it placed zero taro
# genes in the phytoene synthase orthogroup and contributed no evidence to any
# conclusion — so none of that output is read by anything.
#
# The rule applied here: delete what is (a) unused by the current pipeline,
# (b) regenerable from a command, and (c) whose findings are already captured
# in text. Anything failing one of the three stays.
#
# Specifically NOT deleted:
#   - orthofinder_input/*.fa      the 17 proteomes the pipeline reads
#   - the taro GFF3 and .pep      05_report.py reads the GFF3
#   - OrthoFinder summary files   Statistics, Orthogroups.tsv and especially
#                                 Orthogroups_UnassignedGenes.tsv, which is the
#                                 evidence that MCL failed to cluster the PSY
#                                 genes and is cited in LOGBOOK.md
#   - aratha_protein2gene.tsv     every later step joins on it
#   - the exploratory tree files  small, and the logbook cites support values
#                                 from them
#
# Run with --dry-run first. It reports sizes and deletes nothing.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

DRY=0
[ "${1:-}" = "--dry-run" ] && DRY=1

S="$TARO_WORK"
P="$TARO_REFS/proteomes"
MAN="$TARO_WORK/CLEANUP.md"

tot_before=$(du -sb "$(dirname "$S")" 2>/dev/null | cut -f1)
files_before=$(find "$(dirname "$S")" -type f 2>/dev/null | wc -l)

say () { printf "\n%s\n" "$1"; printf '%.0s-' {1..70}; printf "\n"; }
size () { du -sh "$1" 2>/dev/null | cut -f1; }

drop () {   # path, reason
    local p="$1" why="$2"
    [ -e "$p" ] || return 0
    local sz; sz=$(size "$p")
    local n;  n=$(find "$p" -type f 2>/dev/null | wc -l)
    printf "  %-9s %7s files  %s\n" "$sz" "$n" "${p#$TARO_STORE/}"
    printf "            %s\n" "$why"
    if [ "$DRY" -eq 0 ]; then
        rm -rf "$p"
        echo "| \`${p#$TARO_STORE/}\` | $sz | $why |" >> "$MAN.tmp"
    fi
}

echo "store: $TARO_STORE"
echo "before: $files_before files, $(numfmt --to=iec "$tot_before" 2>/dev/null || echo "$tot_before bytes")"
[ "$DRY" -eq 1 ] && echo "MODE: dry run, nothing will be deleted"

[ "$DRY" -eq 0 ] && : > "$MAN.tmp"

# ---------------------------------------------------------------- OrthoFinder
say "OrthoFinder per-orthogroup output"
for R in "$S"/05_orthology/run_*/Results_*; do
    [ -d "$R" ] || continue
    for d in Orthogroup_Sequences MultipleSequenceAlignments Gene_Trees \
             Resolved_Gene_Trees Single_Copy_Orthologue_Sequences \
             WorkingDirectory Phylogenetic_Hierarchical_Orthogroups \
             Putative_Xenologs Orthologues Phylogenetically_Misplaced_Genes; do
        drop "$R/$d" "OrthoFinder is not used by the current method"
    done
done

# ---------------------------------------------------------------- GFF3
say "outgroup GFF3 files"
echo "  The current pipeline reads only the taro GFF3. These were downloaded"
echo "  to map rice and sorghum PSY genes to chromosomes for the subgroup"
echo "  question, which is recorded as not determinable. Regenerable with"
echo "  'datasets download genome accession <ACC> --include gff3'."
echo
for g in anacom aspoff bradis elagui orysat phodac setita sorbic aratha; do
    drop "$S/05_orthology/${g}_gff" "outgroup GFF3, not read by the current pipeline"
done

# ---------------------------------------------------------------- duplicates
say "duplicate proteome copies"
echo "  Each proteome exists as *.primary.faa and again as"
echo "  orthofinder_input/*.fa. The pipeline reads orthofinder_input."
echo
for f in "$P"/*.primary.faa "$P"/colesc.tagged.faa; do
    [ -e "$f" ] || continue
    base=$(basename "$f" .primary.faa)
    base=${base%.tagged}
    if [ -s "$P/orthofinder_input/$base.fa" ]; then
        drop "$f" "duplicate of orthofinder_input/$base.fa"
    else
        echo "  KEEP  $(basename "$f") — no copy in orthofinder_input"
    fi
done

# ---------------------------------------------------------------- orphans
say "orphaned search databases and intermediates"
echo "  The exploratory trees/ directory built its own combined database."
echo "  04_families.sh builds its own under 04_families/ and rebuilds if"
echo "  missing, so these are orphaned."
echo
drop "$S/05_orthology/trees/all.faa"  "orphaned; 04_families.sh rebuilds its own"
drop "$S/05_orthology/trees/all.dmnd" "orphaned; 04_families.sh rebuilds its own"
drop "$S/05_orthology/pathway_homology" "superseded by work/03_homology"
drop "$S/05_orthology/psy_check" "exploratory; findings are in LOGBOOK.md"
for z in "$P"/*.zip; do drop "$z" "download archive, already extracted"; done
for x in "$P"/*_x; do drop "$x" "extraction directory, already processed"; done
drop "$P/.figshare_files.json" "API response cache"

# ---------------------------------------------------------------- report
say "result"
tot_after=$(du -sb "$(dirname "$S")" 2>/dev/null | cut -f1)
files_after=$(find "$(dirname "$S")" -type f 2>/dev/null | wc -l)
echo "  files : $files_before -> $files_after"
echo "  size  : $(numfmt --to=iec "$tot_before" 2>/dev/null) -> $(numfmt --to=iec "$tot_after" 2>/dev/null)"
echo "  freed : $(numfmt --to=iec $((tot_before - tot_after)) 2>/dev/null)"

if [ "$DRY" -eq 1 ]; then
    echo
    echo "  Nothing was deleted. Re-run without --dry-run to apply."
    exit 0
fi

{
    echo "# Store cleanup"
    echo
    echo "Run $(date -Iseconds)."
    echo
    echo "Files: $files_before -> $files_after."
    echo "Size: $(numfmt --to=iec "$tot_before") -> $(numfmt --to=iec "$tot_after")."
    echo
    echo "Everything removed is either unused by the current pipeline or"
    echo "regenerable. Findings that depended on any of it are recorded in"
    echo "\`~/Code/taro/LOGBOOK.md\` and in \`results/tables/archive/\`."
    echo
    echo "| Path | Size | Why |"
    echo "|---|---|---|"
    cat "$MAN.tmp"
    echo
    echo "## Regenerating"
    echo
    echo '```bash'
    echo "# outgroup GFF3, if ever needed again"
    echo "datasets download genome accession <ACC> --include gff3"
    echo
    echo "# OrthoFinder, though the method no longer uses it"
    echo "bash ~/Code/taro/scripts/archive/07_orthofinder.sh"
    echo
    echo "# proteome copies, if an archived script needs them"
    echo "cp refs/proteomes/orthofinder_input/<tag>.fa refs/proteomes/<tag>.primary.faa"
    echo '```'
} > "$MAN"
rm -f "$MAN.tmp"
echo
echo "  manifest: $MAN"
