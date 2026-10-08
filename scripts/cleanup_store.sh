#!/usr/bin/env bash
# cleanup_store.sh — permanent removal of what the current method does not use.
#
# REPLACES the earlier cleanup_store.sh, whose header said 54,365 files. The
# store holds 114,114. The old figure counted only Orthogroup_Sequences and
# MultipleSequenceAlignments, 52,646 between them, and never walked
# WorkingDirectory, where OrthoFinder keeps a second renamed copy of the same
# sequences and alignments plus one tree per orthogroup: another 59,776 files.
# Enumerating directory by directory is how a total stays invisible, so this
# version measures the whole store before proposing anything.
#
# WHERE THE FILES ARE
#
#   113,955   work/05_orthology/run_20261007_2049/   one abandoned run
#       159   everything else, every input and every result included
#
# WHY THE RUN GOES
#
#   OrthoFinder is not part of the method. It placed zero taro genes in the
#   phytoene synthase orthogroup: all three candidates landed in
#   Orthogroups_UnassignedGenes.tsv, with 6.6% of genes unassigned overall.
#   That failure is why the method uses targeted reciprocal best hit, and no
#   per-orthogroup file is read by anything downstream.
#
#   THE RUN USED NINE SPECIES, NOT SEVENTEEN. Its Orthologues directory holds
#   ambtri aratha colesc musacu nelnuc orysat soltub zeamay zosmar. The eight
#   added later, anacom aspoff bradis elagui manesc phodac setita sorbic, are
#   absent. Two consequences, and both are reasons to be careful rather than
#   reasons to keep 113,955 files:
#
#     1. The run is NOT reproducible from the current inputs. Re-running
#        07_orthofinder.sh today clusters seventeen proteomes and produces a
#        different result. So this output cannot be regenerated, only redone.
#
#     2. The recorded failure was observed on nine species. Whether clustering
#        would also fail on the final seventeen is untested. The rescued files
#        must say so, or a later reader takes them for the 17-species analysis.
#        Phase 1 writes a README that states it.
#
#   The evidence is a handful of summary tables. The 113,955 per-orthogroup
#   files are not evidence of anything; they are the intermediate products of a
#   clustering whose summary already records what happened.
#
# THREE CORRECTIONS TO THE PREVIOUS VERSION OF THIS SCRIPT
#
#   a. It deleted *.primary.faa whenever orthofinder_input/<tag>.fa merely
#      EXISTED. Existence is not identity. This version compares the two by a
#      header-independent checksum of their sorted sequences and refuses to
#      delete anything that does not match.
#
#   b. It deleted psy_check/ and pathway_homology/ outright. Those hold the
#      search tables the PSY and pathway findings in LOGBOOK.md are quoted
#      from, and the rewritten 03_homology.py will not reproduce the same
#      numbers because the coverage metric changed from alignment-length over
#      query-length to DIAMOND qcovhsp. Deleting them would leave quoted
#      figures with no surviving source. The four .tsv files are rescued; only
#      the regenerable .dmnd and .faa are removed.
#
#   c. It deleted the *_gff directories including aratha_gff, which now holds
#      nothing but md5sum.txt. That checksum is the only record of which
#      Arabidopsis release produced aratha_protein2gene.tsv, and that mapping
#      file is itself not regenerable without the GFF3. Any checksum or
#      provenance file in a *_gff directory is rescued first.
#
# THE RULE
#
#   Delete only what is (a) unread by the current pipeline, (b) regenerable or
#   genuinely disposable, and (c) whose findings already exist in text or in a
#   rescued table. Anything failing one of the three stays.
#
# NEVER TOUCHED
#
#   refs/proteomes/orthofinder_input/*.fa      the 17 proteomes every step reads
#   refs/proteomes/Colocasia_*.gff3 and .pep   05_report.py reads the GFF3
#   refs/proteomes/manifest.tsv                provenance for all 17
#   work/05_orthology/aratha_protein2gene.tsv  every step joins on it, and it
#                                              cannot be rebuilt without the
#                                              Arabidopsis GFF3, which is gone
#   work/05_orthology/trees/*                  77 small files; LOGBOOK quotes
#                                              support values out of them
#   work/03_homology/{aratha,colesc}.dmnd      the new 03_ reuses these
#
# RUN ORDER
#
#   bash scripts/cleanup_store.sh            measure and plan, delete nothing
#   bash scripts/cleanup_store.sh --apply    rescue, verify, then delete
#
# EXPECT --apply TO BE SLOW. Removing 113,955 files from NTFS through the WSL
# filesystem bridge is a per-file operation. Ten to thirty minutes is normal.
# The script prints elapsed time. Do not interrupt it.
set -euo pipefail
: "${TARO_STORE:?source envs/activate.sh first}"
: "${TARO_WORK:?}"
: "${TARO_REFS:?}"

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

S="$TARO_STORE"
W="$TARO_WORK"
P="$TARO_REFS/proteomes"
PROV="$S/provenance"
MAN="$S/CLEANUP.md"
TMP="$S/.cleanup_rows"

rule () { printf '%.0s-' {1..76}; printf "\n"; }
say  () { printf "\n%s\n" "$1"; rule; }
hsz  () { du -sh "$1" 2>/dev/null | cut -f1 || echo "?"; }
bsz  () { du -sb "$1" 2>/dev/null | cut -f1 || echo 0; }
nfil () { local n; n=$(find "$1" -type f 2>/dev/null | wc -l); [ -f "$1" ] && n=1; echo "$n"; }

# header-independent checksum of the sequence set, so two FASTA files with
# different headers but the same sequences compare equal
seqsum () {
    awk '/^>/{if(s!="")print s; s=""; next}{s=s $0}END{if(s!="")print s}' "$1" \
        | sort | md5sum | cut -d' ' -f1
}

PLAN=()
plan () { [ -e "$1" ] && PLAN+=("$1"$'\t'"$2") || true; }

RESCUE=()      # src<TAB>destination subdirectory
resc () { [ -e "$1" ] && RESCUE+=("$1"$'\t'"$2") || true; }

echo "=========================================================================="
echo "STORE CLEANUP"
echo "=========================================================================="
echo "  store : $S"
echo "  mode  : $([ "$APPLY" -eq 1 ] && echo 'APPLY, deletions are permanent' || echo 'plan only, nothing is deleted')"
echo
echo "  counting files, this walks the whole store and takes a minute"
n_before=$(find "$S" -type f 2>/dev/null | wc -l)
b_before=$(bsz "$S")
echo "  before: $n_before files, $(numfmt --to=iec "$b_before" 2>/dev/null || echo "$b_before B")"

RUN=$(find "$W/05_orthology" -maxdepth 1 -type d -name 'run_*' 2>/dev/null | head -1)
RES=""
[ -n "$RUN" ] && RES=$(find "$RUN" -maxdepth 1 -type d -name 'Results_*' 2>/dev/null | head -1)

# ===================================================================== rescue
say "RESCUE LIST  —  copied into provenance/ and verified before any delete"

if [ -n "$RES" ]; then
    echo "  from the OrthoFinder run, into provenance/orthofinder_20261007/:"
    for k in Citation.txt Log.txt Comparative_Genomics_Statistics Orthogroups \
             Gene_Duplication_Events Species_Tree \
             Phylogenetic_Hierarchical_Orthogroups; do
        if [ -e "$RES/$k" ]; then
            printf "    %-40s %7s  %5s file(s)\n" "$k" "$(hsz "$RES/$k")" "$(nfil "$RES/$k")"
            resc "$RES/$k" "orthofinder_20261007"
        else
            printf "    %-40s absent\n" "$k"
        fi
    done
    echo
    echo "    Orthogroups_UnassignedGenes.tsv is the one that matters. It is the"
    echo "    evidence that MCL left all three taro PSY candidates unclustered."
    echo "    A README recording the nine-species set is written beside it."
fi

echo
echo "  early search tables, into provenance/early_searches/:"
for f in "$W/05_orthology/psy_check"/*.tsv "$W/05_orthology/pathway_homology"/*.tsv; do
    [ -e "$f" ] || continue
    printf "    %-40s %7s\n" "${f#$W/05_orthology/}" "$(hsz "$f")"
    resc "$f" "early_searches"
done
echo
echo "    LOGBOOK.md quotes identity and e-values out of these. The rewritten"
echo "    03_homology.py changes the coverage metric, so it will not reproduce"
echo "    the same numbers and these are the only source for them."

echo
echo "  GFF3 provenance, into provenance/gff_checksums/:"
found_md5=0
for d in "$W/05_orthology"/*_gff; do
    [ -d "$d" ] || continue
    for f in "$d"/md5sum.txt "$d"/*.md5 "$d"/*checksum*; do
        [ -e "$f" ] || continue
        tag=$(basename "$d" _gff)
        printf "    %-40s %7s\n" "$tag/$(basename "$f")" "$(hsz "$f")"
        resc "$f" "gff_checksums/$tag"
        found_md5=1
    done
done
[ "$found_md5" -eq 0 ] && echo "    none found"
echo
echo "    aratha_gff holds nothing but md5sum.txt, and that checksum is the"
echo "    only record of which Arabidopsis release produced"
echo "    aratha_protein2gene.tsv, which cannot be rebuilt without the GFF3."

# ======================================================================= plan
say "1. the OrthoFinder run"
if [ -n "$RUN" ]; then
    plan "$RUN" "OrthoFinder is not part of the method. Nine-species run, not reproducible from current inputs. Summaries rescued"
else
    echo "  no run_* directory; already removed"
fi

say "2. superseded search directories, minus the tables rescued above"
for f in "$W/05_orthology/pathway_homology"/*.dmnd "$W/05_orthology/pathway_homology"/*.faa; do
    plan "$f" "regenerable; superseded by work/03_homology"
done
for f in "$W/05_orthology/psy_check"/*.dmnd "$W/05_orthology/psy_check"/*.faa; do
    plan "$f" "regenerable; the search tables are rescued"
done
plan "$W/05_orthology/trees/all.faa"  "orphaned combined database; 04_families.sh builds its own"
plan "$W/05_orthology/trees/all.dmnd" "orphaned combined database; 04_families.sh builds its own"

say "3. stale caches a later step would silently reuse"
echo "  04_families.sh decides whether to rebuild with 'if [ ! -f all.dmnd ]'."
echo "  A database built from the old nine-species set would be reused against"
echo "  the current seventeen, and the printed species count comes from the"
echo "  variable rather than the database, so nothing would say so. Same shape"
echo "  as the KEGG query returning zero and the single-tip rooting returning"
echo "  all zeros. Removing the cache forces a rebuild; 04_ also needs a check"
echo "  that compares the cache against the species list."
echo
plan "$W/04_families/all.dmnd" "stale database; 04_ reuses it on existence alone"
plan "$W/04_families/all.faa"  "stale concatenation from the nine-species set"
plan "$W/03_homology/q.faa"    "intermediate of the replaced 03_homology.py"
plan "$W/03_homology/t.faa"    "intermediate of the replaced 03_homology.py"

say "4. outgroup GFF3 files"
echo "  The current pipeline reads only the taro GFF3. These were fetched to"
echo "  map rice and sorghum PSY genes to chromosomes for the subgroup"
echo "  question, which is recorded as not determinable. Re-downloadable with"
echo "  'datasets download genome accession <ACC> --include gff3', though a"
echo "  later NCBI release may differ, which is why the checksums are rescued."
echo
for d in "$W/05_orthology"/*_gff; do
    [ -d "$d" ] || continue
    plan "$d" "outgroup GFF3, not read by the current pipeline; checksum rescued"
done

say "5. duplicate proteome copies, checked by sequence content"
echo "  Each proteome may exist as *.primary.faa and again as"
echo "  orthofinder_input/*.fa with different headers. Comparing a checksum of"
echo "  the sorted sequences, so a file that is not actually a duplicate is"
echo "  kept. This reads every pair and takes a moment."
echo
for f in "$P"/*.primary.faa "$P"/colesc.tagged.faa; do
    [ -e "$f" ] || continue
    b=$(basename "$f"); b=${b%.primary.faa}; b=${b%.tagged.faa}
    alt="$P/orthofinder_input/$b.fa"
    if [ ! -s "$alt" ]; then
        printf "  KEEP    %-26s no orthofinder_input/%s.fa\n" "$(basename "$f")" "$b"
        continue
    fi
    n1=$(grep -c '^>' "$f"); n2=$(grep -c '^>' "$alt")
    s1=$(seqsum "$f");       s2=$(seqsum "$alt")
    if [ "$s1" = "$s2" ]; then
        printf "  MATCH   %-26s %6s seqs, identical sequence set\n" "$(basename "$f")" "$n1"
        plan "$f" "sequence-identical duplicate of orthofinder_input/$b.fa (md5 ${s1:0:8}, $n1 seqs)"
    else
        printf "  DIFFER  %-26s %6s vs %6s seqs — KEPT\n" "$(basename "$f")" "$n1" "$n2"
    fi
done

say "6. download leftovers"
for z in "$P"/*.zip "$P"/*_x "$P"/.figshare_files.json; do
    plan "$z" "download leftover, already extracted"
done

# ===================================================================== report
say "WHAT WOULD BE REMOVED"
total=0
printf "  %-9s %9s  %s\n" "size" "files" "path"
rule
for row in "${PLAN[@]}"; do
    p=${row%%$'\t'*}; why=${row#*$'\t'}
    b=$(bsz "$p"); total=$((total + b))
    printf "  %-9s %9s  %s\n" "$(hsz "$p")" "$(nfil "$p")" "${p#$S/}"
    printf "  %-9s %9s      %s\n" "" "" "$why"
done
rule
echo "  reclaimable : $(numfmt --to=iec "$total" 2>/dev/null || echo "$total B")"
echo "  entries     : ${#PLAN[@]}"
echo "  rescued     : ${#RESCUE[@]} items into provenance/"

if [ "$APPLY" -eq 0 ]; then
    echo
    echo "=========================================================================="
    echo "  Nothing was deleted and nothing was copied."
    echo
    echo "  Read both lists. Anything on the removal list that you want kept has"
    echo "  to be said now, because --apply is not reversible."
    echo
    echo "    bash scripts/cleanup_store.sh --apply"
    echo "=========================================================================="
    exit 0
fi

# ==================================================================== phase 1
say "PHASE 1  —  rescue"
for row in "${RESCUE[@]}"; do
    src=${row%%$'\t'*}; sub=${row#*$'\t'}
    mkdir -p "$PROV/$sub"
    cp -r "$src" "$PROV/$sub/"
    printf "  copied  %-54s -> provenance/%s/\n" "$(basename "$src")" "$sub"
done

if [ -n "$RES" ]; then
    cat > "$PROV/orthofinder_20261007/README.md" <<'RM'
# OrthoFinder run, 7 October 2026

Rescued from `work/05_orthology/run_20261007_2049/` before that directory was
deleted. The per-orthogroup sequences, alignments and trees were not kept:
113,955 files of intermediate product, none of it read by any step of the
method.

## The species set was NINE, not seventeen

    ambtri  aratha  colesc  musacu  nelnuc  orysat  soltub  zeamay  zosmar

Absent, because they were added to the project after this run:

    anacom  aspoff  bradis  elagui  manesc  phodac  setita  sorbic

Read every table here as a nine-species result. In particular, the clustering
failure recorded below was observed on nine species. Whether orthogroup
clustering would also fail on the final seventeen is untested, and that
limitation should be stated rather than glossed if the question is ever raised.

## What this records

`Orthogroups/Orthogroups_UnassignedGenes.tsv` holds all three taro phytoene
synthase candidates. MCL placed none of them in an orthogroup, and 6.6% of
genes were unassigned overall. That is the reason the method uses targeted
reciprocal best hit with an explicit outparalog instead of orthogroup
clustering, and this file is the evidence for it.

## Redoing rather than regenerating

Re-running the archived script now would cluster seventeen proteomes and give a
different result, so this output cannot be reproduced:

    bash ~/Code/taro/scripts/archive/07_orthofinder.sh
RM
    echo "  wrote   provenance/orthofinder_20261007/README.md"
fi

miss=0
for row in "${RESCUE[@]}"; do
    src=${row%%$'\t'*}; sub=${row#*$'\t'}
    [ -e "$PROV/$sub/$(basename "$src")" ] || { echo "  RESCUE FAILED: $src"; miss=1; }
done
if [ "$miss" -eq 1 ]; then
    echo
    echo "  Aborting. The audit trail is not safely copied, so nothing is"
    echo "  deleted. Fix the copy and run again."
    exit 1
fi
echo
echo "  verified: $(nfil "$PROV") files in provenance/, $(hsz "$PROV")"

# ==================================================================== phase 2
say "PHASE 2  —  delete"
: > "$TMP"
t0=$SECONDS
for row in "${PLAN[@]}"; do
    p=${row%%$'\t'*}; why=${row#*$'\t'}
    [ -e "$p" ] || continue
    sz=$(hsz "$p"); n=$(nfil "$p")
    printf "  removing %-50s %7s  %6s files ... " "${p#$S/}" "$sz" "$n"
    rm -rf "$p"
    printf "done  [%ds]\n" "$((SECONDS - t0))"
    printf '| `%s` | %s | %s | %s |\n' "${p#$S/}" "$sz" "$n" "$why" >> "$TMP"
done

say "RESULT"
n_after=$(find "$S" -type f 2>/dev/null | wc -l)
b_after=$(bsz "$S")
echo "  files : $n_before -> $n_after   (removed $((n_before - n_after)))"
echo "  size  : $(numfmt --to=iec "$b_before") -> $(numfmt --to=iec "$b_after")"
echo "  freed : $(numfmt --to=iec $((b_before - b_after)))"
echo "  took  : $((SECONDS - t0))s"

{
    echo "# Store cleanup"
    echo
    echo "Applied $(date -Iseconds)."
    echo
    echo "Files: $n_before to $n_after. Size: $(numfmt --to=iec "$b_before") to $(numfmt --to=iec "$b_after")."
    echo
    echo "Everything removed was unread by the current pipeline. Findings that"
    echo "rested on any of it are in \`~/Code/taro/LOGBOOK.md\`, in"
    echo "\`results/tables/archive/\`, or in \`provenance/\` beside this file."
    echo
    echo "## Removed"
    echo
    echo "| Path | Size | Files | Why |"
    echo "|---|---|---|---|"
    cat "$TMP"
    echo
    echo "## Kept in provenance/"
    echo
    echo "Copied and verified before anything was deleted."
    echo
    find "$PROV" -type f -printf '- `%P`\n' 2>/dev/null | sort
    echo
    echo "The OrthoFinder run used nine species, not the seventeen the project"
    echo "now uses, so it cannot be regenerated, only redone differently."
    echo "\`provenance/orthofinder_20261007/README.md\` states this."
    echo
    echo "## Regenerating"
    echo
    echo '```bash'
    echo "# an outgroup GFF3; compare against the rescued checksum"
    echo "datasets download genome accession <ACC> --include gff3"
    echo
    echo "# the 04_ search database rebuilds itself on the next run"
    echo "bash ~/Code/taro/scripts/04_families.sh"
    echo
    echo "# a proteome copy, if an archived script needs one"
    echo "cp refs/proteomes/orthofinder_input/<tag>.fa refs/proteomes/<tag>.primary.faa"
    echo '```'
} > "$MAN"
rm -f "$TMP"
echo
echo "  manifest: $MAN"
echo
echo "  Copy the Removed table into LOGBOOK.md so the record is in git."
