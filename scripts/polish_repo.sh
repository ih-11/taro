#!/usr/bin/env bash
# polish_repo.sh — archive what is superseded, install the rewritten docs.
#
# Supersedes by `git mv`, never `rm`: history survives and the move is
# reversible with one command. Nothing in results/tables/*.tsv that a current
# script wrote is touched.
#
# Run from ~/Code/taro. Plans by default; --apply does it.
set -euo pipefail
: "${TARO_CODE:?source envs/activate.sh first}"
cd "$TARO_CODE"

DL="/mnt/c/Users/ha-ibnu/Downloads"
APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

rule () { printf '%.0s-' {1..74}; printf "\n"; }
say  () { printf "\n%s\n" "$1"; rule; }

# newest match wins: a file delivered twice lands as "name (1).ext"
newest () {
    local pat="$1" f
    f=$(ls -t "$DL"/"${pat%.*}"*."${pat##*.}" 2>/dev/null | head -1 || true)
    [ -n "$f" ] && printf '%s' "$f"
}

install_file () {   # source-basename, destination path
    local src dst="$2"
    src=$(newest "$1")
    if [ -z "$src" ]; then
        printf "  MISSING in Downloads: %-34s -> %s\n" "$1" "$dst"
        return 1
    fi
    printf "  %-38s -> %s\n" "$(basename "$src")" "$dst"
    if [ "$APPLY" -eq 1 ]; then
        mkdir -p "$(dirname "$dst")"
        rm -f "$dst"
        cp "$src" "$dst"
        sed -i 's/\r$//' "$dst"
        # NTFS through WSL hands every file 755. A document is not executable.
        case "$dst" in
            *.sh) chmod 755 "$dst" ;;
            *)    chmod 644 "$dst" ;;
        esac
    fi
}

archive () {        # path, destination directory, reason
    local p="$1" d="$2" why="$3"
    if [ ! -e "$p" ]; then
        printf "  already gone: %s\n" "$p"
        return 0
    fi
    printf "  %-42s -> %s/\n" "$p" "$d"
    printf "  %-42s    %s\n" "" "$why"
    if [ "$APPLY" -eq 1 ]; then
        mkdir -p "$d"
        git mv -f "$p" "$d/" 2>/dev/null || mv -f "$p" "$d/"
    fi
}

echo "=========================================================================="
echo "REPOSITORY POLISH"
echo "  mode: $([ "$APPLY" -eq 1 ] && echo 'APPLY' || echo 'plan only')"
echo "=========================================================================="

say "1. superseded, moved to archive"
archive docs/analysis_plan.md docs/archive \
  "the original plan; METHODS_part1.md is what the method became"
archive scripts/patch_schema.sh scripts/archive \
  "string-replaced columns in scripts since rewritten; would corrupt them now"
archive results/tables/kegg_only_loci.tsv results/tables/archive \
  "the file 02b_ overwrote in place; raw and reviewed now have one writer each"
archive results/tables/pathway_homology.tsv results/tables/archive \
  "coverage over 100%, no query coordinates, sub-threshold hits discarded"
archive results/tables/scope_exclusions.tsv results/tables/archive \
  "early scope list, before the KEGG cross-check"

say "2. rewritten documents"
ok=0
install_file README.md               README.md                              || ok=1
install_file ANALYSIS_PLAN_part2.md  docs/ANALYSIS_PLAN_part2.md            || ok=1
install_file scripts_README.md       scripts/README.md                      || ok=1
install_file scripts_archive_README.md scripts/archive/README.md            || ok=1
install_file notebooks_README.md     notebooks/README.md                    || ok=1
install_file tables_archive_README.md results/tables/archive/README.md      || ok=1

say "4. result"
if [ "$APPLY" -eq 0 ]; then
    echo "  Nothing moved, nothing copied."
    [ "$ok" -ne 0 ] && echo "  Some files are missing from Downloads; fix that first."
    echo
    echo "    bash scripts/polish_repo.sh --apply"
    exit 0
fi

python scripts/06_readme_block.py
echo
echo "  git status:"
git status --short | sed 's/^/    /'
echo
echo "  Read the diff, then:  ga .   gp"
