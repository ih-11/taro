#!/usr/bin/env bash
# extract_provenance.sh — rebuild the record of where the 17 proteomes came
# from, BEFORE the cleanup deletes the last on-disk evidence of it.
#
# WHY THIS HAS TO RUN FIRST
#
# manifest.tsv was supposed to be the provenance record. It is not. It has the
# right header:
#
#   tag  species  source  url_or_accession  file  status  date
#
# and four data rows, for Spirodela polyrhiza, Pistia stratiotes,
# Amorphophallus konjac and Zantedeschia elliottiana. None of those four is in
# the project. All are marked "manual", so the file is a list of Araceae and
# Alismatales genomes that were WANTED and never fetched.
#
# So not one of the seventeen proteomes the inventory is built on has a
# recorded accession, assembly version or annotation release. That is a
# reproducibility hole in Part 1 itself, not a cleanup problem: a gene family
# count depends on which annotation release was counted, and right now nothing
# on disk says which that was.
#
# It is recoverable from two places, and both are on the cleanup's delete list.
#
#   1. The nine NCBI genomic.gff files. An NCBI GFF3 carries its own assembly
#      accession in the header lines beginning #!genome-build-accession and
#      #!annotation-source. Those nine directories are 1.4 GB and are first in
#      line to be deleted, and the header is the only thing in them worth
#      keeping.
#
#   2. The archived download scripts, which have the accessions written into
#      them. Those are in git, so they are not at risk, but reading them lets
#      the two sources be cross-checked against each other.
#
#   3. .figshare_files.json, 4 KB, which is the API response for the taro
#      figshare download and therefore the only record of which figshare
#      article the Sun et al. proteome came from. It is on the delete list as a
#      "download leftover". It is copied out here.
#
# This script reads, copies and writes. It deletes nothing.
set -euo pipefail
: "${TARO_STORE:?source envs/activate.sh first}"
: "${TARO_REFS:?}"
: "${TARO_WORK:?}"
: "${TARO_CODE:?}"

P="$TARO_REFS/proteomes"
G="$TARO_WORK/05_orthology"
OUT="$TARO_STORE/provenance/sources"
mkdir -p "$OUT"

rule () { printf '%.0s-' {1..92}; printf "\n"; }
say  () { printf "\n%s\n" "$1"; rule; }

declare -A SP=(
  [ambtri]="Amborella trichopoda"    [anacom]="Ananas comosus"
  [aratha]="Arabidopsis thaliana"    [aspoff]="Asparagus officinalis"
  [bradis]="Brachypodium distachyon" [colesc]="Colocasia esculenta"
  [elagui]="Elaeis guineensis"       [manesc]="Manihot esculenta"
  [musacu]="Musa acuminata"          [nelnuc]="Nelumbo nucifera"
  [orysat]="Oryza sativa"            [phodac]="Phoenix dactylifera"
  [setita]="Setaria italica"         [soltub]="Solanum tuberosum"
  [sorbic]="Sorghum bicolor"         [zeamay]="Zea mays"
  [zosmar]="Zostera marina"
)

seqsum () {
    awk '/^>/{if(s!="")print s; s=""; next}{s=s $0}END{if(s!="")print s}' "$1" \
        | sort | md5sum | cut -d' ' -f1
}

echo "=============================================================================="
echo "PROVENANCE RECOVERY"
echo "=============================================================================="

# ----------------------------------------------------------------- 1. the GFF3
say "1. assembly accessions from the NCBI GFF3 headers"
declare -A ACC BUILD ANNOT
found=0
for d in "$G"/*_gff; do
    [ -d "$d" ] || continue
    tag=$(basename "$d" _gff)
    f=$(find "$d" -name 'genomic.gff' -o -name '*.gff' -o -name '*.gff3' 2>/dev/null | head -1)
    if [ -z "$f" ] || [ ! -s "$f" ]; then
        printf "  %-8s %-26s no readable GFF3 in this directory\n" "$tag" "${SP[$tag]:-?}"
        continue
    fi
    hdr=$(head -40 "$f")
    a=$(printf '%s\n' "$hdr" | sed -n 's/^#!genome-build-accession *\(.*\)$/\1/p' | head -1)
    b=$(printf '%s\n' "$hdr" | sed -n 's/^#!genome-build *\(.*\)$/\1/p' | head -1)
    n=$(printf '%s\n' "$hdr" | sed -n 's/^#!annotation-source *\(.*\)$/\1/p' | head -1)
    ACC[$tag]="${a#NCBI_Assembly:}"
    BUILD[$tag]="$b"
    ANNOT[$tag]="$n"
    printf "  %-8s %-26s %-20s %s\n" "$tag" "${SP[$tag]:-?}" "${ACC[$tag]:-?}" "${BUILD[$tag]:-?}"
    [ -n "${ANNOT[$tag]}" ] && printf "  %-8s %-26s   annotation: %s\n" "" "" "${ANNOT[$tag]}"
    # keep the whole header, it is a few hundred bytes
    printf '%s\n' "$hdr" | grep '^#' > "$OUT/${tag}_gff_header.txt" || true
    found=$((found + 1))
done
echo
echo "  headers saved for $found species into provenance/sources/*_gff_header.txt"

# -------------------------------------------------- 2. the archived scripts
say "2. accessions written into the archived download scripts"
if ls "$TARO_CODE"/scripts/archive/*.sh >/dev/null 2>&1; then
    grep -hoE 'GC[AF]_[0-9]+\.[0-9]+' "$TARO_CODE"/scripts/archive/*.sh 2>/dev/null \
        | sort -u > "$OUT/accessions_from_scripts.txt" || true
    if [ -s "$OUT/accessions_from_scripts.txt" ]; then
        echo "  found $(wc -l < "$OUT/accessions_from_scripts.txt") distinct accessions:"
        sed 's/^/    /' "$OUT/accessions_from_scripts.txt"
    else
        echo "  no GCA_/GCF_ accessions in scripts/archive/"
    fi
    echo
    echo "  lines that mention a download, for context:"
    grep -nE 'datasets download|GC[AF]_[0-9]|figshare|ndownloader|ftp\.' \
        "$TARO_CODE"/scripts/archive/*.sh 2>/dev/null \
        | sed 's/^/    /' | cut -c1-150 | head -40 > "$OUT/download_lines.txt" || true
    cat "$OUT/download_lines.txt" 2>/dev/null || echo "    none"
else
    echo "  scripts/archive/ is empty"
fi

# ------------------------------------------------------------- 3. figshare
say "3. the taro figshare record"
if [ -s "$P/.figshare_files.json" ]; then
    cp "$P/.figshare_files.json" "$OUT/taro_figshare_files.json"
    echo "  copied .figshare_files.json -> provenance/sources/taro_figshare_files.json"
    python3 - "$P/.figshare_files.json" <<'PY' || echo "  (not parseable as JSON, kept verbatim)"
import json, sys
d = json.load(open(sys.argv[1]))
items = d if isinstance(d, list) else d.get("files", [d])
for it in items if isinstance(items, list) else []:
    if isinstance(it, dict):
        print(f"    {it.get('name','?'):<44} {it.get('size','?')} bytes  "
              f"md5 {it.get('computed_md5') or it.get('supplied_md5') or '?'}")
PY
else
    echo "  .figshare_files.json is absent or empty"
fi

# -------------------------------------------------------- 4. rebuilt manifest
say "4. rebuilt manifest"
echo "  computing a sequence checksum per proteome so content drift is"
echo "  detectable from here on. This reads all 17 files and takes a moment."
echo
MF="$OUT/manifest_rebuilt.tsv"
printf 'tag\tspecies\tassembly_accession\tgenome_build\tannotation_source\tn_proteins\tseq_md5\tevidence\n' > "$MF"
printf "  %-8s %-26s %-18s %8s  %s\n" "tag" "species" "accession" "proteins" "evidence"
rule
for tag in $(printf '%s\n' "${!SP[@]}" | sort); do
    fa="$P/orthofinder_input/$tag.fa"
    [ -s "$fa" ] || { printf "  %-8s %-26s MISSING proteome\n" "$tag" "${SP[$tag]}"; continue; }
    n=$(grep -c '^>' "$fa")
    m=$(seqsum "$fa")
    acc="${ACC[$tag]:-}"
    ev="unrecorded"
    [ -n "$acc" ] && ev="NCBI GFF3 header"
    if [ -z "$acc" ] && [ "$tag" = "colesc" ]; then
        acc=""; ev="figshare, see taro_figshare_files.json"
    fi
    printf "  %-8s %-26s %-18s %8s  %s\n" "$tag" "${SP[$tag]}" "${acc:-?}" "$n" "$ev"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$tag" "${SP[$tag]}" "${acc:-NA}" "${BUILD[$tag]:-NA}" \
        "${ANNOT[$tag]:-NA}" "$n" "$m" "$ev" >> "$MF"
done

# ---------------------------------------------------------- 5. mapping file
say "5. the Arabidopsis protein-to-gene mapping"
AP="$G/aratha_protein2gene.tsv"
if [ -s "$AP" ]; then
    np=$(tail -n +2 "$AP" | cut -f1 | sort -u | wc -l)
    ng=$(tail -n +2 "$AP" | cut -f2 | sort -u | wc -l)
    nl=$(tail -n +2 "$AP" | cut -f4 | grep -c . || true)
    nf=$(grep -c '^>' "$P/orthofinder_input/aratha.fa")
    cov=$(awk 'NR==FNR{if($0~/^>/){a=$0;sub(/^>/,"",a);sub(/ .*/,"",a);split(a,p,"|");k[p[length(p)]]=1};next}
               FNR>1{seen[$1]=1} END{c=0; for(x in k) if(x in seen) c++; print c}' \
          "$P/orthofinder_input/aratha.fa" "$AP")
    {
      echo "file: work/05_orthology/aratha_protein2gene.tsv"
      echo "md5: $(md5sum "$AP" | cut -d' ' -f1)"
      echo "proteins: $np"
      echo "genes: $ng"
      echo "rows with an AT locus: $nl"
      echo "aratha.fa sequences: $nf"
      echo "aratha.fa accessions found in the mapping: $cov"
      echo "source: derived from an Arabidopsis GFF3 by scripts/archive/06_rebuild_primary.sh"
      echo "note: the GFF3 is gone and aratha_gff/md5sum.txt is zero bytes, so the"
      echo "      release that produced this file is not recorded. The checksum above"
      echo "      pins the file itself."
    } > "$OUT/aratha_protein2gene.provenance.txt"
    printf "  proteins in mapping          %s\n" "$np"
    printf "  genes in mapping             %s\n" "$ng"
    printf "  aratha.fa sequences          %s\n" "$nf"
    printf "  of those, present in mapping %s  %s\n" "$cov" \
        "$([ "$cov" = "$nf" ] && echo '(complete)' || echo '(INCOMPLETE, investigate)')"
    echo
    echo "  48,265 proteins against 27,562 genes is expected: the mapping covers"
    echo "  every Arabidopsis isoform, while aratha.fa holds one protein per gene."
    echo "  The number that matters is the last line, whether every sequence in"
    echo "  aratha.fa can be looked up."
    echo
    echo "  recorded: provenance/sources/aratha_protein2gene.provenance.txt"
else
    echo "  MISSING: $AP"
fi

say "DONE"
echo "  written to: $OUT"
ls -1 "$OUT" | sed 's/^/    /'
echo
echo "  Next: copy manifest_rebuilt.tsv into the repo so it is in git, then the"
echo "  GFF3 directories are safe to delete."
echo
echo "    cp \"$MF\" \"\$TARO_CODE/docs/proteome_manifest.tsv\""
echo "    bash scripts/cleanup_store.sh --apply"
