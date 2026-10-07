#!/usr/bin/env bash
# 07_orthofinder.sh — orthogroups across taro + 8 outgroups.
#
# Species set after Yin et al. 2021 (Mol Ecol Resour 21:68), minus the four
# Araceae relatives, none of which has a protein set at NCBI. Spirodela
# polyrhiza — taro's closest sequenced relative, ~73 MYA — is unavailable;
# state as a limitation.
#
# One protein per gene throughout (06_), so orthogroup sizes are gene counts,
# not isoform counts. That distinction is the whole point: it decides whether
# taro appears to have one PSY or several.
#
# Runtime ~2-4 h on 24 threads. Safe to leave.
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

IN="$TARO_REFS/proteomes/orthofinder_input"
OUT="$TARO_WORK/05_orthology"
LOG="$TARO_LOGS/orthofinder_$(date +%Y%m%d_%H%M).log"
mkdir -p "$OUT" "$TARO_LOGS"

echo "input   $IN"
grep -c '^>' "$IN"/*.fa | sed 's|.*/|        |'
echo "output  $OUT"
echo "log     $LOG"
echo "threads ${TARO_THREADS:-24}"
echo

orthofinder \
    -f "$IN" \
    -o "$OUT/run_$(date +%Y%m%d_%H%M)" \
    -t "${TARO_THREADS:-24}" \
    -a 8 \
    -S diamond \
    2>&1 | tee "$LOG"

echo
echo "done. results:"
find "$OUT" -name 'Statistics_Overall.tsv' -newermt '-1 day' | tail -1
