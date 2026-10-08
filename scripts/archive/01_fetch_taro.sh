#!/usr/bin/env bash
# 01_fetch_taro.sh — taro proteome + annotation from Figshare.
# Sun et al. 2026, Sci Data 13:802 — doi 10.6084/m9.figshare.29917034
set -euo pipefail
: "${TARO_REFS:?source envs/activate.sh first}"

DEST="$TARO_REFS/proteomes"
ART=29917034
mkdir -p "$DEST"

for f in Colocasia_esculenta.Genome.V1.pep Colocasia_esculenta.Genome.V1.gff3; do
    [ -s "$DEST/$f" ] && { echo "have  $f"; continue; }
    url=$(curl -sS "https://api.figshare.com/v2/articles/$ART/files" \
          | python -c "import json,sys;print(next(x['download_url'] for x in json.load(sys.stdin) if x['name']=='$f'))")
    echo "get   $f"
    curl -L --fail --progress-bar -o "$DEST/$f" "$url"
done

echo
echo "proteins: $(grep -c '^>' "$DEST/Colocasia_esculenta.Genome.V1.pep")  (expect ~28253)"
ls -lh "$DEST"
