#!/usr/bin/env bash
# 00_fetch_proteomes.sh — Part 1 inputs only.
#
# Part 1 (carotenoid pathway gene set, PSY copy number) is orthology on
# proteomes. It needs NO genome assembly: a few hundred MB of protein FASTA,
# not 7 Gb of references. Assemblies come later, for Parts 2 and 3.
#
#   probe (default)  resolve every source, report availability and size,
#                    download nothing
#   --download       fetch what probe resolved
#
# Usage:
#   source ~/Code/taro/envs/activate.sh
#   bash scripts/00_fetch_proteomes.sh
#   bash scripts/00_fetch_proteomes.sh --download

set -uo pipefail

: "${TARO_REFS:?run: source ~/Code/taro/envs/activate.sh}"
PROT="$TARO_REFS/proteomes"
MAN="$PROT/manifest.tsv"
DL=0; [ "${1:-}" = "--download" ] && DL=1

G='\033[0;32m'; R='\033[0;31m'; Y='\033[0;33m'; B='\033[0;34m'; N='\033[0m'
hdr () { printf "\n${B}=== %s ===${N}\n" "$1"; }
ok  () { printf "  ${G}[ OK ]${N}  %s\n" "$1"; }
wrn () { printf "  ${Y}[MANUAL]${N} %s\n" "$1"; }
bad () { printf "  ${R}[FAIL]${N}  %s\n" "$1"; }
inf () { printf "          %s\n" "$1"; }

mkdir -p "$PROT"
[ -f "$MAN" ] || printf "tag\tspecies\tsource\turl_or_accession\tfile\tstatus\tdate\n" > "$MAN"

note_manifest () { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$1" "$2" "$3" "$4" "$5" "$6" "$(date -I)" >> "$MAN"; }

printf "${B}Part 1 inputs — proteomes for orthology${N}\n"
printf "mode: %s    dest: %s\n" "$([ $DL -eq 1 ] && echo DOWNLOAD || echo PROBE)" "$PROT"

# ============================================================ 1. taro proteome
hdr "Taro proteome (the target)"
inf "Sun et al. 2026, Sci Data 13:802 — Figshare 10.6084/m9.figshare.29917034"

FIG_ID=29917034
FIG_API="https://api.figshare.com/v2/articles/${FIG_ID}/files"
FIG_JSON="$PROT/.figshare_files.json"

if curl -sSf -m 30 "$FIG_API" -o "$FIG_JSON" 2>/dev/null; then
    ok "Figshare API reachable"
    python - "$FIG_JSON" <<'PY'
import json,sys
try:
    files=json.load(open(sys.argv[1]))
except Exception as e:
    print(f"          could not parse: {e}"); sys.exit()
if not files:
    print("          article returned no files"); sys.exit()
for f in files:
    mb=f.get('size',0)/1048576
    print(f"          {f['name']:55s} {mb:8.1f} MB")
PY
    inf ""
    inf "Want: Colocasia_esculenta.Genome.V1.pep   (proteome -> OrthoFinder)"
    inf "Also: Colocasia_esculenta.Genome.V1.gff3  (Part 3 later)"

    if [ $DL -eq 1 ]; then
        python - "$FIG_JSON" "$PROT" <<'PY'
import json,sys,subprocess,os
files=json.load(open(sys.argv[1])); dest=sys.argv[2]
want=('.pep','.gff3','.cds')
for f in files:
    if not f['name'].lower().endswith(want): continue
    out=os.path.join(dest,f['name'])
    if os.path.exists(out):
        print(f"  skip {f['name']} (exists)"); continue
    print(f"  get  {f['name']} ({f.get('size',0)/1048576:.0f} MB)")
    subprocess.run(["aria2c","-c","-x4","-s4","--summary-interval=0",
                    "-d",dest,"-o",f['name'],f['download_url']])
    md5=f.get('computed_md5') or f.get('supplied_md5')
    if md5:
        import hashlib
        h=hashlib.md5()
        with open(out,'rb') as fh:
            for b in iter(lambda: fh.read(1<<20), b''): h.update(b)
        print(f"  md5  {'OK' if h.hexdigest()==md5 else 'MISMATCH'}  {f['name']}")
PY
        note_manifest asm2026_pep "Colocasia esculenta" figshare "$FIG_ID" "Colocasia_esculenta.Genome.V1.pep" downloaded
    fi
else
    bad "Figshare API unreachable"
    inf "fallback: open https://doi.org/10.6084/m9.figshare.29917034 and download .pep manually"
    note_manifest asm2026_pep "Colocasia esculenta" figshare "$FIG_ID" "Colocasia_esculenta.Genome.V1.pep" FAILED
fi

# ============================================================ 2. Araceae outgroups
hdr "Araceae outgroup proteomes"
inf "Same four species Sun et al. used for homology-based prediction,"
inf "so the comparison is anchored to the published annotation's assumptions."

# species | tag | note
SPP=(
"Spirodela polyrhiza|spipol|An et al. 2019 PNAS — closest relative, ~73 MYA"
"Pistia stratiotes|pisstr|Qian et al. 2022 Mol Ecol Resour 22:2732"
"Amorphophallus konjac|amokon|Gao et al. 2022 Comput Struct Biotechnol J 20:1002"
"Zantedeschia elliottiana|zanell|Wang et al. 2023 Sci Data 10:605"
)

for entry in "${SPP[@]}"; do
    IFS='|' read -r sp tag ref <<< "$entry"
    printf "\n  ${B}%s${N}  (%s)\n" "$sp" "$tag"
    inf "$ref"

    SUM=$(datasets summary genome taxon "$sp" --assembly-source all 2>/dev/null)
    CNT=$(echo "$SUM" | python -c "import json,sys;d=json.load(sys.stdin);print(d.get('total_count',0))" 2>/dev/null || echo 0)

    if [ "${CNT:-0}" -gt 0 ]; then
        ACC=$(echo "$SUM" | python - <<'PY'
import json,sys
d=json.load(sys.stdin)
rs=d.get('reports',[])
# prefer RefSeq (GCF) and annotated assemblies
rs.sort(key=lambda r:(not r['accession'].startswith('GCF'),
                      'annotation_info' not in r))
r=rs[0]
print(r['accession'], r.get('organism',{}).get('organism_name','?'),
      'ANNOTATED' if 'annotation_info' in r else 'NO_ANNOTATION', sep='\t')
PY
)
        read -r acc name anno <<< "$ACC"
        if [ "$anno" = "ANNOTATED" ]; then
            ok "$acc  ($name) — has annotation, protein FASTA available"
            if [ $DL -eq 1 ]; then
                ( cd "$PROT" && datasets download genome accession "$acc" \
                    --include protein --filename "${tag}.zip" --no-progressbar \
                  && unzip -oq "${tag}.zip" -d "${tag}_tmp" \
                  && find "${tag}_tmp" -name 'protein.faa' -exec mv {} "${tag}.faa" \; \
                  && rm -rf "${tag}_tmp" "${tag}.zip" \
                  && echo "          -> ${tag}.faa  $(grep -c '^>' "${tag}.faa") proteins" )
                note_manifest "$tag" "$sp" ncbi "$acc" "${tag}.faa" downloaded
            fi
        else
            wrn "$acc  ($name) — assembly present but NO protein set at NCBI"
            inf "get the .pep from the paper's own repository (see ref above)"
            note_manifest "$tag" "$sp" manual "$acc" "${tag}.faa" MANUAL
        fi
    else
        wrn "not found at NCBI under this taxon name"
        inf "likely in NGDC/CNCB or journal supplement — see ref above"
        note_manifest "$tag" "$sp" manual - "${tag}.faa" MANUAL
    fi
done

# ============================================================ 3. reference proteins
hdr "Reference PSY / OR proteins (the anchors)"
inf "Small set. These define what we are searching for."

cat > "$PROT/reference_proteins.txt" <<'EOF'
# Query anchors for the carotenoid pathway search.
# Fetch from UniProt or NCBI protein, concatenate into reference_anchors.faa
#
# PHYTOENE SYNTHASE — the engineering target
Manihot esculenta   MePSY1   Manes.02G081700   # cassava PSY1, the Carvita25 allele
Manihot esculenta   MePSY2   Manes.01G124200   # cassava PSY2
Arabidopsis thaliana  AtPSY   AT5G17230
Oryza sativa          OsPSY1  Os06g0729000     # monocot reference
Zea mays              ZmPSY1  GRMZM2G300348
#
# ORANGE — post-translational control of PSY stability
Arabidopsis thaliana  AtOR        AT5G61670
Arabidopsis thaliana  AtOR-like   AT5G06130
Ipomoea batatas       IbOr        AFQ31814
Cucumis melo          CmOr        XP_008452586
Brassica oleracea     BoOR        ACR61567
#
# downstream / degradation
Arabidopsis thaliana  AtPDS    AT4G14210
Arabidopsis thaliana  AtLCYB   AT3G10230
Arabidopsis thaliana  AtLCYE   AT5G57030
Arabidopsis thaliana  AtCCD4   AT4G19170
Arabidopsis thaliana  AtCCD1   AT3G63520
EOF
ok "wrote reference_proteins.txt"
inf "cassava: https://phytozome-next.jgi.doe.gov  (Mesculenta v8.1)"
inf "Arabidopsis: https://www.arabidopsis.org  or UniProt"

# ============================================================ summary
hdr "Summary"
if [ $DL -eq 0 ]; then
    printf "  Probe only — nothing downloaded.\n\n"
    printf "  To fetch:  bash scripts/00_fetch_proteomes.sh --download\n\n"
    printf "  Anything marked ${Y}MANUAL${N} needs retrieving by hand; the reference\n"
    printf "  is printed above each. Part 1 can proceed with taro + 2 outgroups if\n"
    printf "  one or two prove hard to get — more outgroups sharpen the orthogroups\n"
    printf "  but are not strictly required to count taro PSY paralogs.\n\n"
else
    printf "  Downloaded into %s\n" "$PROT"
    ls -lh "$PROT"/*.faa "$PROT"/*.pep 2>/dev/null || true
    printf "\n  manifest: %s\n\n" "$MAN"
    printf "  Next: bash scripts/01_orthology.sh\n\n"
fi
