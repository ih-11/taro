#!/usr/bin/env python3
"""
02_pathway_set.py — define the pathway scope, and check it for completeness.

Revised after the KEGG cross-check, which did what it was built for: it exposed
a scope and naming problem before any biological search ran.

WHAT THE CHECK FOUND

  Three of our anchors are multi-member Arabidopsis families:
      DXS     AT4G15560  + DXPS1 AT3G21500, DXPS3 AT5G11380    (3)
      NCED3   AT3G14440  + NCED2/5/6/9                          (5)
      GGPPS11 AT4G36810  + GGPS2/3/4/6                          (5)
  No sequences were missed, because a family search seeded by any member pulls
  the rest. What was wrong was the label: "NCED3 copy number" reads as orthologs
  of one gene when the measurement is family size.

  One genuine gap: IDI, isopentenyl diphosphate isomerase. HDR produces IPP and
  DMAPP, IDI interconverts them, GGPS condenses them. Our MEP chain stepped
  straight over it.

  All five anchors KEGG did not return are correct loci. KEGG knows each and
  files it outside ath00906/ath00900, which is expected for ORANGE, ORANGE-like
  and fibrillin, and also holds for CCD1 and NSX (ABA4).

THE SCHEMA CHANGE

  An anchor is a seed for entering a family. It is not the family being
  measured. Conflating the two is what produced "NCED3 copy number". The table
  now separates them:

      pathway  target  anchor_gene  anchor_locus  target_level
      claim_kind  outparalog  scope_reason

TWO ROLES FOR CCD7 AND CCD8

  Biological scope: excluded. They are the strigolactone branch
  (D27 -> CCD7 -> CCD8 -> MAX1 -> LBO1), which is carotenoid-derived but is not
  carotenoid accumulation.

  Phylogenetic reference: included. They sit outside both CCD4 and NCED and are
  the outparalog that makes those family boundaries testable.

  Those are different roles and the record states both.
"""
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
MAP = WORK / "05_orthology" / "aratha_protein2gene.tsv"

# pathway, target, anchor_gene, anchor_locus, target_level, claim_kind,
# outparalog, scope_reason
#
# target_level is "family" throughout: a homology search seeded by one member
# recovers the family, so every count is a family-level count whether or not
# Arabidopsis happens to have one member.
ANCHORS = [
 ("MEP","DXS family","DXS","AT4G15560","family","detection","",
  "3 Arabidopsis members (DXS, DXPS1, DXPS3); no outparalog established, so no copy-number claim"),
 ("MEP","DXR family","DXR","AT5G62790","family","detection","",""),
 ("MEP","MCT family","MCT","AT2G02500","family","detection","",""),
 ("MEP","CMK family","CMK","AT2G26930","family","detection","",""),
 ("MEP","MDS family","MDS","AT1G63970","family","detection","",""),
 ("MEP","HDS family","HDS","AT5G60600","family","detection","",""),
 ("MEP","HDR family","HDR","AT4G34350","family","detection","",
  "no outparalog established, so no copy-number claim"),
 ("MEP","IDI family","IPP1","AT5G16440","family","detection","",
  "added after the KEGG check; HDR makes IPP and DMAPP, IDI interconverts them, "
  "GGPS condenses them, so the chain had a gap"),
 ("backbone","GGPPS family","GGPPS11","AT4G36810","family","copy-number",
  "FPS1:AT5G47770;FPS2:AT4G17190",
  "promoted from detection during pre-analysis scope review: GGPP is the immediate "
  "precursor of phytoene and the family has 5 Arabidopsis members. Outparalog proposed, "
  "not established; tree A tests whether FPS is monophyletic and outside GGPS"),
 ("core","PSY family","PSY","AT5G17230","family","copy-number",
  "SQS1:AT4G34640;SQS2:AT4G34650",""),
 ("core","PDS family","PDS3","AT4G14210","family","detection","",""),
 ("core","Z-ISO family","Z-ISO","AT1G10830","family","detection","",""),
 ("core","ZDS family","ZDS","AT3G04870","family","detection","",""),
 ("core","CRTISO family","CRTISO","AT1G06820","family","detection","",""),
 ("cyclase","LCYB family","LCYB","AT3G10230","family","detection","",""),
 ("cyclase","LCYE family","LCYE","AT5G57030","family","detection","",""),
 ("xanth","CYP97A family","CYP97A3","AT1G31800","family","detection","",""),
 ("xanth","CYP97C family","CYP97C1","AT3G53130","family","detection","",""),
 ("xanth","BCH family","BCH1","AT4G25700","family","detection","",""),
 ("xanth","BCH family 2","BCH2","AT5G52570","family","detection","",""),
 ("xanth","ZEP family","ZEP","AT5G67030","family","detection","",""),
 ("xanth","VDE family","VDE","AT1G08550","family","detection","",""),
 ("xanth","NSX family","NSX","AT1G67080","family","detection","",
  "KEGG names this locus ABA4; it is neoxanthin synthase and the anchor is correct"),
 ("degrade","CCD1 family","CCD1","AT3G63520","family","detection","",""),
 ("degrade","CCD4 family","CCD4","AT4G19170","family","copy-number",
  "CCD7:AT2G44990;CCD8:AT4G32810",
  "outparalogs are the strigolactone branch: excluded from scope, included as "
  "phylogenetic reference"),
 ("degrade","NCED family","NCED3","AT3G14440","family","copy-number",
  "CCD7:AT2G44990;CCD8:AT4G32810",
  "5 Arabidopsis members (NCED2/3/5/6/9); the measurement is family size, not "
  "NCED3 orthologs"),
 ("sink","ORANGE","OR","AT5G61670","gene","presence",""
  ,"acts on PSY protein stability, not a pathway enzyme"),
 ("sink","ORANGE-like","OR-like","AT5G06130","gene","presence","",""),
 ("sink","fibrillin","FBN1a","AT4G04020","gene","presence","",""),
]

# Loci reviewed and deliberately excluded. Recorded so the scope is auditable
# rather than implicit in what is absent.
EXCLUDED = [
 ("strigolactone branch", "D27 AT1G03055; CCD7 AT2G44990; CCD8 AT4G32810; "
  "MAX1/CYP711A1 AT2G26170; LBO1 AT3G21420",
  "carotenoid-derived but not carotenoid accumulation; CCD7 and CCD8 are "
  "retained as phylogenetic outparalogs"),
 ("cytosolic mevalonate", "HMG1 AT1G76490; HMG2 AT2G17370; HMGS AT4G11820; "
  "AACT1 AT5G47720; ACAT2 AT5G48230; MVD1 AT2G38700; MK AT5G27450; "
  "FPS1 AT5G47770; FPS2 AT4G17190",
  "cytosolic MVA, not the plastidial MEP pathway; FPS1/FPS2 are retained as "
  "the proposed GGPPS outparalog"),
 ("abscisic acid catabolism", "CYP707A1 AT4G19230; CYP707A2 AT2G29090; "
  "CYP707A3 AT5G45340; CYP707A4 AT3G19270; ABA2 AT1G52340; AAO3 AT2G27150; "
  "BGLU18 AT1G52400; FLDH AT4G33360",
  "downstream of NCED, outside the inventory's scope"),
 ("protein prenylation", "FTA AT3G59380; FACE2 AT2G36305; FCLY AT5G63910; "
  "STE14A AT5G23320; ATSTE14B AT5G08335; ERA1 AT5G40280; ATSTE24 AT4G01320; "
  "ICME-LIKE1 AT1G26120; ICME-LIKE2 AT3G02410; PCME AT5G15860",
  "protein modification, not isoprenoid biosynthesis"),
 ("dolichol and plastoquinone chains", "all cPT/undecaprenyl loci; "
  "SPS1 AT1G78510; SPS2 AT1G17050",
  "long-chain polyprenyl products, not GGPP"),
]

KEGG = {"ath00906": "carotenoid biosynthesis",
        "ath00900": "terpenoid backbone biosynthesis"}


def kegg_loci(pathway):
    """
    Returns (status, set). Status is 'ok', 'unreachable' or 'empty'.
    An empty response is a failure, not a pathway with no genes: the latter
    does not exist. An earlier version reported an empty result as data, which
    read as KEGG confirming all 28 anchors when the check had not run.
    """
    url = f"https://rest.kegg.jp/link/ath/{pathway}"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            txt = r.read().decode()
    except (urllib.error.URLError, OSError, TimeoutError) as exc:
        return "unreachable", {"_reason": str(exc)}
    out = {p.split("\t")[1][4:].upper() for p in txt.splitlines()
           if len(p.split("\t")) == 2 and p.split("\t")[1].startswith("ath:")}
    return ("ok", out) if out else ("empty", set())


print("=" * 74)
print("PATHWAY SCOPE")
print("=" * 74)
print(f"\n{len(ANCHORS)} targets")

loci_present = {ln.rstrip("\n").split("\t")[3].upper()
                for ln in open(MAP) if not ln.startswith("protein\t")
                and ln.rstrip("\n").split("\t")[3]}

missing = [a for a in ANCHORS if a[3].upper() not in loci_present]
if missing:
    print("\n  anchors that do not resolve in the Arabidopsis proteome:")
    for a in missing:
        print(f"    {a[2]:<10} {a[3]}")
    sys.exit(1)
print("  every anchor resolves in the Arabidopsis proteome")

# outparalogs must resolve too, or the copy-number claim cannot be made
for a in ANCHORS:
    if a[5] != "copy-number":
        continue
    for op in filter(None, a[6].split(";")):
        name, loc = op.split(":")
        if loc.upper() not in loci_present:
            print(f"\n  outparalog {name} ({loc}) for {a[1]} does not resolve")
            sys.exit(1)
print("  every declared outparalog resolves")

print()
print("=" * 74)
print("KEGG COMPLETENESS CHECK")
print("=" * 74)

ours = {a[3].upper() for a in ANCHORS}
kegg_all, failed = set(), []
for pid, name in KEGG.items():
    status, loci = kegg_loci(pid)
    if status == "ok":
        kegg_all |= loci
        print(f"  {pid}  {name}: {len(loci)} loci")
    else:
        print(f"  {pid}  {name}: {status.upper()}")
        failed.append(pid)

if failed:
    print(f"""
  CROSS-CHECK FAILED for {', '.join(failed)}.
  A KEGG pathway cannot contain zero genes, so an empty response means the
  query failed. This exits non-zero rather than reporting the scope as checked.
  Verify by hand:  curl -sS 'https://rest.kegg.jp/link/ath/ath00906' | head
""")
    sys.exit(1)

only_kegg = sorted(kegg_all - ours)
print(f"\n  in both: {len(ours & kegg_all)}   ours only: {len(ours - kegg_all)}"
      f"   KEGG only: {len(only_kegg)}")

# ---------------------------------------------------------------- write
out_dir = CODE / "results" / "tables"
out_dir.mkdir(parents=True, exist_ok=True)

anchors_f = out_dir / "pathway_anchors.tsv"
with open(anchors_f, "w") as fh:
    fh.write("pathway\ttarget\tanchor_gene\tanchor_locus\ttarget_level\t"
             "claim_kind\toutparalog\tin_kegg\tscope_reason\n")
    for p, t, g, l, lvl, k, op, why in ANCHORS:
        fh.write(f"{p}\t{t}\t{g}\t{l}\t{lvl}\t{k}\t{op}\t"
                 f"{'yes' if l.upper() in kegg_all else 'no'}\t{why}\n")

# raw list, never overwritten by the annotation step
raw_f = out_dir / "kegg_only_loci_raw.tsv"
with open(raw_f, "w") as fh:
    fh.write("at_locus\n")
    for loc in only_kegg:
        fh.write(loc + "\n")

excl_f = out_dir / "scope_exclusions.tsv"
with open(excl_f, "w") as fh:
    fh.write("group\tloci\treason\n")
    for grp, loci, why in EXCLUDED:
        fh.write(f"{grp}\t{loci}\t{why}\n")

print()
print("=" * 74)
print("CLAIM KINDS, fixed before any search")
print("=" * 74)
kinds = {}
for a in ANCHORS:
    kinds.setdefault(a[5], []).append(a[1])
for k in ("copy-number", "detection", "presence"):
    v = kinds.get(k, [])
    print(f"\n  {k}  ({len(v)})\n    {', '.join(v)}")

print(f"""
  Every target_level is "family" except the three presence claims, because a
  homology search seeded by one member recovers the family. A count is a
  family-level count whether or not Arabidopsis has one member.

  DXS and HDR are detection, not copy-number. Both are multi-member families
  and neither has an established outparalog, so neither can satisfy the
  phylogeny limb. Saying copy-number in the table while skipping them in 04_
  would make the machine-readable method and the narrative disagree.

  GGPPS is copy-number with FPS1/FPS2 proposed as outparalog. Proposed, not
  established: tree A tests whether FPS is monophyletic and falls outside GGPS,
  and 04_ makes no copy-number claim if it does not. A monophyletic outgroup
  shows FPS sits outside GGPS; it does not show FPS is the NEAREST outparalog,
  and too distant an outgroup draws the family boundary too wide.

written: {anchors_f}
written: {raw_f}       (input to 02b_, never overwritten)
written: {excl_f}      (scope decisions, auditable)
next:    python scripts/02b_annotate_kegg.py
""")
