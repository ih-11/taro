#!/usr/bin/env python3
"""
05_report.py — the master table, with decision rules applied as written.

Copy number requires both limbs:

    family assignment  (phylogeny, from 04_)
  + independent locus  (genomic coordinates, here)
  = copy number

A tree establishes what a sequence is. It cannot establish that two proteins
are two genomic copies rather than one gene split across two models, two
haplotypes of one locus, or a duplicated assembly region. Taro is heterozygous
at 1.83% in this cultivar, so allelic sequence retained separately is a live
possibility rather than a formality.

Criteria for an independent locus, from docs/METHODS_part1.md:

  - a different anchored chromosome from its putative paralog, OR
  - the same chromosome, non-overlapping coordinates, AND distinct flanking
    gene content
  - and not on an unanchored scaffold

Flanking gene content is read from the GFF3 and reported, not assumed.

Output is one machine-readable table. Part 2 modifies exactly these columns, so
this is the join between the parts.
"""
import csv
import os
import re
from collections import defaultdict
from pathlib import Path

REFS = Path(os.environ["TARO_REFS"])
WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
GFF = REFS / "proteomes" / "Colocasia_esculenta.Genome.V1.gff3"
FAM = WORK / "04_families"
TAB = CODE / "results" / "tables"

ALRT_MIN, BOOT_MIN = 80.0, 95
FLANK_N = 3          # genes either side used for the flanking comparison

# ------------------------------------------------------------------ inputs
cands = list(csv.DictReader(open(TAB / "homology_candidates.tsv"), delimiter="\t"))
sens = {r["gene"]: r for r in
        csv.DictReader(open(TAB / "homology_sensitivity.tsv"), delimiter="\t")}

# ------------------------------------------------------------------ GFF3
genes_by_scaf = defaultdict(list)
coord = {}
print("reading the taro annotation")
with open(GFF) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        m = re.search(r"ID=(?:gene-)?([^;]+)", f[8])
        if not m:
            continue
        gid, s, e, strand = m.group(1), int(f[3]), int(f[4]), f[6]
        genes_by_scaf[f[0]].append((s, e, gid, strand))
        coord[gid] = (f[0], s, e, strand)
for k in genes_by_scaf:
    genes_by_scaf[k].sort()
order = {sc: {g[2]: i for i, g in enumerate(v)} for sc, v in genes_by_scaf.items()}
print(f"  {len(coord)} genes on {len(genes_by_scaf)} sequences")


def flanking(gid, n=FLANK_N):
    if gid not in coord:
        return []
    sc = coord[gid][0]
    i = order[sc].get(gid)
    if i is None:
        return []
    lo, hi = max(0, i - n), min(len(genes_by_scaf[sc]), i + n + 1)
    return [g[2] for g in genes_by_scaf[sc][lo:hi] if g[2] != gid]


# ------------------------------------------------------------------ trees
def support(node):
    raw = node.name or (str(node.confidence) if node.confidence is not None else "")
    m = re.match(r"^([\d.]+)/(\d+)$", str(raw).strip())
    return (float(m.group(1)), int(m.group(2))) if m else (None, None)


family_members = {}   # gene -> {taro accession: (alrt, boot)}
for fa in sorted(FAM.glob("*_B.faa")):
    gene = fa.name[:-6]
    members = set()
    for ln in open(fa):
        if ln.startswith(">") and ln[1:].startswith("colesc_"):
            members.add(ln[1:].strip().split("_", 1)[1])
    sup = {}
    tf = FAM / f"{gene}_B_iq.treefile"
    if tf.exists():
        try:
            from Bio import Phylo
            t = Phylo.read(str(tf), "newick")
            t.root_at_midpoint()
            for tip in t.get_terminals():
                if not tip.name or not tip.name.startswith("colesc_"):
                    continue
                acc = tip.name.split("_", 1)[1]
                best = (None, None)
                for nd in reversed(t.get_path(tip)[:-1]):
                    a, b = support(nd)
                    if a is not None:
                        best = (a, b)
                        break
                sup[acc] = best
        except Exception as exc:
            print(f"  could not read {tf.name}: {exc}")
    family_members[gene] = {m: sup.get(m, (None, None)) for m in members}

print(f"  {len(family_members)} families with trees: "
      f"{', '.join(sorted(family_members)) or 'none'}")

# ------------------------------------------------------------------ rows
rows = []
for c in cands:
    g = c["taro_gene"]
    if not g:
        continue
    gene, kind = c["gene"], c["claim_kind"]
    sc, s, e, strand = coord.get(g, ("", "", "", ""))
    aa = ""
    assigned = ""
    alrt = boot = ""

    if gene in family_members:
        if g in family_members[gene]:
            assigned = gene
            a, b = family_members[gene][g]
            if a is not None:
                alrt, boot = f"{a:.1f}", str(b)
        else:
            assigned = "outside family"
    elif kind == "presence":
        assigned = "n/a (presence claim)"
    else:
        assigned = "RBH-only"

    rows.append(dict(
        family=gene, at_locus=c["at_locus"], taro_gene=g, claim_kind=kind,
        assignment=assigned, identity=c["identity"],
        coverage_pct=c["coverage_pct"], aln_aa=c["aln"], call=c["call"],
        scaffold=sc, start=s, end=e, strand=strand,
        tree_alrt=alrt, tree_boot=boot,
        candidate_set_stable=sens.get(gene, {}).get("stable", ""),
    ))

# ---------------------------------------------- independent locus + status
by_fam = defaultdict(list)
for r in rows:
    if r["assignment"] == r["family"]:
        by_fam[r["family"]].append(r)

print()
print("=" * 78)
print("INDEPENDENT-LOCUS TEST")
print("=" * 78)

for fam, members in sorted(by_fam.items()):
    print(f"\n{fam}  —  {len(members)} family member(s)")
    for r in members:
        anchored = str(r["scaffold"]).lower().startswith(("chr", "superscaffold"))
        r["anchored"] = "yes" if anchored else "no"
        print(f"  {r['taro_gene']:<12} {r['scaffold']:<18} "
              f"{r['start']}-{r['end']} {r['strand']}")

    if len(members) == 1:
        members[0]["independent_locus"] = "n/a (single member)"
        members[0]["status"] = "single copy"
        continue

    for i, r in enumerate(members):
        verdicts = []
        for j, o in enumerate(members):
            if i == j:
                continue
            if r["scaffold"] != o["scaffold"]:
                verdicts.append("different sequence")
                continue
            overlap = not (int(r["end"]) < int(o["start"])
                           or int(o["end"]) < int(r["start"]))
            if overlap:
                verdicts.append("OVERLAPPING")
                continue
            fr, fo = set(flanking(r["taro_gene"])), set(flanking(o["taro_gene"]))
            shared = fr & fo
            gap = (int(o["start"]) - int(r["end"]) if int(o["start"]) > int(r["end"])
                   else int(r["start"]) - int(o["end"]))
            if shared:
                verdicts.append(f"shares {len(shared)} flanking genes")
            else:
                verdicts.append(f"same sequence, {gap:,} bp apart, distinct flanks")

        bad = any("OVERLAP" in v or "shares" in v for v in verdicts)
        if r["anchored"] == "no":
            r["independent_locus"] = "unresolved (unanchored sequence)"
        elif bad:
            r["independent_locus"] = "no"
        else:
            r["independent_locus"] = "yes"
        print(f"    {r['taro_gene']:<12} -> {r['independent_locus']}"
              f"   [{'; '.join(verdicts)}]")

    n_ind = sum(1 for r in members if r["independent_locus"] == "yes")
    for k, r in enumerate(members, 1):
        r["status"] = (f"copy {k} of {n_ind}" if r["independent_locus"] == "yes"
                       else "not an independent locus")

# -------------------------------------------------------------- evidence
for r in rows:
    r.setdefault("independent_locus", "")
    r.setdefault("status", "")
    r.setdefault("anchored", "")
    if r["claim_kind"] == "presence":
        r["evidence"] = ("presence" if r["call"].startswith("ortholog")
                         else "presence (weak coverage)")
    elif r["claim_kind"] == "copy-number":
        tree_ok = (r["tree_boot"] != "" and float(r["tree_alrt"] or 0) >= ALRT_MIN
                   and int(r["tree_boot"] or 0) >= BOOT_MIN)
        if r["assignment"] != r["family"]:
            r["evidence"] = "not in family"
        elif r["independent_locus"] == "yes" and tree_ok:
            r["evidence"] = "copy-number verified"
        elif r["independent_locus"] == "yes":
            r["evidence"] = "locus verified, tree support below threshold"
        else:
            r["evidence"] = "family assigned, locus unresolved"
    else:
        r["evidence"] = "candidate (RBH)"

# ------------------------------------------------------------------ out
cols = ["family", "at_locus", "taro_gene", "claim_kind", "assignment",
        "identity", "coverage_pct", "aln_aa", "call", "scaffold", "start",
        "end", "strand", "anchored", "independent_locus", "tree_alrt",
        "tree_boot", "candidate_set_stable", "evidence", "status"]
out = TAB / "part1_inventory.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in cols})

print()
print("=" * 78)
print("SUMMARY BY EVIDENCE LEVEL")
print("=" * 78)
tot = defaultdict(int)
for r in rows:
    tot[r["evidence"]] += 1
for k in sorted(tot, key=lambda x: -tot[x]):
    print(f"  {tot[k]:>3}  {k}")

print()
print("  copy-number claims that pass both limbs:")
any_ = False
for r in rows:
    if r["evidence"] == "copy-number verified":
        print(f"    {r['family']:<9} {r['taro_gene']:<12} {r['status']:<14} "
              f"{r['scaffold']} {r['tree_alrt']}/{r['tree_boot']}")
        any_ = True
if not any_:
    print("    none")

print(f"""
  Counts in the 'candidate (RBH)' rows are lower bounds on family size, not
  copy numbers. Reciprocal best hit misses a paralog divergent enough that its
  best Arabidopsis hit is a different gene, so a detection count cannot be read
  as a family size.

  Any figure drawn from this table encodes the evidence column. A figure that
  displays every count identically is making a copy-number claim for every
  gene, whatever the caption says.

written: {out}
""")
