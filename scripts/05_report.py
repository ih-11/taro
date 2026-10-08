#!/usr/bin/env python3
"""
05_report.py — genomic coordinates, locus independence, and the master table.

FOUR THINGS WRONG WITH THE FIRST RUN. Three changed the counts.

1. "NOTHING ANNOTATED BETWEEN THEM" IS NOT EVIDENCE IN THIS GENOME.

   Ces13250 and Ces13251 are 75,073 bp apart on Superscaffold7, same strand,
   with no annotated gene between, and the first version called that one locus
   reported as two models. It is not. Taro carries 28,253 genes in 2,318 Mb, so
   the mean spacing between genes is about 82 kb. Two genes 75 kb apart with
   nothing between them is the single most ordinary arrangement in this
   annotation, not a signal.

   And the sequence evidence said so already. A split model's two pieces cover
   DIFFERENT parts of the anchor: Ces12496 and Ces12497 cover query residues
   126-238 and 281-427, zero overlap. Ces13250 and Ces13251 cover 69% and 88%
   of the same anchor, which cannot be two non-overlapping pieces of one gene.
   03_'s split detector correctly never flagged them; only this script's cruder
   rule did.

   So the call now needs both halves, which is what was intended all along:

     sequence : query ranges overlap by at most 25% of the shorter, AND run in
                the same order as the genomic positions given the strand
     genome   : same strand, no annotated gene between, and the combined span
                reported against the family's own gene spans

   Distance alone is not used as a criterion, because in a genome this sparse it
   carries almost no information.

2. EVERY EXON COUNT WAS ZERO.

   Exons name their parent mRNA, "rna-Ces00001.1", which strips to
   "Ces00001.1". Genes are keyed "Ces00001". So every count was filed under an
   identifier nothing looked up. The mRNA-to-gene map was built and then never
   used. Exon counts now fold through it.

3. SUPERSCAFFOLD12 WAS CALLED A SHORT SEQUENCE.

   The 90% cumulative-length rule cut between Superscaffold13 at 113 Mb and
   Superscaffold12 at 105 Mb, so a 105 Mb sequence carrying 1,963 genes was
   treated as unplaceable, and Ces27429 and Ces28200 were struck from the NCED
   count for sitting on it.

   The assembly declares its own answer and I invented a threshold instead of
   reading it. The 79 sequences carrying genes are 14 named Superscaffold1 to
   Superscaffold14 and 65 named unanchor*. Taro is 2n = 28, so 14 chromosomes,
   and the names line up exactly. Anchored is now "not named unanchor", with the
   14-against-n=14 check printed so the reasoning is visible rather than
   buried in a constant.

4. THE OVER-MERGED CHECK COMPARED A GENE TO A MEDIAN THAT INCLUDED ITS OWN
   FRAGMENTS.

   It flagged Ces24605, the full-length PSY, for spanning 7,233 bp against a
   "target median" of 2,365 bp. That median was taken over Ces24605, Ces12496
   and Ces12497, so the full gene was being compared against the two fragments
   of a different locus. Circular. The length flag now comes from protein
   length against the family's reference proteins, which is what 04b_ already
   computes honestly, and genomic span and exon count are reported as context
   rather than as the trigger.

   A length-flagged model also no longer reports as a verified copy number.
   Ces14428 came out "copy number, verified: 1" while simultaneously being
   under suspicion of spanning several genes. Those cannot both be printed.

   And a diagnostic is added that settles it: the flagged protein is searched
   against Arabidopsis keeping all HSPs, and the subject gene of each query
   region is reported. One gene across the whole protein means a genuinely long
   protein. Different genes in different regions means a merge.

THE GOVERNING RULE, which this script completes

  Copy number = family assignment (phylogeny) + independent genomic locus
                (coordinates). Presence = reciprocal homology with adequate
                coverage. A detection count is not a copy number.
"""
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REFS = Path(os.environ["TARO_REFS"])
WORK = Path(os.environ["TARO_WORK"])
CODE = Path(os.environ["TARO_CODE"])
TBL = CODE / "results" / "tables"
D = REFS / "proteomes" / "orthofinder_input"

SPLIT_MAX_QOVERLAP = 0.25     # same rule as 03_
LENGTH_FLAG = 2.0             # x the family's reference median
FLANK = 3
TARO_CHROMS = 14              # 2n = 28

TBL.mkdir(parents=True, exist_ok=True)


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def read_fasta(p):
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].strip(); d[k] = []
        elif k is not None:
            d[k].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


search = []
for d in (REFS / "proteomes", REFS / "annotation", REFS / "assembly"):
    if d.exists():
        for pat in ("*.gff3", "*.gff", "*.gff3.gz", "*.gff.gz"):
            search.extend(sorted(d.glob(pat)))
GFF = next((p for p in search if "olocasia" in p.name or "colesc" in p.name.lower()),
           search[0] if search else None)
if GFF is None:
    sys.exit(f"no taro GFF3 found under {REFS/'proteomes'}")

print("=" * 108)
print("GFF3")
print("=" * 108)
print(f"  {GFF}")

id_re = re.compile(r"(?:^|;)ID=([^;]+)")
par_re = re.compile(r"(?:^|;)Parent=([^;]+)")


def clean(x):
    for pre in ("gene-", "gene:", "rna-", "mRNA:", "transcript:"):
        if x.startswith(pre):
            return x[len(pre):]
    return x


genes, by_seq = {}, defaultdict(list)
mrna2gene, exon_by_mrna = {}, defaultdict(int)
feat = defaultdict(int)
with opener(GFF) as fh:
    for ln in fh:
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        seq, typ, start, end, strand, attr = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        feat[typ] += 1
        if typ in ("gene", "pseudogene"):
            m = id_re.search(attr)
            if not m:
                continue
            g = clean(m.group(1))
            genes[g] = dict(gene=g, seq=seq, start=start, end=end, strand=strand)
            by_seq[seq].append(g)
        elif typ in ("mRNA", "transcript"):
            m, p = id_re.search(attr), par_re.search(attr)
            if m and p:
                mrna2gene[clean(m.group(1))] = clean(p.group(1))
        elif typ == "exon":
            p = par_re.search(attr)
            if p:
                exon_by_mrna[clean(p.group(1))] += 1

# fold exon counts through the mRNA-to-gene map: fix 2
exons = defaultdict(int)
for mrna, n in exon_by_mrna.items():
    exons[mrna2gene.get(mrna, mrna)] = max(exons[mrna2gene.get(mrna, mrna)], n)

for s in by_seq:
    by_seq[s].sort(key=lambda g: genes[g]["start"])
order = {g: i for s in by_seq for i, g in enumerate(by_seq[s])}
seqlen = {s: max(genes[g]["end"] for g in by_seq[s]) for s in by_seq}

print(f"  features: " + ", ".join(f"{k}={v}" for k, v in
                                  sorted(feat.items(), key=lambda kv: -kv[1])[:6]))
print(f"  genes {len(genes)}   sequences {len(by_seq)}   "
      f"mRNA->gene {len(mrna2gene)}   genes with exon counts "
      f"{sum(1 for g in genes if exons.get(g))}")
tot = sum(seqlen.values())
print(f"  assembly carrying genes: {tot/1e6:,.0f} Mb, "
      f"mean spacing {tot/len(genes)/1e3:,.0f} kb per gene")

# ------------------------------------------------- anchored by name: fix 3
print()
print("=" * 108)
print("ANCHORED SEQUENCES  —  taken from the assembly's own naming")
print("=" * 108)
anchored = {s for s in by_seq if not s.lower().startswith("unanchor")}
unanch = set(by_seq) - anchored
print(f"\n  {len(anchored)} sequences not named 'unanchor', "
      f"{len(unanch)} named 'unanchor'")
print(f"  taro is 2n = 28, so {TARO_CHROMS} chromosomes expected: "
      + ("the counts agree" if len(anchored) == TARO_CHROMS
         else f"MISMATCH, {len(anchored)} anchored against {TARO_CHROMS} expected"))
print()
for s in sorted(anchored, key=lambda x: -seqlen[x]):
    print(f"    anchored  {s:<22} {seqlen[s]/1e6:>7.1f} Mb  "
          f"{len(by_seq[s]):>5} genes")
print(f"    unanchored: {len(unanch)} sequences, "
      f"{sum(seqlen[s] for s in unanch)/1e6:,.0f} Mb, "
      f"{sum(len(by_seq[s]) for s in unanch):,} genes")
print()
print("  A gene on an unanchored sequence is reported but does not contribute a")
print("  copy, because in a 1.83% heterozygous assembly a haplotig of a gene")
print("  already counted cannot be told from a paralog by coordinates alone.")

# ---------------------------------------------------------------- candidates
hits = list(csv.DictReader(open(TBL / "homology_hits_all.tsv"), delimiter="\t"))
anchors = {r["target"]: r for r in
           csv.DictReader(open(TBL / "pathway_anchors.tsv"), delimiter="\t")}
assign = {}
if (TBL / "family_assignment.tsv").exists():
    for r in csv.DictReader(open(TBL / "family_assignment.tsv"), delimiter="\t"):
        assign[r["taro_gene"]] = r

MIN_ALN, MIN_COV = 100, 50.0
best = {}
for h in hits:
    if h["rbh_family"] != "1":
        continue
    k = (h["target"], h["taro_gene"])
    if k not in best or float(h["identity"]) > float(best[k]["identity"]):
        best[k] = h
want = defaultdict(set)
for (tgt, g) in best:
    want[tgt].add(g)
allwant = {g for _, g in best}

print()
print("=" * 108)
print("LOCATING THE CANDIDATES")
print("=" * 108)
missing = sorted(g for g in allwant if g not in genes)
print(f"\n  to locate {len(allwant)}, located {len(allwant)-len(missing)}")
if missing:
    print(f"  NOT LOCATED: {' '.join(missing[:20])}")
    sys.exit(2)

# --------------------------------------------- one locus or two: fix 1
print()
print("=" * 108)
print("ONE LOCUS OR TWO")
print("=" * 108)
print("\n  A split model needs sequence evidence AND genome evidence. Genomic")
print(f"  distance is not a criterion: mean gene spacing here is "
      f"{tot/len(genes)/1e3:,.0f} kb, so two")
print("  genes tens of kb apart with nothing between is unremarkable.\n")

same_locus, notes = {}, defaultdict(list)
for tgt, gl in sorted(want.items()):
    gl = sorted(gl)
    fam_spans = [genes[g]["end"] - genes[g]["start"] + 1 for g in gl
                 if int(best[(tgt, g)]["aln_aa"]) >= MIN_ALN
                 and float(best[(tgt, g)]["qcov_pct"]) >= MIN_COV]
    for i in range(len(gl)):
        for j in range(i + 1, len(gl)):
            a, b = gl[i], gl[j]
            A, B = genes[a], genes[b]
            if A["seq"] != B["seq"]:
                continue
            ha, hb = best[(tgt, a)], best[(tgt, b)]
            qa = (int(ha["q_start"]), int(ha["q_end"]))
            qb = (int(hb["q_start"]), int(hb["q_end"]))
            la, lb = qa[1] - qa[0] + 1, qb[1] - qb[0] + 1
            qov = max(0, min(qa[1], qb[1]) - max(qa[0], qb[0]) + 1)
            qov_frac = qov / min(la, lb) if min(la, lb) else 1.0
            gov = max(0, min(A["end"], B["end"]) - max(A["start"], B["start"]) + 1)
            i0, j0 = sorted((order[a], order[b]))
            btw = max(0, j0 - i0 - 1)
            same_strand = A["strand"] == B["strand"]
            # collinear: genomic order matches query order, given the strand
            first_genomic = a if A["start"] < B["start"] else b
            first_query = a if qa[0] < qb[0] else b
            collinear = ((first_genomic == first_query) if A["strand"] == "+"
                         else (first_genomic != first_query))
            comb = max(A["end"], B["end"]) - min(A["start"], B["start"]) + 1
            fam_max = max(fam_spans) if fam_spans else None

            print(f"\n  {tgt}:  {a} and {b} on {A['seq']}")
            print(f"    {a}: {A['start']:>11,}-{A['end']:<11,} {A['strand']}  "
                  f"{exons.get(a,0):>2} exons   anchor residues {qa[0]}-{qa[1]} "
                  f"(qcov {ha['qcov_pct']}%)")
            print(f"    {b}: {B['start']:>11,}-{B['end']:<11,} {B['strand']}  "
                  f"{exons.get(b,0):>2} exons   anchor residues {qb[0]}-{qb[1]} "
                  f"(qcov {hb['qcov_pct']}%)")
            print(f"    query overlap {qov} aa ({qov_frac*100:.0f}% of the shorter)"
                  f"   genomic overlap {gov} bp   {btw} genes between   "
                  f"{'same' if same_strand else 'opposite'} strand")

            if gov > 0:
                print("    -> genomic coordinates overlap: not two loci")
                same_locus[b] = a
                notes[b].append(f"overlaps {a}")
                continue
            seq_ok = qov_frac <= SPLIT_MAX_QOVERLAP
            if not seq_ok:
                print(f"    -> the two cover the SAME part of the anchor "
                      f"({qov_frac*100:.0f}% overlap), so they cannot be two")
                print(f"       pieces of one gene. Two independent loci.")
                continue
            if not same_strand:
                print("    -> opposite strands: two independent loci")
                continue
            if btw > 0:
                print(f"    -> {btw} annotated gene(s) between: two independent loci")
                continue
            if not collinear:
                print("    -> query order does not match genomic order on this")
                print("       strand, so they are not collinear pieces of one gene")
                continue
            print(f"    combined span {comb:,} bp"
                  + (f", against {fam_max:,} bp for the largest complete member "
                     f"of this family" if fam_max else ""))
            print("    -> non-overlapping collinear anchor coverage, same strand,")
            print("       nothing annotated between: ONE LOCUS, two models")
            same_locus[b] = a
            notes[b].append(f"same locus as {a}, split model")
            notes[a].append(f"same locus as {b}, split model")

# --------------------------------------------- length anomalies: fix 4
print()
print("=" * 108)
print("GENE MODEL LENGTH")
print("=" * 108)
taro_fa = read_fasta(D / "colesc.fa")
plen = {k.split("|")[-1]: len(v) for k, v in taro_fa.items()}
flagged = {}
for tgt, gl in sorted(want.items()):
    # The family's reference median comes from tree B's input, which holds the
    # other sixteen species' orthologs and no taro. Comparing a gene against a
    # median that included its own fragments is what flagged the full-length PSY
    # last run.
    med = None
    fa = WORK / "04_families" / f"{anchors[tgt]['family_group'] or tgt}_B.faa"
    if fa.exists():
        L = sorted(len(v) for k, v in read_fasta(fa).items()
                   if not k.startswith("colesc_"))
        if L:
            med = L[len(L) // 2]
    for g in sorted(gl):
        n = plen.get(g)
        if med and n and n > LENGTH_FLAG * med:
            print(f"\n  {tgt}: {g} protein {n} aa against a family reference "
                  f"median of {med} aa ({n/med:.1f}x)")
            gi = genes[g]
            print(f"    {gi['seq']} {gi['start']:,}-{gi['end']:,} {gi['strand']}, "
                  f"span {gi['end']-gi['start']+1:,} bp, {exons.get(g,0)} exons")
            i = order[g]
            nb = by_seq[gi["seq"]][max(0, i-FLANK):i+FLANK+1]
            print(f"    neighbours: "
                  + " ".join(f"{x}{'*' if x == g else ''}" for x in nb))
            flagged[g] = (tgt, n, med)

if flagged:
    print()
    print("  DOMAIN DIAGNOSTIC  —  one domain, or several genes in one record?")
    print("  " + "-" * 104)
    print("  The first version of this test counted DISTINCT Arabidopsis genes")
    print("  among the HSPs and called more than one a merge. That is wrong: for")
    print("  Ces14428 all fourteen matches landed on residues 1907-2232, the same")
    print("  stretch of protein, because fourteen members of one superfamily all")
    print("  match one prenyltransferase domain. Fourteen genes hitting ONE region")
    print("  is a domain. Different genes hitting DIFFERENT regions is a merge.")
    print("  So HSPs are merged into query regions first, and what matters is how")
    print("  much of the protein matches nothing at all.")
    db = WORK / "03_homology" / "aratha"
    mp = WORK / "05_orthology" / "aratha_protein2gene.tsv"
    sym = {}
    for ln in open(mp):
        if ln.startswith("protein\t"):
            continue
        pr, gg, sy, loc = ln.rstrip("\n").split("\t")
        sym[pr] = sy or loc
    for g in flagged:
        key = next(k for k in taro_fa if k.endswith("|" + g))
        L = len(taro_fa[key])
        q = WORK / "04_families" / f"{g}_probe.faa"
        q.write_text(f">{g}\n{taro_fa[key]}\n")
        out = WORK / "04_families" / f"{g}_probe.tsv"
        try:
            subprocess.run(["diamond", "blastp", "-q", q, "-d", db, "-o", out,
                            "--outfmt", "6", "qseqid", "sseqid", "pident",
                            "length", "qstart", "qend", "evalue", "bitscore",
                            "--max-hsps", "20", "--max-target-seqs", "60",
                            "--evalue", "1e-5", "--quiet"], check=True)
        except Exception as e:
            print(f"    {g}: probe failed ({e})")
            continue
        hs = []
        for ln in open(out):
            f = ln.rstrip("\n").split("\t")
            hs.append((int(f[4]), int(f[5]), sym.get(f[1].split("|")[-1], f[1]),
                       float(f[7])))
        if not hs:
            print(f"\n    {g} ({L} aa): no Arabidopsis match at e < 1e-5 anywhere")
            continue
        # merge HSPs into query regions regardless of which gene they came from
        hs.sort()
        regions = []
        for qs, qe, sy, bits in hs:
            if regions and qs <= regions[-1]["end"]:
                r = regions[-1]
                r["end"] = max(r["end"], qe)
                r["hits"].append((sy, bits))
            else:
                regions.append(dict(start=qs, end=qe, hits=[(sy, bits)]))
        covered = sum(r["end"] - r["start"] + 1 for r in regions)
        print(f"\n    {g} ({L} aa): {len(regions)} matched region(s), "
              f"{covered} aa matched, {L - covered} aa ({(L-covered)/L*100:.0f}%) "
              f"matching nothing")
        for r in regions:
            top = sorted(set(r["hits"]), key=lambda t: -t[1])[:4]
            n = len({h[0] for h in r["hits"]})
            print(f"      {r['start']:>5}-{r['end']:<5} "
                  f"({r['end']-r['start']+1:>4} aa)  {n:>2} Arabidopsis genes, "
                  f"best: " + ", ".join(t[0] for t in top))
        if len(regions) == 1:
            r = regions[0]
            print(f"    -> ONE matched region. The {n if False else len({h[0] for h in r['hits']})} "
                  f"Arabidopsis genes matching it are one")
            print(f"       superfamily matching one domain, not several genes merged.")
            print(f"       What needs explaining is the {L - covered} aa "
                  f"({(L-covered)/L*100:.0f}%) that matches nothing:")
            print(f"       either a lineage-specific extension or an annotation that")
            print(f"       ran one gene into adjacent sequence. Either way this model")
            print(f"       is not one clean copy of anything, and the copy number is")
            print(f"       withheld rather than guessed.")
        else:
            disjoint = all(not ({h[0] for h in regions[i]["hits"]}
                                & {h[0] for h in regions[j]["hits"]})
                           for i in range(len(regions))
                           for j in range(i + 1, len(regions)))
            print(f"    -> {len(regions)} separate regions"
                  + (", matching DIFFERENT Arabidopsis genes: consistent with"
                     " several genes merged into one record" if disjoint
                     else ", sharing matches, so more likely repeated domains"
                          " than merged genes"))
else:
    print("\n  no candidate protein exceeds twice its family's reference median")

# ---------------------------------------------------------------- master table
rows = []
for (tgt, g), h in sorted(best.items()):
    a, gi = anchors[tgt], genes[g]
    counted = int(int(h["aln_aa"]) >= MIN_ALN and float(h["qcov_pct"]) >= MIN_COV)
    asg = assign.get(g)
    rep = same_locus.get(g)
    anch = gi["seq"] in anchored
    if rep:
        indep, locus = 0, f"shared with {rep}"
    elif not anch:
        indep, locus = 0, "unanchored sequence"
    else:
        indep, locus = 1, "independent"

    if a["claim_kind"] == "presence":
        ev = "present" if counted else "present, below the counting floor"
    elif a["claim_kind"] == "detection":
        ev = "candidate (RBH)" if counted else "sub-threshold hit, recorded not counted"
    elif g in flagged:
        ev = "family assigned; gene model length anomalous, copy number not called"
    elif asg is None:
        ev = "no tree; locus evidence only"
    elif asg["assignment"] == "NOT DETERMINABLE":
        ev = "family assignment not determinable"
    elif not indep:
        ev = f"family assigned; {locus}"
    elif str(asg.get("treeA_support_ok")) == "1":
        ev = "copy-number verified"
    else:
        ev = "locus verified, family support below threshold"

    rows.append(dict(
        pathway=a["pathway"], target=tgt, claim_kind=a["claim_kind"],
        anchor_locus=a["anchor_locus"], family_group=a["family_group"],
        taro_gene=g, protein_aa=plen.get(g, ""), identity=h["identity"],
        qcov_pct=h["qcov_pct"], aln_aa=h["aln_aa"],
        anchor_from=h["q_start"], anchor_to=h["q_end"], counted=counted,
        rbh_anchor=h["rbh_anchor"], rbh_family=h["rbh_family"],
        assignment=(asg or {}).get("assignment", ""),
        member_alrt=(asg or {}).get("member_alrt", ""),
        member_boot=(asg or {}).get("member_boot", ""),
        tree_alrt=(asg or {}).get("treeA_alrt", ""),
        tree_boot=(asg or {}).get("treeA_boot", ""),
        treeB_clade=(asg or {}).get("treeB_split", ""),
        seq=gi["seq"], start=gi["start"], end=gi["end"], strand=gi["strand"],
        span_bp=gi["end"] - gi["start"] + 1, exons=exons.get(g, 0),
        anchored=int(anch), independent_locus=indep, locus_note=locus,
        length_flag=int(g in flagged), notes="; ".join(notes.get(g, [])),
        evidence=ev))

out = TBL / "part1_inventory.tsv"
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

print()
print("=" * 108)
print("THE CALLS")
print("=" * 108)
print("\n  COPY-NUMBER TARGETS")
print("  " + "-" * 104)
for tgt, a in anchors.items():
    if a["claim_kind"] != "copy-number":
        continue
    rr = [r for r in rows if r["target"] == tgt]
    ver = [r for r in rr if r["evidence"] == "copy-number verified"]
    low = [r for r in rr if r["evidence"].startswith("locus verified")]
    oth = [r for r in rr if r not in ver and r not in low]
    print(f"\n  {tgt}")
    print(f"    copy number, verified         : {len(ver)}"
          + (f"   {' '.join(r['taro_gene'] for r in ver)}" if ver else ""))
    if low:
        print(f"    locus verified, support below : {len(low)}"
              f"   {' '.join(r['taro_gene'] for r in low)}")
    for r in oth:
        print(f"    {r['taro_gene']:<12} {r['evidence']}")

print("\n  DETECTION TARGETS, lower bounds")
print("  " + "-" * 104)
for tgt, a in anchors.items():
    if a["claim_kind"] != "detection":
        continue
    c = [r for r in rows if r["target"] == tgt and r["counted"]]
    s = [r for r in rows if r["target"] == tgt and not r["counted"]]
    loci = len({r["taro_gene"] for r in c if r["independent_locus"]})
    print(f"    {tgt:<16} {len(c):>2} counted, {loci:>2} independent loci"
          + (f", {len(s)} below the floor" if s else ""))

print("\n  PRESENCE TARGETS")
print("  " + "-" * 104)
for tgt, a in anchors.items():
    if a["claim_kind"] != "presence":
        continue
    rr = [r for r in rows if r["target"] == tgt]
    print(f"    {tgt:<16} {'present' if rr else 'ABSENT'}"
          + (f"   {' '.join(r['taro_gene'] for r in rr)}" if rr else ""))

print(f"""
{'=' * 108}

written: {out}   ({len(rows)} rows)

next:  update README.md and docs/, then rebuild
       notebooks/part1_landscape.ipynb against this table
""")
