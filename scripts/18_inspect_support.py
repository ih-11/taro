#!/usr/bin/env python3
"""
Which nodes are well supported, and do the claims depend on them?

IQ-TREE writes internal labels as "SH-aLRT/UFBoot". The conventional
threshold is SH-aLRT >= 80 and UFBoot >= 95. Reporting an overall count is
not enough: what matters is whether the specific nodes a claim rests on are
among the supported ones.
"""
import os, re
from pathlib import Path
from Bio import Phylo

W = Path(os.environ["TARO_WORK"]) / "05_orthology" / "trees"

SP = {"colesc": "taro", "manesc": "cassava", "aratha": "Arabidopsis",
      "orysat": "rice", "zeamay": "maize", "musacu": "banana",
      "soltub": "potato", "nelnuc": "lotus", "ambtri": "Amborella",
      "zosmar": "Zostera"}

def label(l):
    if not l.name: return "?"
    tag, acc = l.name.split("_", 1)
    return f"{SP.get(tag, tag)}:{acc}"

for fam, focus in [("psy", "colesc"), ("ccd", "colesc")]:
    f = W / f"{fam}_iq.treefile"
    if not f.exists():
        print(f"{fam}: missing"); continue

    t = Phylo.read(str(f), "newick")
    print("\n" + "=" * 68)
    print(f"{fam.upper()}  —  path from each taro gene to the root")
    print("=" * 68)

    taro = [l for l in t.get_terminals() if l.name and l.name.startswith(focus)]
    for tip in taro:
        print(f"\n{label(tip)}")
        path = t.get_path(tip)
        for node in reversed(path[:-1]):
            tips = node.get_terminals()
            conf = node.confidence
            # Biopython stores "80/95" style labels in .confidence or .name
            raw = node.name if node.name else (str(conf) if conf else "")
            m = re.match(r"([\d.]+)/(\d+)", str(raw))
            if m:
                alrt, boot = float(m.group(1)), int(m.group(2))
                ok = "OK  " if (alrt >= 80 and boot >= 95) else "WEAK"
                sup = f"{alrt:5.1f}/{boot:3d}"
            else:
                ok, sup = "  - ", "    -   "
            sp = sorted({l.name.split('_')[0] for l in tips})
            print(f"   {ok} {sup}  {len(tips):3d} tips  "
                  f"{' '.join(SP.get(s, s) for s in sp)}")
