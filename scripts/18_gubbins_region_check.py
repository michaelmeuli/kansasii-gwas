"""Final Gubbins run: is the 26 kb region flagged, on which isolates, and where do the carriers sit in the recombination-filtered tree?
NOTE: iteration 4 (interim) showed a 20-isolate monophyletic clade with one shared block; the FINAL run does not. This script uses final files only."""
from __future__ import annotations

from typing import Any
import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Phylo

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"; GB = G / "gubbins"
LO, HI = 1944065, 1970050; L = HI - LO + 1
def grp(pat: str, text: str) -> str:
    m = re.search(pat, text)
    assert m is not None, pat
    return m.group(1)
rows: list[dict[str, Any]] = []
for l in open(GB / "kansasii_gubbins.recombination_predictions.gff"):
    if l.startswith("#") or not l.strip(): continue
    f = l.rstrip("\n").split("\t")
    rows.append(dict(start=int(f[3]), end=int(f[4]), taxa=grp(r'taxa="([^"]*)"', f[8]).split(),
                     snps=int(grp(r'snp_count="(\d+)"', f[8])), node=grp(r'node="([^"]*)"', f[8])))
d = pd.DataFrame(rows); d["len"] = d.end - d.start + 1
d["ov"] = [max(0, min(e, HI) - max(s, LO) + 1) for s, e in zip(d.start, d.end)]
print(f"{len(d)} recombination blocks genome-wide (median {int(d.len.median())} bp; {int(d.len.sum()):,} bp summed over branches)")
o = d[d.ov > 0].copy(); o["n_taxa"] = o.taxa.map(len)
print(f"blocks overlapping the region: {len(o)}; on branches shared by >=20 isolates: {(o.n_taxa >= 20).sum()}; single-isolate blocks: {(o.n_taxa == 1).sum()}")
h = pd.read_csv(G / "lineage/kansasii_region_haplotype.csv", index_col=0); H880 = set(h.index[h.hap == "H~880"])
def covfrac(t: str) -> float:
    iv = sorted((max(s, LO), min(e, HI)) for s, e, tx in zip(o.start, o.end, o.taxa) if t in tx); tot = 0; cur: list[int] | None = None
    for a, b in iv:
        if cur is None or a > cur[1] + 1:
            if cur: tot += cur[1] - cur[0] + 1
            cur = [a, b]
        else: cur[1] = max(cur[1], b)
    return (tot + (cur[1] - cur[0] + 1 if cur else 0)) / L
h["covered"] = [covfrac(t) for t in h.index]
print("fraction of the region covered by blocks, by haplotype class:\n", h.groupby(h.hap.astype(str)).covered.agg(["count", "mean", "min", "max"]).round(2).to_string())
T = Phylo.read(GB / "kansasii_gubbins.final_tree.tre", "newick")  # type: ignore[attr-defined,no-untyped-call]
m = T.common_ancestor(sorted(H880)); desc = {t.name for t in m.get_terminals()}
print(f"\nfinal tree: MRCA of the {len(H880)} haplotype-880 isolates spans {len(desc)} tips ({len(desc - H880)} non-carriers inside)")
Hl = sorted(H880); iu = np.triu_indices(len(Hl), 1)
dc = np.array([[T.distance(a, b) for b in Hl] for a in Hl])[iu]
rest = [t for t in h.index if t not in H880]; rng = np.random.default_rng(1)
dn = [T.distance(*rng.choice(rest, 2, replace=False)) for _ in range(500)]
print(f"filtered-tree pairwise distance (SNPs): carrier pairs median {np.median(dc):.0f} (min {dc.min():.0f}, max {dc.max():.0f}); non-carrier pairs median {np.median(dn):.0f}")
h.to_csv(GB / "region_block_coverage_by_isolate.csv")
