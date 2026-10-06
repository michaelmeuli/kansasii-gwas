"""Association controlled for the clonal frame: among isolates with an amikacin 1 mg/L MGIT result, find nearest relatives in the Gubbins
recombination-filtered tree and ask whether pairs that differ in region haplotype also differ in phenotype."""
from __future__ import annotations

from typing import Any
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Phylo
from scipy.stats import binomtest

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"
T = Phylo.read(G / "gubbins/kansasii_gubbins.final_tree.tre", "newick")  # type: ignore[attr-defined,no-untyped-call]
h = pd.read_csv(G / "lineage/kansasii_region_haplotype.csv", index_col=0)
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = h[h.NR.isin(ak.index)].copy(); x["ak"] = ak.reindex(x.NR).values; x["R"] = x.ak == "R"; x["H880"] = x.hap == "H~880"
ids = list(x.index); D = pd.DataFrame([[T.distance(a, b) for b in ids] for a in ids], index=ids, columns=ids)
print(f"{len(ids)} phenotyped kansasii: R {int(x.R.sum())}, S/I {int((~x.R).sum())}; haplotype-880 {int(x.H880.sum())}")
print("\nphylogenetic signal of the phenotype (filtered-tree distances, SNPs):")
R_ = [i for i in ids if x.R[i]]; N_ = [i for i in ids if not x.R[i]]
def w(a: list[str], b: list[str]) -> float:
    return float(np.mean([D.loc[i, j] for i in a for j in b if i != j]))
print(f"  R-R {w(R_, R_):.0f} | S/I-S/I {w(N_, N_):.0f} | R-S/I {w(R_, N_):.0f}")
print("\nnearest phenotyped relative of each isolate:")
rows: list[dict[str, Any]] = []
for i in ids:
    j = D.loc[i].drop(i).idxmin(); rows.append(dict(isolate=i, hap=str(x.hap[i]), ak=x.ak[i], nearest=j, dist=round(D.loc[i, j]), nearest_hap=str(x.hap[j]), nearest_ak=x.ak[j]))
nn = pd.DataFrame(rows); print(nn.to_string(index=False))
pairs = {tuple(sorted((i_, j_))) for i_, j_ in zip(nn.isolate.tolist(), nn.nearest.tolist())}
disc = [(a, b) for a, b in pairs if x.H880[a] != x.H880[b]]
print(f"\nunique nearest-neighbour pairs: {len(pairs)}; pairs discordant for haplotype-880: {len(disc)}")
conc = 0
for a, b in disc:
    c, n = (a, b) if x.H880[a] else (b, a)
    ok = bool(x.R[c]) and not x.R[n]; conc += ok
    print(f"  carrier {c} ({x.ak[c]}) vs non-carrier {n} ({x.ak[n]}), tree distance {D.loc[a, b]:.0f}: {'carrier R, non-carrier S/I (supports)' if ok else 'does not support'}")
if disc: print(f"  supporting {conc}/{len(disc)}; sign test (H0: 50%) p = {binomtest(conc, len(disc), 0.5).pvalue:.3f}")
x.to_csv(G / "gubbins/amikacin_pairs_table.csv")
