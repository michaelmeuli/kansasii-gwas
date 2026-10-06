"""Region FCNCLM_01713..01736 (NC_022663.1:1,944,065-1,970,050). Count variants vs the reference in this 26 kb region per kansasii genome,
classify divergent-haplotype (>50 variants) vs reference-like, and test against amikacin in MGIT (35 kansasii) and in the independent MIC cohort."""
from __future__ import annotations

import os
from typing import cast
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); G = ROOT / "output/gwas"; RES = ROOT / "runs/mkan329/assembly/results"
LO, HI = 1944065, 1970050
lin = pd.read_csv(G / "lineage/kansasii_clades.csv", index_col=0)
names = list(lin.index)
V: dict[str, set[tuple[int, str, str]]] = {}
for s in names:
    t = pd.read_csv(RES / s / "5_typing/kansasii_snippy/snippy_out/snps.tab", sep="\t", dtype=str, usecols=["CHROM", "POS", "REF", "ALT"])
    t["POS"] = t.POS.astype(int); t = t[(t.CHROM == "NC_022663.1") & (t.POS >= LO) & (t.POS <= HI)]
    V[s] = set(zip(t.POS, t.REF, t.ALT))
lin["n_var_region"] = [len(V[s]) for s in names]
lin["haplotype"] = np.where(lin.n_var_region > 50, "divergent", "ref-like")
print("distribution of region variant counts (all 82 kansasii):", np.percentile(lin.n_var_region, [0, 10, 25, 50, 75, 90, 100]).astype(int).tolist())
print("values between 1 and 800 (intermediate states):", lin.n_var_region[(lin.n_var_region > 0) & (lin.n_var_region < 800)].tolist())
print(pd.crosstab(lin.clade, lin.haplotype).to_string())
lin.to_csv(G / "lineage/kansasii_region_haplotype.csv")
div = [s for s in names if lin.haplotype[s] == "divergent"]
if len(div) > 1:
    j = [len(V[a] & V[b]) / len(V[a] | V[b]) for i, a in enumerate(div) for b in div[i + 1:]]
    print(f"pairwise Jaccard among divergent-haplotype isolates: mean {np.mean(j):.2f}, min {np.min(j):.2f}")

mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = lin[lin.NR.isin(ak.index)].copy(); x["R"] = ak.reindex(x.NR).values == "R"
t = pd.crosstab(x.haplotype, x.R); print("\nMGIT amikacin 1 mg/L, all kansasii (R = True):\n", t.to_string())
tab = np.array([[int(cast(int, t.loc["divergent", True])), int(cast(int, t.loc["divergent", False]))], [int(cast(int, t.loc["ref-like", True])), int(cast(int, t.loc["ref-like", False]))]])
print(f"  Fisher p = {fisher_exact(tab)[1]:.4f}")
for c in ("A", "B"):
    y = x[x.clade == c]; tt = pd.crosstab(y.haplotype, y.R); print(f"  clade {c}:\n{tt.to_string()}")

ph = pd.read_csv(G / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded & (ph.antibiotic == "Amikacin")].drop_duplicates("NR")
ph["haplotype"] = ph.NR.map(lin.set_index("NR").haplotype); ph = ph.dropna(subset=["haplotype"])
d, r = ph[ph.haplotype == "divergent"], ph[ph.haplotype == "ref-like"]
print(f"\nMIC cohort (independent, kansasii): amikacin log2 MIC divergent n={len(d)} median {d.log2_mic.median()} | ref-like n={len(r)} median {r.log2_mic.median()}; Mann-Whitney p = {mannwhitneyu(d.log2_mic, r.log2_mic).pvalue:.3f}")
print(ph.groupby("haplotype").log2_mic.value_counts().unstack(fill_value=0).to_string())
# same test for the other drugs, for context
ph2 = pd.read_csv(G / "geno_pheno_table_long.csv"); ph2 = ph2[ph2.has_ngs & ~ph2.excluded]
ph2["haplotype"] = ph2.NR.map(lin.set_index("NR").haplotype); ph2 = ph2.dropna(subset=["haplotype"])
rows = []
for drug, g in ph2.groupby("antibiotic"):
    g = g.drop_duplicates("NR"); a, b = g[g.haplotype == "divergent"].log2_mic, g[g.haplotype == "ref-like"].log2_mic
    rows.append((drug, len(a), len(b), a.median(), b.median(), mannwhitneyu(a, b).pvalue))
print("\nall drugs, MIC cohort kansasii, region haplotype:\n", pd.DataFrame(rows, columns=["drug", "n_div", "n_ref", "med_div", "med_ref", "p_unadj"]).sort_values("p_unadj").round(3).to_string(index=False))
