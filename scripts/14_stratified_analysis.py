"""Lineage-stratified re-analysis. Clades A (n=20) and B (n=60) = the two main groups of the 4-way average-linkage cut of core-SNP
distances (script 13); the remaining 2 isolates are single long-branch genomes. (1) amikacin 1 mg/L MGIT by clade;
(2) gene-level Fisher tests within clade B only; (3) the same within-clade question for ANI; (4) MIC cohort: clade vs MIC, all drugs."""
from __future__ import annotations

import os
from typing import Any
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import fcluster
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); G = ROOT / "output/gwas"; O = G / "lineage"
D = pd.read_csv(O / "core_snp_distance_matrix.csv", index_col=0); names = D.index.tolist()
cl4 = pd.Series(fcluster(np.load(O / "linkage.npy"), 4, "maxclust"), index=names)
size = cl4.value_counts()
clade = cl4.map(lambda c: "B" if c == size.index[0] else "A" if c == size.index[1] else "single").rename("clade")
print("clade sizes:", clade.value_counts().to_dict())
meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "gtdb_ani"]).set_index("PROBENNUMMER")
nr = meta.NR; lin = pd.DataFrame({"clade": clade, "NR": nr.reindex(names).values, "ani": meta.gtdb_ani.reindex(names).values}, index=names)
lin.to_csv(O / "kansasii_clades.csv")

mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = lin[lin.NR.isin(ak.index)].copy(); x["ak"] = ak.reindex(x.NR).values; x["R"] = x.ak == "R"
print("\n(1) amikacin 1 mg/L by clade (R / S+I):\n", x.groupby("clade").R.agg(R="sum", n="size").assign(nonR=lambda d: d.n - d.R).to_string())
a, b = x[x.clade == "A"], x[x.clade == "B"]
print(f"   clade A vs B (R vs S+I): Fisher p = {fisher_exact(np.array([[a.R.sum(), (~a.R).sum()], [b.R.sum(), (~b.R).sum()]]))[1]:.3f}")

# (2) gene-level within clade B
M = pd.read_csv(G / "mgit/gene_matrix.csv", index_col=0)
xb = b[b.NR.isin(M.index)]; y = xb.R.to_numpy(); sub = M.loc[xb.NR]
print(f"\n(2) gene-level, within clade B: {y.sum()} R vs {(~y).sum()} S+I")
pats: dict[tuple[Any, ...], list[tuple[Any, ...]]] = {}
for col in sub.columns:
    v = sub[col].to_numpy(); c1, c0 = int(v[y].sum()), int(v[~y].sum())
    t = [[c1, int(y.sum()) - c1], [c0, int((~y).sum()) - c0]]
    if c1 + c0 < 2 or (y.sum() - c1) + ((~y).sum() - c0) < 2: continue
    pats.setdefault(tuple(v), []).append((col, *t[0], *t[1], fisher_exact(np.array(t))[1]))
p = np.array([v[0][5] for v in pats.values()]); n = len(p); o = np.argsort(p); q = np.empty(n)
q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
rows = [dict(feature=f, R_car=a1, R_non=a2, nonR_car=b1, nonR_non=b2, p=pp, q=min(qq, 1), block=len(l)) for (pat, l), qq in zip(pats.items(), q) for f, a1, a2, b1, b2, pp in l]
res = pd.DataFrame(rows).sort_values("p"); res.to_csv(O / "cladeB_gene_tests.csv", index=False)
print(f"   patterns tested: {n}; min q = {res.q.min():.3f}")
pd.set_option("display.width", 220)
print(res.head(10).round(4).to_string(index=False))
kn = res[res.feature.str.contains("rrs|rpsL|eis|whiB7", case=False)]
print("   aminoglycoside loci:\n", kn.head(6).round(3).to_string(index=False) if len(kn) else "   none passed the minimum-carrier filter")

# (3) ANI within clade B
print(f"\n(3) ANI within clade B: R mean {b.ani[b.R].mean():.2f} vs S+I mean {b.ani[~b.R].mean():.2f}; MW p = {mannwhitneyu(b.ani[b.R], b.ani[~b.R]).pvalue:.3f}")

# (4) MIC cohort: clade vs MIC
ph = pd.read_csv(G / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded]
ph["clade"] = ph.NR.map(lin.set_index("NR").clade)
ph = ph[ph.clade.isin(["A", "B"])]
rows = []
for drug, d in ph.groupby("antibiotic"):
    d = d.drop_duplicates("NR"); A_, B_ = d[d.clade == "A"].log2_mic, d[d.clade == "B"].log2_mic
    rows.append(dict(drug=drug, nA=len(A_), nB=len(B_), medA=A_.median(), medB=B_.median(), p=mannwhitneyu(A_, B_).pvalue))
r = pd.DataFrame(rows).sort_values("p"); r["q_bh"] = np.minimum(1, (r.p * len(r) / (np.arange(len(r)) + 1))[::-1].cummin()[::-1])
r.to_csv(O / "mic_cohort_clade_vs_mic.csv", index=False)
print("\n(4) MIC cohort kansasii: clade A vs clade B, all drugs:\n", r.round(3).to_string(index=False))
