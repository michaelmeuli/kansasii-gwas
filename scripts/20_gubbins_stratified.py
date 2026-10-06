"""Re-do the lineage stratification with the Gubbins recombination-filtered tree (final run) instead of raw SNP distances.
The filtered tree is shallow and star-like, so instead of two big clades we use (a) clusters of close relatives as strata (Mantel-Haenszel),
(b) Mantel-type tests of phenotype distance vs tree distance, (c) how the region haplotype is distributed on the tree, and
(d) gene-level tests within the main cluster. Mkan329-068 is excluded: its 28,959-SNP terminal branch is a tree artifact (raw nearest-neighbour distance
is normal, no private Gubbins blocks)."""
from __future__ import annotations

import os
from typing import Any
import numpy.typing as npt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.contingency_tables import StratifiedTable

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); G = ROOT / "output/gwas"; OUT = G / "gubbins"
rng = np.random.default_rng(7)
D = pd.read_csv(OUT / "filtered_tree_distance_matrix.csv", index_col=0).drop(index="Mkan329-068", columns="Mkan329-068")
names = D.index.tolist(); print("isolates in tree analyses:", len(names), "(Mkan329-068 excluded)")
Z = linkage(squareform(D.values, checks=False), "average")
cl = {h: pd.Series(fcluster(Z, h, "distance"), index=names) for h in (40, 60, 100)}
for h, c in cl.items(): print(f"  cut {h} SNPs: {c.nunique()} clusters, largest {c.value_counts().iloc[0]}")
hap = pd.read_csv(G / "lineage/kansasii_region_haplotype.csv", index_col=0).loc[names]
meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER"]).set_index("PROBENNUMMER").NR
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = hap[hap.NR.isin(ak.index)].copy(); x["R"] = ak.reindex(x.NR).values == "R"; x["H880"] = x.hap == "H~880"
for h in cl: x[f"c{h}"] = cl[h].reindex(x.index)
print(f"\nphenotyped for amikacin 1 mg/L: {len(x)} (R {int(x.R.sum())}, S/I {int((~x.R).sum())}); haplotype-880: {int(x.H880.sum())}")

# (a) stratified (Mantel-Haenszel) test: haplotype-880 vs amikacin R, strata = clusters of close relatives
print("\n(a) Mantel-Haenszel, strata = clusters of close relatives")
for h in (40, 60, 100):
    tabs, info = [], []
    for cid, g in x.groupby(f"c{h}"):
        t = np.array([[(g.H880 & g.R).sum(), (g.H880 & ~g.R).sum()], [(~g.H880 & g.R).sum(), (~g.H880 & ~g.R).sum()]])
        if t[0].sum() > 0 and t[1].sum() > 0 and t[:, 0].sum() > 0 and t[:, 1].sum() > 0:   # informative: exposure and outcome both vary
            tabs.append(t); info.append((cid, t.tolist()))
    print(f"  cut {h}: strata {x[f'c{h}'].nunique()}, informative strata {len(tabs)}")
    if tabs:
        print("  (the asymptotic Mantel-Haenszel p is unreliable with one stratum of 5-13 isolates, so exact Fisher p per stratum is reported)")
        for cid, t in info: print(f"     stratum {cid}: carriers R/S-I {t[0]}, non-carriers R/S-I {t[1]}; exact Fisher p = {fisher_exact(np.array(t))[1]:.3f}")
        if len(tabs) > 1:
            st = StratifiedTable([t + 0.0 for t in tabs]); print(f"     pooled OR {st.oddsratio_pooled:.2f}, MH p = {st.test_null_odds(correction=False).pvalue:.3f}")
    else: print("     no informative stratum, the test cannot be done")

# (b) Mantel-type: phenotype distance vs tree distance
def mantel(ids: list[str], y: npt.NDArray[np.floating[Any]], binary: bool = False, n: int = 5000) -> tuple[float, float]:
    iu = np.triu_indices(len(ids), 1); dd = D.loc[ids, ids].values[iu]
    py = (np.abs(y[:, None] - y[None, :]) > 0 if binary else np.abs(y[:, None] - y[None, :]))[iu].astype(float)
    obs = spearmanr(dd, py).statistic
    perm: list[float] = []
    for _ in range(n):
        p = rng.permutation(len(ids)); yy = y[p]
        py2 = (np.abs(yy[:, None] - yy[None, :]) > 0 if binary else np.abs(yy[:, None] - yy[None, :]))[iu].astype(float)
        perm.append(spearmanr(dd, py2).statistic)
    return obs, (1 + sum(abs(v) >= abs(obs) for v in perm)) / (n + 1)
print("\n(b) Mantel-type test: do closely related isolates have more similar phenotypes? (Spearman of |phenotype difference| vs tree distance; +ve = lineage signal)")
ph = pd.read_csv(G / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded]; ph = ph[ph.PROBENNUMMER.isin(names)]
rows: list[dict[str, Any]] = []
for drug, g in ph.groupby("antibiotic"):
    g = g.drop_duplicates("PROBENNUMMER"); ids = g.PROBENNUMMER.tolist()
    if len(ids) < 8: continue
    rho_d, p_d = mantel(ids, g.log2_mic.to_numpy(dtype=float)); rows.append(dict(drug=drug, n=len(ids), rho=round(rho_d, 3), p_perm=round(p_d, 4)))
r = pd.DataFrame(rows).sort_values("p_perm"); r["q_bh"] = np.minimum(1, (r.p_perm * len(r) / (np.arange(len(r)) + 1))[::-1].cummin()[::-1]).round(3)
print(r.to_string(index=False)); r.to_csv(OUT / "mantel_mic_vs_filtered_tree.csv", index=False)
rho, p = mantel(x.index.tolist(), x.R.astype(float).to_numpy(dtype=float), binary=True)
print(f"  MGIT amikacin 1 mg/L (R vs S/I, n={len(x)}): rho = {rho:.3f}, permutation p = {p:.4f}")

# (c) how the haplotype is distributed on the tree
carr = hap.index[hap.hap == "H~880"]
print(f"\n(c) haplotype-880 carriers: {len(carr)} of {len(names)}")
for h in (40, 60):
    cc = cl[h].reindex(carr); print(f"  cut {h}: carriers fall in {cc.nunique()} different clusters (cluster sizes of those: {sorted(cl[h].value_counts()[cc.unique()].tolist(), reverse=True)[:10]})")
nn = D.where(~np.eye(len(D), dtype=bool)).idxmin(axis=1)
obs = np.mean([nn[c] in set(carr) for c in carr])
exp = np.mean([np.mean([nn[c] in set(rng.choice(names, len(carr), replace=False)) for c in rng.choice(names, len(carr), replace=False)]) for _ in range(2000)])
print(f"  fraction of carriers whose nearest relative is also a carrier: {obs:.2f} (random expectation about {exp:.2f})")

# (d) gene-level, within the main cluster
M = pd.read_csv(G / "mgit/gene_matrix.csv", index_col=0)
main = cl[100].value_counts().index[0]
xm = x[x.c100 == main]; xm = xm[xm.NR.isin(M.index)]
y = xm.R.to_numpy(); sub = M.loc[xm.NR]
print(f"\n(d) gene-level within the main cluster (cut 100; {int((cl[100] == main).sum())} isolates): phenotyped {len(xm)} (R {y.sum()}, S/I {(~y).sum()})")
pats: dict[tuple[Any, ...], list[tuple[Any, ...]]] = {}
for col in sub.columns:
    v = sub[col].to_numpy(); c1, c0 = int(v[y].sum()), int(v[~y].sum()); t2 = [[c1, int(y.sum()) - c1], [c0, int((~y).sum()) - c0]]
    if c1 + c0 < 2 or (y.sum() - c1) + ((~y).sum() - c0) < 2: continue
    pats.setdefault(tuple(v), []).append((col, *t2[0], *t2[1], fisher_exact(np.array(t2))[1]))
if pats:
    pa = np.array([v[0][5] for v in pats.values()]); n = len(pa); o = np.argsort(pa); q = np.empty(n)
    q[o] = np.minimum.accumulate((pa[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    res = pd.DataFrame([dict(feature=f, R_car=a, R_non=b, nonR_car=c, nonR_non=d, p=pp, q=min(qq, 1), block=len(l)) for (pat, l), qq in zip(pats.items(), q) for f, a, b, c, d, pp in l]).sort_values("p")
    res.to_csv(OUT / "main_cluster_gene_tests.csv", index=False)
    print(f"  patterns tested {n}; min q = {res.q.min():.3f}"); print(res.head(6).round(4).to_string(index=False))
