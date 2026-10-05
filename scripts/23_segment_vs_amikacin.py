"""Segment classes (script 22, persicum-like threshold relaxed to >= 0.75; none <= 0.10) vs amikacin: MGIT 1 mg/L, MIC cohort, and tree-aware checks."""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import fisher_exact, mannwhitneyu, binomtest

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"
d = pd.read_csv(G / "boundaries/segment_classes.csv", index_col=0)
def cls(r):
    L, R = r.frac_persicum_left, r.frac_persicum_right
    if L >= .75 and R >= .75: return "full"
    if L >= .75 and R <= .10: return "left-only"
    if L <= .10 and R <= .10: return "none"
    return "partial/mosaic"
d["cls"] = d.apply(cls, axis=1); d.to_csv(G / "boundaries/segment_classes_v2.csv")
print("segment classes (82 kansasii):", d.cls.value_counts().to_dict())
print("old haplotype class vs new:\n", pd.crosstab(d.hap.astype(str), d.cls).to_string())
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = d[d.NR.isin(ak.index)].copy(); x["R"] = ak.reindex(x.NR).values == "R"
print(f"\nMGIT amikacin 1 mg/L, kansasii (n={len(x)}; R {int(x.R.sum())}):\n", pd.crosstab(x.cls, x.R).rename(columns={True: "R", False: "S/I"}).to_string())
def fisher(mask, label):
    a, b, c, e = int((mask & x.R).sum()), int((mask & ~x.R).sum()), int((~mask & x.R).sum()), int((~mask & ~x.R).sum())
    print(f"  {label}: R {a}/{a+b} vs {c}/{c+e}; Fisher p = {fisher_exact([[a, b], [c, e]])[1]:.4f}")
fisher(x.cls == "full", "full segment vs rest")
fisher(x.cls.isin(["full", "left-only"]), "any persicum-like piece (full or left-only) vs none/partial")
fisher(x.cls.isin(["full"]) | (x.frac_persicum_right >= .75), "right part persicum-like vs rest")
fisher(x.cls == "left-only", "left-only vs rest")
ph = pd.read_csv(G / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded & (ph.antibiotic == "Amikacin")].drop_duplicates("NR")
ph["cls"] = ph.NR.map(d.reset_index().drop_duplicates("NR").set_index("NR").cls); ph = ph.dropna(subset=["cls"])
print("\nMIC cohort (21 kansasii) amikacin log2 MIC by class:\n", ph.groupby("cls").log2_mic.agg(["count", "median", "min", "max"]).to_string())
a, b = ph[ph.cls == "full"].log2_mic, ph[ph.cls != "full"].log2_mic
print(f"  full (n={len(a)}) vs others (n={len(b)}): Mann-Whitney p = {mannwhitneyu(a, b).pvalue:.3f}" if len(a) and len(b) else "  too few")
a2 = ph[ph.cls.isin(["full", "left-only"])].log2_mic; b2 = ph[~ph.cls.isin(["full", "left-only"])].log2_mic
print(f"  full or left-only (n={len(a2)}) vs others (n={len(b2)}): Mann-Whitney p = {mannwhitneyu(a2, b2).pvalue:.3f}")

# tree-aware: nearest-relative discordant pairs, Gubbins filtered tree (Mkan329-068 excluded)
D = pd.read_csv(G / "gubbins/filtered_tree_distance_matrix.csv", index_col=0).drop(index="Mkan329-068", columns="Mkan329-068")
xx = x[x.index.isin(D.index)].copy(); ids = xx.index.tolist(); Dm = D.loc[ids, ids]
nn = {i: Dm.loc[i].drop(i).idxmin() for i in ids}; pairs = {tuple(sorted((i, j))) for i, j in nn.items()}
for label, flag in (("full segment", xx.cls == "full"), ("any persicum-like piece", xx.cls.isin(["full", "left-only"]))):
    disc = [(a, b) for a, b in pairs if flag[a] != flag[b]]; sup = sum((flag[a] and xx.R[a] and not xx.R[b]) or (flag[b] and xx.R[b] and not xx.R[a]) for a, b in disc)
    rev = sum((flag[a] and not xx.R[a] and xx.R[b]) or (flag[b] and not xx.R[b] and xx.R[a]) for a, b in disc)
    print(f"\nnearest-relative pairs discordant for '{label}': {len(disc)}; carrier R and non-carrier S/I: {sup}; carrier S/I and non-carrier R: {rev}")
    for a_, b_ in disc: print(f"    {a_} ({xx.cls[a_]}, {'R' if xx.R[a_] else 'S/I'}) vs {b_} ({xx.cls[b_]}, {'R' if xx.R[b_] else 'S/I'}), distance {Dm.loc[a_, b_]:.0f}")
    if len(disc): print(f"    sign test p = {binomtest(sup, sup + rev, .5).pvalue:.3f}" if sup + rev else "    no resolvable pairs")
Z = linkage(squareform(D.values, checks=False), "average")
for cut in (60, 100):
    c = pd.Series(fcluster(Z, cut, "distance"), index=D.index); xx[f"c{cut}"] = c.reindex(xx.index)
    for label, flag in (("full", xx.cls == "full"), ("full or left-only", xx.cls.isin(["full", "left-only"]))):
        inf = []
        for k, g in xx.groupby(f"c{cut}"):
            f = flag[g.index]; r = g.R
            if f.any() and (~f).any() and r.any() and (~r).any(): inf.append((k, int((f & r).sum()), int((f & ~r).sum()), int((~f & r).sum()), int((~f & ~r).sum())))
        print(f"\ncut {cut} SNPs, '{label}': informative strata {len(inf)}", end="")
        for k, a, b, c_, e in inf: print(f"\n    stratum {k}: carriers R/S-I {a}/{b}, non-carriers {c_}/{e}; exact Fisher p = {fisher_exact([[a, b], [c_, e]])[1]:.3f}", end="")
        print()
