"""Gene-specific version of the amikacin hypothesis: persicum-like RsmA or RsmI allele (vs kansasii reference-like) against amikacin, MGIT and MIC cohorts.
Alleles from script 24 (any non-reference allele counts as persicum-like; all non-reference alleles derive from the persicum-like haplotype)."""
import os
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); G = ROOT / "output/gwas"; B = G / "boundaries"
seg = pd.read_csv(B / "segment_classes_v2.csv", index_col=0)
for g in ("RsmA", "RsmI"): seg[g] = pd.read_csv(B / f"{g}_allele_per_isolate.csv", index_col=0).iloc[:, 0].ne("REF-like")
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0) & mg.erg.isin(["R", "I", "S"])].drop_duplicates("NR").set_index("NR").erg
x = seg[seg.NR.isin(ak.index)].copy(); x["R"] = ak.reindex(x.NR).values == "R"
ph = pd.read_csv(G / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded & (ph.antibiotic == "Amikacin")].drop_duplicates("NR")
m = seg.reset_index().drop_duplicates("NR").set_index("NR"); ph = ph[ph.NR.isin(m.index)].copy()
for g in ("RsmA", "RsmI"): ph[g] = ph.NR.map(m[g])
print(f"MGIT amikacin 1 mg/L, kansasii (n={len(x)}):")
for g in ("RsmA", "RsmI"):
    t = pd.crosstab(x[g], x.R).reindex(index=[True, False], columns=[True, False], fill_value=0)
    print(f"  persicum-like {g}: R {t.loc[True, True]}/{t.loc[True].sum()} vs ref-like R {t.loc[False, True]}/{t.loc[False].sum()}; Fisher p = {fisher_exact(t.values)[1]:.4f}")
print(f"\nMIC cohort amikacin log2 MIC (kansasii, n={len(ph)}):")
for g in ("RsmA", "RsmI"):
    a, b = ph[ph[g]].log2_mic, ph[~ph[g]].log2_mic
    print(f"  persicum-like {g} (n={len(a)}) median {a.median()} vs ref-like (n={len(b)}) median {b.median()}; Mann-Whitney p = {mannwhitneyu(a, b).pvalue:.3f}" if len(a) and len(b) else f"  {g}: too few")
