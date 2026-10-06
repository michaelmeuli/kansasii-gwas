"""Follow-up to script 11. The top MGIT amikacin-1 mg/L pattern (a block of genes identical across isolates) is a divergence-from-
reference signal. Characterise it (block membership, ANI, variant-set distances) and test ANI against amikacin MIC in the
independent MIC cohort (no isolate overlap for kansasii)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas/mgit"; RES = ROOT / "runs/mkan329/assembly/results"
KEY = "FCNCLM_04057||prom"
M = pd.read_csv(OUT / "gene_matrix.csv", index_col=0)
s = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "species", "gtdb_ani", "gtdb_reference"]).set_index("NR")
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv")
ak = mg[(mg.antibiotic == "Amikacin") & (mg.concentration_mg_l == 1.0)].drop_duplicates("NR").set_index("NR").erg
kan = [n for n in M.index if s.species[n] == "kansasii"]
block = [c for c in M.columns if (M[c] == M[KEY]).all()]
print(f"block of identical-pattern gene features among all 47 isolates: {len(block)}")
d = pd.DataFrame({"sample": s.PROBENNUMMER[kan], "ani": s.gtdb_ani[kan], "block": M.loc[kan, KEY], "amikacin_1mg": ak.reindex(kan)})
print("ANI reference genome:", s.gtdb_reference[kan].unique().tolist())
print(d.groupby("block").ani.describe().round(2).to_string())
k = d[d.amikacin_1mg.isin(["R", "I", "S"])].copy(); k["R"] = k.amikacin_1mg == "R"
print(f"ANI, amikacin 1 mg/L R vs S+I (n={k.R.sum()}/{(~k.R).sum()}): Mann-Whitney p = {mannwhitneyu(k.ani[k.R], k.ani[~k.R]).pvalue:.4f}")
d.to_csv(OUT / "kansasii_mgit_block_ani_amikacin.csv")

# variant-set distances among all non-excluded kansasii genomes
allk = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "species"])
allk = allk[(allk.species == "kansasii") & ~allk.NR.isin([35, 113, 106, 10, 21])]
S = {nr: set(map(tuple, pd.read_csv(RES / p / "5_typing/kansasii_snippy/snippy_out/snps.tab", sep="\t", dtype=str, usecols=["CHROM", "POS", "ALT"]).values))
     for nr, p in zip(allk.NR, allk.PROBENNUMMER)}
ids = list(S); D = np.array([[len(S[a] ^ S[b]) for b in ids] for a in ids]); ix = {n: i for i, n in enumerate(ids)}
cl = fcluster(linkage(squareform(D, checks=False), "average"), 1000, "distance")
print(f"kansasii genomes {len(ids)}: {len(set(cl))} clusters at 1000 differing sites; largest {pd.Series(cl).value_counts().max()}")
c = [n for n in kan if M.loc[n, KEY] == 1]; nc = [n for n in kan if M.loc[n, KEY] == 0]
def w(a: list[int], b: list[int]) -> float:
    return float(np.mean([D[ix[x], ix[y]] for x in a for y in b if x != y]))
print(f"mean differing sites: carriers-carriers {w(c, c):.0f}, non-non {w(nc, nc):.0f}, carriers-non {w(c, nc):.0f}")

# independent MIC cohort
ph = pd.read_csv(OUT.parent / "geno_pheno_table_long.csv"); ph = ph[ph.has_ngs & ~ph.excluded & (ph.gtdb_species == "kansasii")].copy()
ph["ani"] = ph.NR.map(s.gtdb_ani)
print(f"\nMIC cohort kansasii isolates: {ph.NR.nunique()}, also in MGIT: {ph[ph.NR.isin(set(mg.NR))].NR.nunique()}")
rows = []
for drug, g in ph.groupby("antibiotic"):
    g = g.drop_duplicates("NR"); rho, p = spearmanr(g.ani, g.log2_mic); rows.append((drug, len(g), round(rho, 2), round(p, 3)))
r = pd.DataFrame(rows, columns=["drug", "n", "spearman_ANI_vs_log2MIC", "p_unadj"]).sort_values("p_unadj"); r.to_csv(OUT / "mic_cohort_ani_vs_mic.csv", index=False)
print(r.to_string(index=False))
a = ph[ph.antibiotic == "Amikacin"].drop_duplicates("NR"); lo, hi = a[a.ani < 99.67], a[a.ani >= 99.67]
print(f"Amikacin log2 MIC: ANI<99.67 n={len(lo)} median {lo.log2_mic.median()} vs ANI>=99.67 n={len(hi)} median {hi.log2_mic.median()}; MW p = {mannwhitneyu(lo.log2_mic, hi.log2_mic).pvalue:.3f}")
