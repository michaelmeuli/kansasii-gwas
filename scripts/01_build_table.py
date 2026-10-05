"""Join MIC phenotypes with species/NGS metadata; flag exclusions. No patient-identifying columns kept."""
import pandas as pd
from pathlib import Path

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii")
OUT = ROOT / "output/gwas"
EXCLUDE = {35, 113, 106, 10, 21}      # mixed / contaminated cultures
REVIEW = {8, 28, 46}                  # species calls disagree
CONTROL_SPECIES = {"Mycobacterium tuberculosis", "Mycobacterium avium", "Mycobacterium abscessus",
                   "Mycobacterium xenopi", "Mycobacterium marinum"}

mic = pd.read_csv(ROOT / "output/mic/mhk/mic_parsed.csv")
meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv",
                   usecols=["NR", "PROBENNUMMER", "species", "contamination_flag", "checkm_contamination"])
meta = meta.rename(columns={"species": "gtdb_species"})

df = mic.merge(meta.drop(columns="PROBENNUMMER"), on="NR", how="left")
df["has_ngs"] = df["gtdb_species"].notna()
df["excluded"] = df["NR"].isin(EXCLUDE) | df["gtdb_species"].isin(CONTROL_SPECIES)
df["review"] = df["NR"].isin(REVIEW)
keep = ["NR", "PROBENNUMMER", "TNR", "antibiotic", "mhk_parsed", "log2_mic", "censored", "int_erg",
        "gtdb_species", "checkm_contamination", "has_ngs", "excluded", "review"]
df = df[keep]
df.to_csv(OUT / "geno_pheno_table_long.csv", index=False)

use = df[df.has_ngs & ~df.excluded]
print("rows", len(df), "| isolates with MIC", df.NR.nunique(), "| usable (NGS, not excluded)", use.NR.nunique())
print(use.drop_duplicates("NR").gtdb_species.value_counts().to_string())
print(use.groupby("antibiotic").agg(n=("NR", "nunique"), min_log2=("log2_mic", "min"),
                                    max_log2=("log2_mic", "max"), sd=("log2_mic", "std")).round(2).to_string())
