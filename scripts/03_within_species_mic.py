"""Within-species carrier vs non-carrier MIC comparison for candidate-gene variants.
Features: codon-level amino-acid changes (multi-base snippy calls split per base and merged per codon),
rRNA nucleotide changes, and per-gene frameshift/indel flags. Carrier = differs from reference at that feature;
polarity is irrelevant here because only sites segregating inside a species are tested."""
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import mannwhitneyu

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii")
OUT = ROOT / "output/gwas"
RES = ROOT / "runs/mkan329/assembly/results"
REFFNA = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.fna"
CONTIG = {"contig_1": "NC_022663.1", "contig_2": "NC_022654.1"}
genome = {CONTIG[r.id]: str(r.seq) for r in SeqIO.parse(REFFNA, "fasta")}
genes = pd.read_csv(OUT / "candidate_gene_coords.csv")
genes["id"] = genes.gene + "|" + genes.locus
RELEVANT = {"Ciprofloxacin": ["gyrA", "gyrB"], "Moxifloxacin": ["gyrA", "gyrB"], "Rifampicin": ["rpoB"],
            "Rifabutin": ["rpoB"], "Ethambutol": ["embA", "embB", "embC"], "Isoniazid": ["katG", "inhA"],
            "Ethionamid": ["ethA", "ethR", "inhA"], "Streptomycin": ["rpsL", "rrs"], "Amikacin": ["rrs"],
            "Clarithromycin": ["rrl"], "Linezolid": ["rplC", "rrl"], "Sulfamethoxazole/Trimethoprim": ["folP"]}

def per_base(row):
    """Yield (pos, ref, alt) single-base substitutions; None for length-changing calls."""
    if len(row.REF) != len(row.ALT):
        return None
    return [(row.POS + i, r, a) for i, (r, a) in enumerate(zip(row.REF, row.ALT)) if r != a]

feats = {}   # sample -> set of feature strings
for tab in sorted(RES.glob("Mkan329-*/5_typing/kansasii_snippy/snippy_out/snps.tab")):
    sample = tab.parts[-5]
    t = pd.read_csv(tab, sep="\t", dtype=str, usecols=["CHROM", "POS", "TYPE", "REF", "ALT"])
    t["POS"] = t.POS.astype(int)
    fs = set()
    for g in genes.itertuples():
        sub = t[(t.CHROM == g.chrom) & (t.POS >= g.start) & (t.POS <= g.end)]
        if sub.empty:
            continue
        codon_subs = {}
        for v in sub.itertuples():
            subs = per_base(v)
            if subs is None:
                if g.type == "CDS":
                    fs.add(f"{g.id}:indel")
                else:
                    fs.add(f"{g.id}:indel@{v.POS}")
                continue
            for pos, r, a in subs:
                if not g.start <= pos <= g.end:   # multi-base call overhanging the gene
                    continue
                if g.type == "rRNA":
                    fs.add(f"{g.id}:{r}{pos}{a}")
                else:
                    off = pos - g.start if g.strand == "+" else g.end - pos
                    if g.strand == "-":
                        r, a = str(Seq(r).reverse_complement()), str(Seq(a).reverse_complement())
                    codon_subs.setdefault(off // 3, []).append((off % 3, a))
        if g.type == "CDS" and codon_subs:
            cds = genome[g.chrom][g.start - 1:g.end]
            cds = cds if g.strand == "+" else str(Seq(cds).reverse_complement())
            for ci, lst in codon_subs.items():
                c0 = cds[ci * 3:ci * 3 + 3]
                c1 = list(c0)
                for k, a in lst:
                    c1[k] = a
                a0, a1 = str(Seq(c0).translate(table=11)), str(Seq("".join(c1)).translate(table=11))
                if a0 != a1:
                    fs.add(f"{g.id}:{a0}{ci + 1}{a1}")
    feats[sample] = fs

pheno = pd.read_csv(OUT / "geno_pheno_table_long.csv")
pheno = pheno[pheno.has_ngs & ~pheno.excluded]
all_feats = sorted({f for s in feats.values() for f in s})
print("samples", len(feats), "| distinct candidate features", len(all_feats))

rows = []
for sp_label, sdf in [("kansasii", pheno[pheno.gtdb_species == "kansasii"]),
                      ("persicum", pheno[pheno.gtdb_species == "persicum"]),
                      ("pooled_species_centred", pheno)]:
    for drug, d in sdf.groupby("antibiotic"):
        d = d.drop_duplicates("NR").copy()
        if sp_label.startswith("pooled"):
            d["y"] = d.log2_mic - d.groupby("gtdb_species").log2_mic.transform("mean")
        else:
            d["y"] = d.log2_mic
        for f in all_feats:
            carr = d.PROBENNUMMER.map(lambda s: f in feats.get(s, set())).astype(bool)
            n1, n0 = int(carr.sum()), int((~carr).sum())
            if n1 < 2 or n0 < 2:
                continue
            y1, y0 = d.y[carr], d.y[~carr]
            p = mannwhitneyu(y1, y0, alternative="two-sided").pvalue
            gene = f.split("|")[0]
            rows.append(dict(subset=sp_label, drug=drug, feature=f, n_carrier=n1, n_noncarrier=n0,
                             median_carrier=y1.median(), median_noncarrier=y0.median(),
                             delta_log2=y1.mean() - y0.mean(), p=p,
                             relevant_gene=gene in RELEVANT.get(drug, [])))
res = pd.DataFrame(rows)

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.minimum(q, 1)

res["q_bh_subset"] = res.groupby("subset").p.transform(lambda s: bh(s))
res = res.sort_values("p")
res.to_csv(OUT / "within_species_mic_tests.csv", index=False)
pd.DataFrame([(s, sorted(v)) for s, v in feats.items()], columns=["sample", "features"]).to_csv(OUT / "candidate_features_per_sample.csv", index=False)
print("tests:", len(res), res.groupby("subset").size().to_dict())
pd.set_option("display.width", 220)
cols = ["subset", "drug", "feature", "n_carrier", "n_noncarrier", "delta_log2", "p", "q_bh_subset", "relevant_gene"]
print("\nTop 25 overall:\n", res[cols].head(25).round(4).to_string(index=False))
print("\nTop hits in biologically relevant gene-drug pairs:\n",
      res[res.relevant_gene][cols].head(20).round(4).to_string(index=False))
