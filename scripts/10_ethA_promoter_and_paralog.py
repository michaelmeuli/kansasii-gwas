"""(a) Verify 054's paralog (FCNCLM_03714) calls against its own assembly. (b) Test the ethA upstream change at 3145420
(A>C) against ethionamide MIC across all isolates with MIC and NGS."""
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Seq import Seq
from scipy.stats import fisher_exact, mannwhitneyu

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas"; RES = ROOT / "runs/mkan329/assembly/results"
REFDIR = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation"
genome = {r.id: str(r.seq) for r in SeqIO.parse(REFDIR / "GCF_000157895.3_query.fna", "fasta")}["contig_1"]
meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "species", "NR"]).set_index("PROBENNUMMER")
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5

# (a) paralog protein from own assembly: 054 vs other persicum
s = genome[4126237:4127734]; refp = str(Seq(str(Seq(s).reverse_complement())).translate(table=11))
pers = [x for x in meta.index if meta.loc[x, "species"] == "persicum" and meta.loc[x, "NR"] not in (35, 113, 106, 10, 21)]
prots = {}
for x in pers:
    best, bs = None, -1
    for rec in SeqIO.parse(RES / x / f"2_annotation/{x}.faa", "fasta"):
        if 400 < len(rec.seq) < 600:
            sc = al.score(refp, str(rec.seq))
            if sc > bs: best, bs = str(rec.seq), sc
    prots[x] = best
def diffs(a, b):
    aln = al.align(a, b)[0]; d = []
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        d += [f"{a[a0+i]}{a0+i+1}{b[b0+i]}" for i in range(a1 - a0) if a[a0 + i] != b[b0 + i]]
    return d
print("paralog protein lengths:", {k[-3:]: len(v) for k, v in prots.items()})
others = [x for x in pers if x != "Mkan329-054"]
ref_other = prots["Mkan329-009"]
print("054 vs 009 (persicum) paralog protein differences:", diffs(ref_other, prots["Mkan329-054"]))
print("002 vs 009 for scale:", diffs(ref_other, prots["Mkan329-002"]))
print("kansasii-ref paralog vs 009:", len(diffs(refp, ref_other)), "differences")

# (b) 3145420 A>C across all genomes with MIC
pheno = pd.read_csv(OUT / "geno_pheno_table_long.csv"); pheno = pheno[pheno.has_ngs & ~pheno.excluded & pheno.antibiotic.eq("Ethionamid")].drop_duplicates("NR")
carr = {}
for x in pheno.PROBENNUMMER:
    t = pd.read_csv(RES / x / "5_typing/kansasii_snippy/snippy_out/snps.tab", sep="\t", dtype=str, usecols=["CHROM", "POS", "REF", "ALT"])
    t["POS"] = t.POS.astype(int); h = t[(t.CHROM == "NC_022663.1") & (t.POS == 3145420)]
    carr[x] = ("; ".join(f"{r.REF}>{r.ALT}" for r in h.itertuples()) or "ref")
pheno["state_3145420"] = pheno.PROBENNUMMER.map(carr)
print("\nstate at 3145420 (ethA upstream) among isolates with ethionamide MIC:\n", pheno.groupby(["gtdb_species", "state_3145420"]).agg(n=("NR", "size"), median_log2=("log2_mic", "median")).to_string())
c = pheno[pheno.state_3145420 == "A>C"]; n = pheno[pheno.state_3145420 == "ref"]
print("\ncarriers:", c[["PROBENNUMMER", "gtdb_species", "mhk_parsed", "log2_mic"]].to_string(index=False))
print("non-carriers n =", len(n), "median log2", n.log2_mic.median())
print("Mann-Whitney p (all species, unadjusted):", mannwhitneyu(c.log2_mic, n.log2_mic).pvalue)
pp = pheno[pheno.gtdb_species == "persicum"]
print("persicum only:\n", pp[["PROBENNUMMER", "mhk_parsed", "log2_mic", "state_3145420"]].to_string(index=False))
pheno.to_csv(OUT / "ethionamide_ethA_upstream_3145420.csv", index=False)
