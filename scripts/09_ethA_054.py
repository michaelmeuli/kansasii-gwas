"""Which reference locus is the true EthA ortholog, and what does Mkan329-054 carry there (read evidence + own assembly)?"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Seq import Seq

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas"
RES = ROOT / "runs/mkan329/assembly/results"
REFDIR = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation"
genome = {r.id: str(r.seq) for r in SeqIO.parse(REFDIR / "GCF_000157895.3_query.fna", "fasta")}["contig_1"]
G = pd.read_csv(OUT / "candidate_gene_coords.csv"); E = G[(G.gene == "ethA") & (G.chrom == "NC_022663.1")]
mtb = "".join(l.strip() for l in open(OUT / "ref/Mtb_ethA.faa") if not l.startswith(">"))
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
def ident(a, b):
    aln = al.align(a, b)[0]; n = m = 0
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        n += a1 - a0; m += sum(a[a0 + i] == b[b0 + i] for i in range(a1 - a0))
    return m / len(a), aln
prot = {}
print(f"Mtb EthA length {len(mtb)}")
for r in E.itertuples():
    s = genome[r.start - 1:r.end]; s = s if r.strand == "+" else str(Seq(s).reverse_complement())
    p = str(Seq(s).translate(table=11)); prot[r.locus] = (p, r)
    i, _ = ident(mtb, p)
    print(f"  {r.locus} {r.start}-{r.end} {r.strand} len {len(p)} aa: identity to Mtb EthA over Mtb length = {i:.2f}")
best = max(prot, key=lambda k: ident(mtb, prot[k][0])[0])
print("true EthA ortholog:", best)
refp, r = prot[best]
_, aln = ident(mtb, refp)
m = {}
for (a0, a1), (b0, b1) in zip(*aln.aligned):
    for i in range(a1 - a0): m[b0 + i + 1] = a0 + i + 1   # kansasii -> Mtb

# 1. read evidence in snps.tab for 054 and the other persicum genomes across the ortholog
meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "species", "NR"]).set_index("PROBENNUMMER")
pers = [s for s in meta.index if meta.loc[s, "species"] == "persicum" and meta.loc[s, "NR"] not in (35, 113, 106, 10, 21)]
print(f"\npersicum genomes (non-excluded): {len(pers)}")
calls = {}
for s in pers:
    t = pd.read_csv(RES / s / "5_typing/kansasii_snippy/snippy_out/snps.tab", sep="\t", dtype=str)
    t["POS"] = t.POS.astype(int)
    calls[s] = t[(t.CHROM == "NC_022663.1") & (t.POS >= r.start) & (t.POS <= r.end)]
t54 = calls["Mkan329-054"]
print(f"\nMkan329-054 calls inside {best} ({len(t54)}):")
print(t54[["POS", "TYPE", "REF", "ALT", "EVIDENCE"]].to_string(index=False))
near = t54[(t54.POS >= 4126230) & (t54.POS <= 4126270)]
print("\nsame positions in the other persicum genomes:")
for s in pers:
    if s == "Mkan329-054": continue
    c = calls[s]; c = c[(c.POS >= 4126230) & (c.POS <= 4126270)]
    print(" ", s, "|", "; ".join(f"{x.POS}:{x.REF}>{x.ALT}[{x.EVIDENCE}]" for x in c.itertuples()) or "no calls (= reference-like)")

# 2. own assembly: protein of each persicum genome best matching the ortholog
print("\nEthA protein from each genome's own Bakta annotation vs the kansasii reference ortholog:")
rows = []
for s in pers:
    best_rec, bs = None, -1
    for rec in SeqIO.parse(RES / s / f"2_annotation/{s}.faa", "fasta"):
        if 400 < len(rec.seq) < 600 and "monooxygenase" in rec.description.lower():
            sc = al.score(refp, str(rec.seq))
            if sc > bs: best_rec, bs = rec, sc
    if best_rec is None: rows.append((s, None)); continue
    p = str(best_rec.seq); i, a2 = ident(refp, p)
    mm = []
    for (a0, a1), (b0, b1) in zip(*a2.aligned):
        for k in range(a1 - a0):
            if refp[a0 + k] != p[b0 + k]: mm.append(f"{refp[a0+k]}{a0+k+1}{p[b0+k]}")
    rows.append((s, (len(p), round(i, 3), mm)))
for s, v in rows:
    print(" ", s, "no match" if v is None else f"len {v[0]} id {v[1]} diffs vs ref ({len(v[2])}): {', '.join(v[2][:15])}")
