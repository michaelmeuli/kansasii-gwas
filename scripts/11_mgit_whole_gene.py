"""MGIT whole-gene association. Per gene: carries >=1 protein-altering variant (aa change, frameshift, inframe indel, stop),
any rRNA/ncRNA variant, or any variant in the 150 bp upstream region. Fisher exact within species; BH over distinct
carrier patterns; block size = number of genes sharing a pattern (large block = lineage marker, not a gene-specific effect)."""
from __future__ import annotations

import os
import re
from collections.abc import Iterator
from typing import Any, cast
import numpy as np
import numpy.typing as npt
import pandas as pd
from bisect import bisect_right
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import fisher_exact

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); OUT = ROOT / "output/gwas/mgit"; OUT.mkdir(exist_ok=True)
RES = ROOT / "runs/mkan329/assembly/results"
REFDIR = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation"
CONTIG = {"contig_1": "NC_022663.1", "contig_2": "NC_022654.1"}
EXCL = {35, 113, 106, 10, 21}
CTRL = {"Mycobacterium tuberculosis", "Mycobacterium avium", "Mycobacterium abscessus", "Mycobacterium xenopi", "Mycobacterium marinum"}
UP = 150
genome = {CONTIG[r.id]: str(r.seq) for r in SeqIO.parse(REFDIR / "GCF_000157895.3_query.fna", "fasta")}

# ---- annotation
gene_rows: list[dict[str, Any]] = []
for l in open(REFDIR / "GCF_000157895.3_query.gff3", encoding="utf-8"):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("CDS", "rRNA", "tRNA", "ncRNA"): continue
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    gene_rows.append(dict(chrom=CONTIG[f[0]], start=int(f[3]), end=int(f[4]), strand=f[6], ftype=f[2],
                      gene=attrs.get("gene", ""), locus=attrs.get("locus_tag", ""), product=attrs.get("product", "")))
genes = pd.DataFrame(gene_rows).sort_values(["chrom", "start"]).reset_index(drop=True)
cds_seq: dict[int, str] = {}
for i_, g in genes.iterrows():
    i = cast(int, i_)
    if g.ftype == "CDS":
        s = genome[g.chrom][g.start - 1:g.end]; cds_seq[i] = s if g.strand == "+" else str(Seq(s).reverse_complement())
idx = {c: (genes[genes.chrom == c].index.to_numpy(), genes[genes.chrom == c].start.to_numpy()) for c in genes.chrom.unique()}
def coding_hit(chrom: str, pos: int) -> int | None:
    ids, starts = idx[chrom]; k = bisect_right(starts, pos) - 1
    for j in range(k, max(k - 4, -1), -1):
        g = genes.loc[ids[j]]
        if g.start <= pos <= g.end: return int(ids[j])
    return None
prom: defaultdict[str, list[tuple[int, int, int]]] = defaultdict(list)   # chrom -> (lo, hi, gene index) for CDS upstream regions
for i_, g in genes.iterrows():
    i = cast(int, i_)
    if g.ftype != "CDS": continue
    lo, hi = (g.start - UP, g.start - 1) if g.strand == "+" else (g.end + 1, g.end + UP)
    prom[g.chrom].append((lo, hi, i))
for ch in prom: prom[ch].sort()
plo = {c: np.array([x[0] for x in v]) for c, v in prom.items()}
def prom_hit(chrom: str, pos: int) -> int | None:
    k = bisect_right(plo[chrom], pos) - 1
    for j in range(k, max(k - 6, -1), -1):
        lo, hi, i = prom[chrom][j]
        if lo <= pos <= hi: return i
    return None

# ---- per-isolate gene calls
def gene_calls(snps: pd.DataFrame) -> set[tuple[int, str]]:
    """snps: DataFrame CHROM,POS,TYPE,REF,ALT -> set of (gene_idx, kind) with kind in prot/rna/prom"""
    out: set[tuple[int, str]] = set()
    codon: defaultdict[int, defaultdict[int, list[tuple[int, str]]]] = defaultdict(lambda: defaultdict(list))
    for chrom, pos, typ, ref, alt in snps.itertuples(index=False):
        h = coding_hit(chrom, pos)
        if h is None:
            p = prom_hit(chrom, pos)
            if p is not None: out.add((p, "prom"))
            continue
        g = genes.loc[h]
        if g.ftype != "CDS": out.add((h, "rna")); continue
        if len(ref) != len(alt):
            out.add((h, "prot")); continue          # frameshift or inframe indel: both protein-altering
        for i, (r, a) in enumerate(zip(ref, alt)):
            q = pos + i
            if r == a or not g.start <= q <= g.end: continue
            off = q - g.start if g.strand == "+" else g.end - q
            a2 = a if g.strand == "+" else str(Seq(a).reverse_complement())
            codon[h][off // 3].append((off % 3, a2))
    for h, d in codon.items():
        for ci, lst in d.items():
            c0 = cds_seq[h][ci * 3:ci * 3 + 3]
            if len(c0) < 3: continue
            c1 = list(c0)
            for k, a in lst: c1[k] = a
            if str(Seq(c0).translate(table=11)) != str(Seq("".join(c1)).translate(table=11)):
                out.add((h, "prot")); break
    return out

meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "species"])
mg = pd.read_csv(ROOT / "output/mic/mgit/mgit_parsed.csv").merge(meta[["NR", "species"]], on="NR", how="left")
mg = mg[mg.species.notna() & ~mg.NR.isin(EXCL) & ~mg.species.isin(CTRL)].copy()
samples = meta.set_index("NR").loc[sorted(mg.NR.unique()), "PROBENNUMMER"].to_dict()
print("usable MGIT isolates:", len(samples), flush=True)
calls: dict[int, set[tuple[int, str]]] = {}
for nr, s in samples.items():
    t = pd.read_csv(RES / s / "5_typing/kansasii_snippy/snippy_out/snps.tab", sep="\t", dtype=str, usecols=["CHROM", "POS", "TYPE", "REF", "ALT"])
    t["POS"] = t.POS.astype(int)
    calls[nr] = gene_calls(t[["CHROM", "POS", "REF", "TYPE", "ALT"]].rename(columns={"TYPE": "T"})[["CHROM", "POS", "T", "REF", "ALT"]])
    print(" ", s, len(t), "variants ->", len(calls[nr]), "gene calls", flush=True)

feat: dict[tuple[int, str], set[int]] = {}   # (gene_idx, kind) -> {nr: 1}
for nr, cs in calls.items():
    for gk in cs: feat.setdefault(gk, set()).add(nr)
nrs = sorted(samples)
M = pd.DataFrame({f"{genes.locus[g]}|{genes.gene[g]}|{k}": [int(n in v) for n in nrs] for (g, k), v in feat.items()}, index=nrs)
M.to_csv(OUT / "gene_matrix.csv")
print("gene-level features:", M.shape, flush=True)

# ---- phenotypes
mg["R"] = (mg.erg == "R").astype(int)
mg = mg[mg.erg.isin(["S", "I", "R"])]
def tests() -> Iterator[tuple[str, float, str, str, npt.NDArray[Any], npt.NDArray[Any]]]:
    for key, d in mg.groupby(["antibiotic", "concentration_mg_l"]):
        drug, conc = cast("tuple[str, float]", key)
        for sp in ["kansasii", "persicum"]:
            x = d[d.species == sp].drop_duplicates("NR")
            for name, y in [("R_vs_SI", x.R), ("R_vs_S", x.R[x.erg != "I"])]:
                n1, n0 = int(y.sum()), int((y == 0).sum())
                if n1 >= 3 and n0 >= 3: yield drug, conc, sp, name, x.loc[y.index, "NR"].to_numpy(), y.to_numpy()
rows, ph = [], []
KNOWN = re.compile(r"eis|whiB7|WhiB7|aminoglycoside|16S|rrs|rpsL|streptomycin|acetyltransferase", re.I)
for drug, conc, sp, name, ids, y in tests():
    sub = M.loc[ids]
    ph.append(dict(drug=drug, conc=conc, species=sp, contrast=name, n_R=int(y.sum()), n_nonR=int((y == 0).sum())))
    pats: dict[tuple[Any, ...], list[tuple[Any, ...]]] = {}
    for col in sub.columns:
        v = sub[col].to_numpy(); c1 = int(v[y == 1].sum()); c0 = int(v[y == 0].sum())
        a, b, c, d_ = c1, int((y == 1).sum()) - c1, c0, int((y == 0).sum()) - c0
        if a + c < 2 or b + d_ < 2: continue
        _, p = fisher_exact(np.array([[a, b], [c, d_]]))
        pats.setdefault(tuple(v), []).append((col, a, b, c, d_, p))
    allp = np.array([v[0][5] for v in pats.values()])
    o = np.argsort(allp); q = np.empty(len(allp)); n = len(allp)
    q[o] = np.minimum.accumulate((allp[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    for (pat, lst), qq in zip(pats.items(), np.minimum(q, 1)):
        for col, a, b, c, d_, p in lst:
            rows.append(dict(drug=drug, conc=conc, species=sp, contrast=name, feature=col, R_carriers=a, R_noncarriers=b,
                             nonR_carriers=c, nonR_noncarriers=d_, p=p, q_bh_patterns=qq, block_size=len(lst), n_patterns=n))
pd.DataFrame(ph).to_csv(OUT / "mgit_tests_run.csv", index=False)
res = pd.DataFrame(rows).sort_values("p"); res.to_csv(OUT / "mgit_gene_tests.csv", index=False)
print("\ntests run:\n", pd.DataFrame(ph).to_string(index=False))
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 50)
for key, gr in res.groupby(["drug", "conc", "species", "contrast"]):
    print(f"\n=== {key}: patterns tested {gr.n_patterns.iloc[0]}, min q {gr.q_bh_patterns.min():.3f}")
    print(gr.sort_values("p").head(8)[["feature", "R_carriers", "R_noncarriers", "nonR_carriers", "nonR_noncarriers", "p", "q_bh_patterns", "block_size"]].round(4).to_string(index=False))
    kn = gr[gr.feature.str.contains(KNOWN)].sort_values("p").head(6)
    if len(kn): print("  known amikacin/aminoglycoside-related genes:\n", kn[["feature", "R_carriers", "R_noncarriers", "nonR_carriers", "nonR_noncarriers", "p", "block_size"]].round(4).to_string(index=False))
