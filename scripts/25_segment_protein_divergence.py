"""Is the divergence of RsmA/RsmI persicum-like alleles unusual? Protein identity between the kansasii reference and a full-segment carrier (Mkan329-001)
for all 59 genes of the segment, plus a check of the SAM-binding (GxGxG) motif region in each allele."""
from __future__ import annotations

from typing import Any
import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"; RES = ROOT / "runs/mkan329/assembly/results"
REFFAA = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.faa"
ref = {r.id: str(r.seq) for r in SeqIO.parse(REFFAA, "fasta")}
genes = pd.read_csv(G / "boundaries/segment_genes.csv")
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62"); al.open_gap_score, al.extend_gap_score = -10, -0.5
iso = {r.id: str(r.seq) for r in SeqIO.parse(RES / "Mkan329-001/2_annotation/Mkan329-001.faa", "fasta")}
def ident(a: str, b: str) -> float:
    aln = al.align(a, b)[0]; m = sum(a[a0 + i] == b[b0 + i] for (a0, a1), (b0, b1) in zip(*aln.aligned) for i in range(a1 - a0))
    return m / max(len(a), len(b))
rows: list[dict[str, Any]] = []
r: Any
for r in genes.itertuples():
    a = ref.get(r.locus)
    if not a: continue
    cands = [(al.score(a, p), p) for p in iso.values() if 0.5 * len(a) < len(p) < 1.5 * len(a)]
    if not cands: continue
    rows.append(dict(locus=r.locus, gene=r.gene if isinstance(r.gene, str) else "", length=len(a), identity=ident(a, max(cands)[1])))
d = pd.DataFrame(rows); d["pct"] = (100 * d.identity).round(1)
print(f"{len(d)} segment proteins compared (kansasii reference vs Mkan329-001): median identity {d.pct.median():.1f}%, quartiles {d.pct.quantile(.25):.1f}-{d.pct.quantile(.75):.1f}%, min {d.pct.min():.1f}%")
for lt in ("FCNCLM_01703", "FCNCLM_01716"):
    r = d[d.locus == lt].iloc[0]; rank = int((d.identity < r.identity).sum()) + 1
    print(f"  {lt} ({r.gene}): {r.pct}% identity; {rank} of {len(d)} when ranked from most to least divergent; {100*(d.identity > r.identity).mean():.0f}% of segment proteins are more conserved")
d.to_csv(G / "boundaries/segment_protein_identity.csv", index=False)
print("most divergent segment proteins:\n", d.sort_values("identity").head(6)[["locus", "gene", "length", "pct"]].to_string(index=False))
# NOTE: a core-motif check (GxGxG-type SAM-binding pattern) was attempted but the pattern matched neither protein, so it is not reported.
