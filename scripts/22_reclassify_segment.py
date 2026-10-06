"""Reclassify all 82 kansasii genomes on the BLAST-derived extent of the persicum-derived segment (NC_022663.1:1,929,793-1,994,501; script 21).
Diagnostic sites = positions in the segment where all contig-complete carriers share a non-reference base and both reference-like controls equal the reference.
Left part 1,929,793-1,950,330 (persicum-like also in the ~211 haplotype), right part 1,950,331-1,994,501 (only in the full haplotype)."""
from __future__ import annotations

from typing import Any
import numpy.typing as npt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import fisher_exact

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"; OUT = G / "boundaries"; RES = ROOT / "runs/mkan329/assembly/results"
LO, HI, MID = 1_929_793, 1_994_501, 1_950_330
CARR = ["Mkan329-001", "Mkan329-116", "Mkan329-003", "Mkan329-013", "Mkan329-020", "Mkan329-066", "Mkan329-101"]; CTRL = ["Mkan329-011", "Mkan329-014"]
def read_chr(s: str) -> npt.NDArray[np.uint8]:
    seq: list[str] = []; on = False
    for l in open(RES / s / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa"):
        if l.startswith(">"):
            if on: break
            on = l[1:].split()[0] == "NC_022663.1"
        elif on: seq.append(l.strip().upper())
    return np.frombuffer("".join(seq).encode(), dtype=np.uint8)[LO - 1:HI]
lin = pd.read_csv(G / "lineage/kansasii_region_haplotype.csv", index_col=0)
A = {s: read_chr(s) for s in lin.index}
ref = np.frombuffer(next(l for l in open(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa") if False) if False else b"", dtype=np.uint8) if False else None
from Bio import SeqIO
refseq = np.frombuffer(str(next(SeqIO.parse(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa", "fasta")).seq).upper().encode(), dtype=np.uint8)[LO - 1:HI]
acgt = np.frombuffer(b"ACGT", dtype=np.uint8)
carr = np.stack([A[s] for s in CARR]); ctrl = np.stack([A[s] for s in CTRL])
ok = np.isin(carr, acgt).all(0) & np.isin(ctrl, acgt).all(0)
share = (carr == carr[0]).all(0) & (carr[0] != refseq) & (ctrl == refseq).all(0) & ok
pos = np.flatnonzero(share) + LO
print(f"diagnostic sites: {share.sum()} over {HI - LO + 1} bp ({share.sum() / ((HI - LO + 1) / 1000):.1f}/kb); left {int((pos <= MID).sum())}, right {int((pos > MID).sum())}")
idx = np.flatnonzero(share); alt = carr[0][idx]; refb = refseq[idx]; left = pos <= MID
rows: list[dict[str, Any]] = []
for s in lin.index:
    b = A[s][idx]; called = np.isin(b, acgt)
    def part(m: npt.NDArray[np.bool_]) -> tuple[float, float]:
        c = called & m; n = int(c.sum())
        return (float((b[c] == alt[c]).mean()) if n else np.nan), n / max(1, int(m.sum()))
    fl, cl = part(left); fr, cr = part(~left)
    rows.append(dict(sample=s, frac_persicum_left=fl, called_left=cl, frac_persicum_right=fr, called_right=cr))
d = pd.DataFrame(rows).set_index("sample")
def cls(r: Any) -> str:
    L = r.frac_persicum_left >= 0.9 and r.called_left >= 0.5; R = r.frac_persicum_right >= 0.9 and r.called_right >= 0.5
    Ln = r.frac_persicum_left <= 0.1 and r.called_left >= 0.5; Rn = r.frac_persicum_right <= 0.1 and r.called_right >= 0.5
    if L and R: return "full"
    if L and Rn: return "left-only"
    if Ln and R: return "right-only"
    if Ln and Rn: return "none"
    return "ambiguous"
d["segment_class"] = d.apply(cls, axis=1)
d = d.join(lin[["NR", "clade", "hap", "n_var_region"]]); d.to_csv(OUT / "segment_classes.csv")
print("\nnew classification vs the old region-based haplotype class:\n", pd.crosstab(d.hap.astype(str), d.segment_class).to_string())
print("\nambiguous isolates:\n", d[d.segment_class == "ambiguous"][["frac_persicum_left", "called_left", "frac_persicum_right", "called_right", "hap"]].round(2).to_string())
# genes inside the segment
grows: list[tuple[str | None, int, int, str, str]] = []
for l in open(ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.gff3"):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[0] != "contig_1" or f[2] != "CDS": continue
    a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    grows.append((a.get("locus_tag"), int(f[3]), int(f[4]), a.get("gene", ""), a.get("product", "")))
gn = pd.DataFrame(grows, columns=["locus", "start", "end", "gene", "product"])
inside = gn[(gn.end >= LO) & (gn.start <= HI)]; inside = inside.assign(part=np.where(inside.start <= MID, "left", "right"))
inside.to_csv(OUT / "segment_genes.csv", index=False)
print(f"\ngenes overlapping the segment: {len(inside)} ({inside.locus.iloc[0]} .. {inside.locus.iloc[-1]}); left part {int((inside.part == 'left').sum())}, right part {int((inside.part == 'right').sum())}")
print("genes spanning the left boundary:", inside[(inside.start < LO + 1500) | (inside.start <= LO)].tail(2)[["locus", "gene", "start", "end"]].values.tolist(),
      "| right boundary:", inside[inside.end >= HI - 1500][["locus", "gene", "start", "end"]].values.tolist())
