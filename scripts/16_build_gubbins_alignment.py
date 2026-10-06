"""Full-length chromosome alignment (NC_022663.1) of the 82 non-excluded kansasii genomes + the reference, from snippy snps.aligned.fa
('-' = no coverage, 'N' = low coverage; Gubbins treats both as missing). Writes output/gwas/gubbins/kansasii_chr_aln.fasta."""
from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas/gubbins"; RES = ROOT / "runs/mkan329/assembly/results"
EXCL = {35, 113, 106, 10, 21}
m = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "species"])
m = m[(m.species == "kansasii") & ~m.NR.isin(EXCL)].sort_values("NR")

def read_chr(path: Path) -> str:
    seq: list[str] = []; on = False
    for l in open(path):
        if l.startswith(">"):
            if on: break
            on = l[1:].split()[0] == "NC_022663.1"
        elif on: seq.append(l.strip().upper())
    return "".join(seq)

ref = str(next(SeqIO.parse(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa", "fasta")).seq).upper()
rows: list[tuple[str, float]] = []
with open(OUT / "kansasii_chr_aln.fasta", "w") as f:
    f.write(f">Reference_ATCC12478\n{ref}\n")
    for n in m.PROBENNUMMER:
        s = read_chr(RES / n / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa")
        assert len(s) == len(ref), (n, len(s), len(ref))
        a = np.frombuffer(s.encode(), dtype=np.uint8)
        miss = 1 - np.isin(a, np.frombuffer(b"ACGT", dtype=np.uint8)).mean()
        rows.append((n, round(100 * miss, 2))); f.write(f">{n}\n{s}\n")
d = pd.DataFrame(rows, columns=["sample", "pct_missing"]); d.to_csv(OUT / "missing_per_isolate.csv", index=False)
print("sequences:", len(d) + 1, "| length", len(ref))
print(f"missing per isolate (%): min {d.pct_missing.min():.2f}, median {d.pct_missing.median():.2f}, max {d.pct_missing.max():.2f}; isolates over 25%: {(d.pct_missing > 25).sum()}")
