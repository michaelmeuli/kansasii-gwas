"""Core-SNP distances, NJ tree and sub-lineages within M. kansasii sensu stricto (GTDB 'kansasii', former subtype I).
Uses snippy snps.aligned.fa (full-length alignment to NC_022663.1). Core = sites with an unambiguous base (ACGT) in every genome.
Distances are Hamming counts over variable core sites (no recombination filtering; an approximation, not an ML tree)."""
from __future__ import annotations

import numpy.typing as npt
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas/lineage"; RES = ROOT / "runs/mkan329/assembly/results"
EXCL = {35, 113, 106, 10, 21}
m = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["NR", "PROBENNUMMER", "species"])
m = m[(m.species == "kansasii") & ~m.NR.isin(EXCL)].sort_values("NR")
names = m.PROBENNUMMER.tolist(); print("kansasii genomes:", len(names), flush=True)

def read_chr(path: Path) -> npt.NDArray[np.uint8]:
    seq: list[str] = []; on = False
    for l in open(path):
        if l.startswith(">"):
            if on: break
            on = l[1:].split()[0] == "NC_022663.1"
        elif on: seq.append(l.strip().upper())
    return np.frombuffer("".join(seq).encode(), dtype=np.uint8)

A = np.stack([read_chr(RES / n / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa") for n in names])
print("alignment", A.shape, flush=True)
valid = np.isin(A, np.frombuffer(b"ACGT", dtype=np.uint8))
core = valid.all(axis=0)
print(f"core sites (unambiguous in all {len(names)}): {core.sum()} of {A.shape[1]} ({100*core.mean():.1f}%)", flush=True)
C = A[:, core]
var = (C != C[0]).any(axis=0)
V = C[:, var]; print("variable core sites:", V.shape[1], flush=True)
pos = np.flatnonzero(core)[var] + 1
pd.DataFrame({"pos": pos}).to_csv(OUT / "core_snp_positions.csv", index=False)
with open(OUT / "kansasii_core_snps.fasta", "w") as f:
    for n, row in zip(names, V): f.write(f">{n}\n{row.tobytes().decode()}\n")
D = np.array([(V != V[i]).sum(axis=1) for i in range(len(names))])
pd.DataFrame(D, index=names, columns=names).to_csv(OUT / "core_snp_distance_matrix.csv")
iu = np.triu_indices(len(names), 1)
print("pairwise core-SNP distance: min %d, median %d, max %d" % (D[iu].min(), np.median(D[iu]), D[iu].max()))
print("nearest-neighbour distance quantiles:", np.percentile(np.sort(D + np.eye(len(D)) * 10**9, axis=1)[:, 0], [0, 25, 50, 75, 100]).astype(int).tolist(), flush=True)

dm = DistanceMatrix(names, [[int(D[i, j]) for j in range(i + 1)] for i in range(len(names))])
tree = DistanceTreeConstructor().nj(dm); Phylo.write(tree, str(OUT / "kansasii_nj.nwk"), "newick")  # type: ignore[attr-defined,no-untyped-call]
Z = linkage(squareform(D, checks=False), "average"); np.save(OUT / "linkage.npy", Z)
print("\ncluster counts and sizes at cut heights (average linkage, core SNPs):")
for h in (500, 1000, 2000, 4000, 6000, 8000, 10000, 15000):
    cl = fcluster(Z, h, "distance"); s = sorted(pd.Series(cl).value_counts().tolist(), reverse=True)
    print(f"  h={h:6d}: {len(set(cl)):3d} clusters, sizes {s[:10]}", flush=True)
# gaps in merge heights: large jumps suggest natural lineage boundaries
hts = Z[:, 2]; gaps = np.diff(hts); top = np.argsort(gaps)[-6:][::-1]
print("largest jumps between successive merge heights (from -> to, clusters left after the lower merge):")
for i in top: print(f"  {hts[i]:.0f} -> {hts[i+1]:.0f}  (clusters at that level: {len(names) - 1 - i})")
