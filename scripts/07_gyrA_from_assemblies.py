"""gyrA QRDR residues (Mtb numbering) read from each isolate's own Bakta protein annotation (no reference mapping)."""
from __future__ import annotations

import os
from typing import Any
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); OUT = ROOT / "output/gwas"
RES = ROOT / "runs/mkan329/assembly/results"
mtb = "".join(l.strip() for l in open(OUT / "ref/Mtb_gyrA.faa", encoding="utf-8") if not l.startswith(">"))
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5
sites = [74, 88, 89, 90, 91, 94]   # QRDR residues, Mtb numbering (wt A90, D94 verified earlier)
sp = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "species", "NR"]).set_index("PROBENNUMMER")
rows: list[dict[str, Any]] = []
for d in sorted(RES.glob("Mkan329-*")):
    faa = d / f"2_annotation/{d.name}.faa"
    if not faa.exists():
        continue
    cands = [r for r in SeqIO.parse(faa, "fasta") if "gyrase subunit A" in r.description or "(ATP-hydrolyzing) subunit A" in r.description]
    if not cands:
        rows.append(dict(sample=d.name, status="no gyrA annotated")); continue
    best = max(cands, key=lambda r: al.score(mtb[:300], str(r.seq)[:400]))
    p = str(best.seq)
    aln = al.align(mtb, p)[0]
    m: dict[int, int] = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for i in range(a1 - a0):
            m[a0 + i + 1] = b0 + i + 1
    row = dict(sample=d.name, status="ok", copies=len(cands), prot_len=len(p))
    for s in sites:
        row[f"r{s}"] = (mtb[s - 1] + "->" + (p[m[s] - 1] if s in m else "-"))
    rows.append(row)
df = pd.DataFrame(rows).merge(sp, left_on="sample", right_index=True, how="left")
df.to_csv(OUT / "gyrA_QRDR_from_assemblies.csv", index=False)
print("annotated gyrA in", (df.status == "ok").sum(), "of", len(df), "genomes")
cols = [f"r{s}" for s in sites]
sub = df[df.status == "ok"]
print("\nresidue patterns (Mtb wt -> isolate):")
print(sub.groupby(["species"] + cols, dropna=False).size().to_string())
print("\nnot annotated:", df[df.status != "ok"]["sample"].tolist())
