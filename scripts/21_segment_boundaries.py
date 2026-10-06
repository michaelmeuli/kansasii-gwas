"""Redefine the extent of the persicum-derived segment around FCNCLM_01713-01736 (NC_022663.1:1,944,065-1,970,050).
For each isolate, extract 1.80-2.10 Mb (reference coordinates) from its OWN assembly, tile it into 1 kb windows (500 bp step) and BLAST each window against
(a) the reference chromosome and (b) the 12 public M. persicum genomes. A window is persicum-like when identity to persicum >= 99.3% and at least 2 points
above identity to the reference; kansasii-like when reference identity >= 98.5% and persicum identity at least 2 points lower. Window positions are
expressed in reference coordinates through the window's best reference hit. Controls: isolates with the reference-like and ~211-variant haplotypes."""
from __future__ import annotations

from collections.abc import Sequence
from typing import Any
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas/boundaries"; OUT.mkdir(exist_ok=True)
RES = ROOT / "runs/mkan329/assembly/results"; PUB = ROOT / "data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii"
SCAN_LO, SCAN_HI = 1_800_000, 2_100_000; W, STEP = 1000, 500
ref = str(next(SeqIO.parse(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa", "fasta")).seq).upper()
(OUT / "ref.fa").write_text(">NC_022663.1\n" + ref + "\n")
with open(OUT / "persicum.fa", "w") as f:
    for fa in sorted((PUB / "persicum").glob("*.fna")):
        for r in SeqIO.parse(fa, "fasta"): f.write(f">{fa.stem}|{r.id}\n{str(r.seq)}\n")
for fname, db in (("ref.fa", "refdb"), ("persicum.fa", "persdb")):
    if not (OUT / f"{db}.nsq").exists() and not (OUT / f"{db}.nin").exists():
        subprocess.run(["makeblastdb", "-in", str(OUT / fname), "-dbtype", "nucl", "-out", str(OUT / db)], check=True, stdout=subprocess.DEVNULL)
hap = pd.read_csv(ROOT / "output/gwas/lineage/kansasii_region_haplotype.csv", index_col=0)
H880 = ["Mkan329-001", "Mkan329-003", "Mkan329-012", "Mkan329-013", "Mkan329-020", "Mkan329-022", "Mkan329-101", "Mkan329-116", "Mkan329-066", "Mkan329-072"]
CTRL = ["Mkan329-011", "Mkan329-014", "Mkan329-005", "Mkan329-015"]
print("hap classes:", {s: str(hap.hap.get(s)) for s in H880 + CTRL})

def blast(q: Path, db: str, fields: str, extra: Sequence[str] = ()) -> list[list[str]]:
    r = subprocess.run(["blastn", "-query", str(q), "-db", str(OUT / db), "-outfmt", f"6 {fields}", "-evalue", "1e-10", "-num_threads", "4", *extra],
                       check=True, capture_output=True, text=True).stdout
    return [l.split("\t") for l in r.strip().split("\n")] if r.strip() else []

def windows_for(sample: str) -> tuple[str, str, int, tuple[int, int, bool]]:
    fa = RES / sample / f"1_unicycler/{sample}.fasta"; tmp = OUT / f"tmp_{sample}"; tmp.mkdir(exist_ok=True)
    (tmp / "q.fa").write_text(">scan\n" + ref[SCAN_LO - 1:SCAN_HI] + "\n")
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(tmp / "db")], check=True, stdout=subprocess.DEVNULL)
    h = pd.DataFrame(blast_local(tmp / "q.fa", tmp / "db"), columns=["c", "qs", "qe", "ss", "se", "len"]).astype({"qs": int, "qe": int, "ss": int, "se": int, "len": int})
    best = str(h.groupby("c").len.sum().idxmax()); d = h[h.c == best]
    contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(fa, "fasta")}; seq = contigs[best]
    minus = bool(d.sort_values("len", ascending=False).iloc[0].se < d.sort_values("len", ascending=False).iloc[0].ss)
    a, b = int(min(d.ss.min(), d.se.min())), int(max(d.ss.max(), d.se.max()))
    seg = seq[a - 1:b]; seg = str(Seq(seg).reverse_complement()) if minus else seg
    return seg, best, len(seq), (a, b, minus)

def blast_local(q: Path, db: Path) -> list[list[str]]:
    r = subprocess.run(["blastn", "-query", str(q), "-db", str(db), "-evalue", "1e-10", "-outfmt", "6 sseqid qstart qend sstart send length"],
                       check=True, capture_output=True, text=True).stdout
    return [l.split("\t") for l in r.strip().split("\n")]

rows: list[dict[str, Any]] = []
for s in H880 + CTRL:
    seg, contig, clen, span = windows_for(s)
    qf = OUT / f"win_{s}.fa"
    with open(qf, "w") as f:
        starts = list(range(0, max(1, len(seg) - W + 1), STEP))
        for i in starts: f.write(f">w{i}\n{seg[i:i + W]}\n")
    def best(db: str, with_ref: bool) -> dict[str, tuple[float, float, int, int, int]]:
        out: dict[str, tuple[float, float, int, int, int]] = {}
        for r in blast(qf, db, "qseqid pident length sstart send bitscore"):
            w = r[0]; bs = float(r[5])
            if w not in out or bs > out[w][0]: out[w] = (bs, float(r[1]), int(r[2]), int(r[3]), int(r[4]))
        return out
    R = best("refdb", True); P = best("persdb", False)
    for i in starts:
        w = f"w{i}"; r = R.get(w); p = P.get(w)
        rows.append(dict(sample=s, contig=contig, contig_len=clen, win_start=i, ref_pid=r[1] if r and r[2] >= 0.8 * W else np.nan,
                         ref_pos=(min(r[3], r[4]) if r and r[2] >= 0.8 * W else np.nan), pers_pid=p[1] if p and p[2] >= 0.8 * W else np.nan))
    print(f"{s}: contig {contig} ({clen} bp), extracted {len(seg)} bp, windows {len(starts)}", flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT / "window_identities.csv", index=False)
def call(r: Any) -> str:
    if np.isnan(r.pers_pid) and np.isnan(r.ref_pid): return "unaligned"
    pp, rp = (0 if np.isnan(r.pers_pid) else r.pers_pid), (0 if np.isnan(r.ref_pid) else r.ref_pid)
    if pp >= 99.3 and pp - rp >= 2: return "persicum-like"
    if rp >= 98.5 and rp - pp >= 2: return "kansasii-like"
    return "ambiguous"
df["call"] = df.apply(call, axis=1); df.to_csv(OUT / "window_calls.csv", index=False)
print("\nwindow calls per isolate:\n", pd.crosstab(df["sample"], df.call).to_string())
