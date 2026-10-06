"""rRNA hotspot check in E. coli numbering: rrl A2058/A2059 (macrolides), rrs A1408 (aminoglycosides)."""
from __future__ import annotations

from typing import Any
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Seq import Seq

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas"
ec = "".join(l.strip() for l in open(OUT / "ref/Ecoli_rrnA_region.fa") if not l.startswith(">"))
s16 = ec.find("AAATTGAAGAGTTTGATCATGGCTCAG"); s23 = ec.find("GGTTAAGCGACTAAGCGTACACGGTGGATGCC")
ec16, ec23 = ec[s16:s16 + 1542], ec[s23 + 1:s23 + 1 + 2904]  # +1: motif hit sits one base 5-prime of the numbered 23S start (verified by 2057-2063 == GAAAGAC)
assert ec23[2056:2063] == "GAAAGAC", ec23[2056:2063]
print("E. coli 16S start/23S start in region:", s16, s23, "| 23S 2057-2063:", ec23[2056:2063], "| 16S 1405-1411:", ec16[1404:1411])
genome = {r.id: str(r.seq) for r in SeqIO.parse(ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.fna", "fasta")}["contig_1"]
G = pd.read_csv(OUT / "candidate_gene_coords.csv"); G = G[(G.chrom == "NC_022663.1") & G.gene.isin(["rrl", "rrs"])].set_index("gene")
al = Align.PairwiseAligner(); al.mode = "global"; al.match_score, al.mismatch_score, al.open_gap_score, al.extend_gap_score = 2, -1, -4, -0.5
feat = pd.read_csv(OUT / "candidate_features_per_sample.csv"); feat["features"] = feat.features.apply(eval)
sp = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "species"]).set_index("PROBENNUMMER").species
HOT = {"rrl": (ec23, {2057: "G", 2058: "A", 2059: "A", 2503: "A", 2611: "C"}), "rrs": (ec16, {1408: "A", 1409: "A", 1491: "G"})}
for gene, (ecseq, sites) in HOT.items():
    r: Any = G.loc[gene]; kan = genome[r.start - 1:r.end]; kan = kan if r.strand == "+" else str(Seq(kan).reverse_complement())
    aln = al.align(ecseq, kan)[0]
    m: dict[int, int] = {}
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for i in range(a1 - a0):
            m[a0 + i + 1] = b0 + i + 1
    ident = sum(ecseq[a - 1] == kan[b - 1] for a, b in m.items()) / len(m)
    print(f"\n{gene}: kansasii len {len(kan)} vs E. coli {len(ecseq)}; aligned identity {ident:.2f}")
    for ep, wt in sites.items():
        kp = m.get(ep)
        if kp is None:
            print(f"  E.coli {ep}: no aligned base"); continue
        gpos = r.start + kp - 1 if r.strand == "+" else r.end - kp + 1
        hits = [(s, f) for s, fs in zip(feat["sample"], feat["features"]) for f in fs if f.startswith(f"{gene}|") and f.split(":")[1][1:-1] == str(gpos) or (f"indel@{gpos}" in f)]
        print(f"  E.coli {ep} (expected {wt}, actual {ecseq[ep-1]}) -> kansasii {gene} pos {kp} base {kan[kp-1]} (genome {gpos}); isolates with a call there: {len(hits)}")
        for s, f in hits[:8]:
            print("     ", s, sp.get(s), f)
