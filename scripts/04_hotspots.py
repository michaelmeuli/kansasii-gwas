"""Hotspot check. Map kansasii (Bakta) residues/bases to M. tuberculosis (proteins) or E. coli (rRNA) numbering by
alignment, verify the expected wild-type residue at each hotspot, then tabulate each isolate's state and MIC."""
from __future__ import annotations

import os
from typing import Any
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Seq import Seq

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
OUT = ROOT / "output/gwas"
REF = OUT / "ref"
REFFNA = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.fna"
genome = {r.id: str(r.seq) for r in SeqIO.parse(REFFNA, "fasta")}["contig_1"]
G = pd.read_csv(OUT / "candidate_gene_coords.csv")
G = G[G.chrom == "NC_022663.1"]

# Mtb numbering; wild-type residue in Mtb checked against the downloaded sequence, not assumed
HOT = {"gyrA": {90: "A", 91: "S", 94: "D", 88: "G", 89: "D"},
       "gyrB": {461: "D", 464: "G", 499: "N", 501: "E"},
       "rpoB": {435: "D", 445: "H", 450: "S", 452: "L", 432: "Q"},
       "embB": {306: "M", 406: "G", 497: "Q"},
       "katG": {315: "S"}, "rpsL": {43: "K", 88: "K"}, "rplC": {154: "C"}}
# clinical (Mtb) numbering vs UniProt sequence index: rpoB is +6 (clinical S450 = UniProt S456)
UNIPROT_OFFSET = {"rpoB": 6}
DRUG = {"gyrA": ["Ciprofloxacin", "Moxifloxacin"], "gyrB": ["Ciprofloxacin", "Moxifloxacin"],
        "rpoB": ["Rifampicin", "Rifabutin"], "embB": ["Ethambutol"], "katG": ["Isoniazid"],
        "rpsL": ["Streptomycin"], "rplC": ["Linezolid"]}

def kan_prot(name: str) -> tuple[str, Any]:
    r = G[G.gene == name].iloc[0]
    s = genome[r.start - 1:r.end]
    s = s if r.strand == "+" else str(Seq(s).reverse_complement())
    return str(Seq(s).translate(table=11)), r

al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
al.open_gap_score, al.extend_gap_score = -10, -0.5

def pos_map(a: str, b: str) -> dict[int, int]:
    """map positions of a -> positions of b (1-based) from the best alignment"""
    aln = al.align(a, b)[0]
    m: dict[int, int] = {}
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i in range(e1 - s1):
            m[s1 + i + 1] = s2 + i + 1
    return m

feat = pd.read_csv(OUT / "candidate_features_per_sample.csv")
feat["features"] = feat.features.apply(eval)
pheno = pd.read_csv(OUT / "geno_pheno_table_long.csv")
sp = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "species"]).set_index("PROBENNUMMER").species
check: list[dict[str, Any]] = []
rows: list[dict[str, Any]] = []
for gene, hs in HOT.items():
    mtb = "".join(l.strip() for l in open(REF / f"Mtb_{gene}.faa", encoding="utf-8") if not l.startswith(">"))
    kan, r = kan_prot(gene)
    m = pos_map(mtb, kan)  # Mtb -> kansasii
    ident = sum(mtb[a - 1] == kan[b - 1] for a, b in m.items()) / len(m)
    off = UNIPROT_OFFSET.get(gene, 0)
    for mp, wt in hs.items():
        kp = m.get(mp + off)
        check.append(dict(gene=gene, mtb_pos=mp, mtb_wt_expected=wt, mtb_wt_actual=mtb[mp + off - 1],
                          kan_pos=kp, kan_ref_residue=kan[kp - 1] if kp else None, aln_identity=round(ident, 2)))
        if kp is None:
            continue
        for s, fs in zip(feat["sample"], feat["features"]):
            alts = [f.split(":")[1] for f in fs if f.startswith(f"{gene}|{r.locus}:") and re.fullmatch(rf"[A-Z*]{kp}[A-Z*]", f.split(":")[1])]
            state = alts[0] if alts else f"{kan[kp - 1]}{kp} (ref)"
            rows.append(dict(sample=s, gene=gene, mtb_pos=mp, kan_pos=kp, kan_ref=kan[kp - 1], state=state, variant=bool(alts)))
chk = pd.DataFrame(check)
chk.to_csv(OUT / "hotspot_numbering_check.csv", index=False)
pd.set_option("display.width", 200)
print(chk.to_string(index=False))
ho = pd.DataFrame(rows)
ho.to_csv(OUT / "hotspot_genotypes_per_isolate.csv", index=False)
print("\nvariant calls at hotspots (all 131 genomes):")
v = ho[ho.variant].merge(sp.rename("gtdb"), left_on="sample", right_index=True, how="left")
print(v.groupby(["gene", "mtb_pos", "kan_ref", "kan_pos", "state", "gtdb"], dropna=False).size().to_string() or "none")
