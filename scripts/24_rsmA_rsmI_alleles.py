"""Compare RsmA (FCNCLM_01703) and RsmI (FCNCLM_01716) alleles: kansasii reference vs the persicum-like alleles, across all 82 cohort kansasii genomes
(own Bakta proteins), related to the segment class from script 22/23, and checked in public persicum and kansasii genomes (tblastn).
Conservation is judged against the M. tuberculosis orthologs (KsgA/Rv1010, RsmI/Rv2372c): a position is 'conserved' when kansasii reference and Mtb agree."""
from __future__ import annotations

from typing import Any
import subprocess
import pandas as pd
from collections import Counter
from pathlib import Path
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); G = ROOT / "output/gwas"; OUT = G / "boundaries"; RES = ROOT / "runs/mkan329/assembly/results"
REFFAA = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation/GCF_000157895.3_query.faa"
PUB = ROOT / "data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii"
ref = {r.id: str(r.seq) for r in SeqIO.parse(REFFAA, "fasta")}
GENES = {"RsmA": ("FCNCLM_01703", "ksgA"), "RsmI": ("FCNCLM_01716", "rsmI")}
B62 = substitution_matrices.load("BLOSUM62")
al = Align.PairwiseAligner(); al.mode = "global"; al.substitution_matrix = B62; al.open_gap_score, al.extend_gap_score = -10, -0.5
def diffs(a: str, b: str) -> list[tuple[int, str, str]]:
    """substitutions/indels of b relative to a, using a's numbering"""
    aln = al.align(a, b)[0]; out: list[tuple[int, str, str]] = []; (ra, rb) = aln.aligned
    for (a0, a1), (b0, b1) in zip(ra, rb):
        out += [(a0 + i + 1, a[a0 + i], b[b0 + i]) for i in range(a1 - a0) if a[a0 + i] != b[b0 + i]]
    covered = sum(a1 - a0 for a0, a1 in ra)
    if covered < len(a) or sum(b1 - b0 for b0, b1 in rb) < len(b): out.append((0, "indel", f"aligned {covered}/{len(a)} vs {sum(b1-b0 for b0,b1 in rb)}/{len(b)}"))
    return out
seg = pd.read_csv(OUT / "segment_classes_v2.csv", index_col=0)
rows: list[dict[str, Any]] = []
for s in seg.index:
    prots = {r.id: str(r.seq) for r in SeqIO.parse(RES / s / f"2_annotation/{s}.faa", "fasta")}
    row = {"sample": s, "cls": seg.cls[s]}
    for name, (lt, _) in GENES.items():
        cand = [(al.score(ref[lt], p), p) for p in prots.values() if 0.8 * len(ref[lt]) < len(p) < 1.2 * len(ref[lt])]
        row[name] = max(cand)[1] if cand else None
    rows.append(row)
A = pd.DataFrame(rows).set_index("sample")
for name, (lt, mg) in GENES.items():
    mtb = "".join(l.strip() for l in open(G / f"ref/Mtb_{mg}.faa") if not l.startswith(">"))
    mtbmap: dict[int, str] = {}
    aln = al.align(ref[lt], mtb)[0]
    for (a0, a1), (b0, b1) in zip(*aln.aligned):
        for i in range(a1 - a0): mtbmap[a0 + i + 1] = mtb[b0 + i]
    print(f"\n================ {name} ({lt}, {len(ref[lt])} aa) ================")
    A[f"{name}_allele"] = [("ref" if p == ref[lt] else p) if p else None for p in A[name]]
    alleles = Counter(p for p in A[name] if p)
    ids = {p: ("REF-like" if p == ref[lt] else f"alleleP{k}") for k, p in enumerate([p for p, _ in alleles.most_common() if p != ref[lt]], 1)}
    ids[ref[lt]] = "REF-like"
    A[f"{name}_id"] = A[name].map(ids)
    print("allele counts across 82 kansasii genomes:", dict(Counter(A[f"{name}_id"])))
    print("allele by segment class:\n", pd.crosstab(A[f"{name}_id"], A.cls).to_string())
    for p, a in ids.items():
        if a == "REF-like": continue
        d = diffs(ref[lt], p)
        conserved = [(pos, x, y) for pos, x, y in d if pos and mtbmap.get(pos) == x]
        radical = [(pos, x, y) for pos, x, y in d if pos and B62[x][y] < 0]
        print(f"  {a} (n={alleles[p]}): {len(d)} differences vs kansasii reference; at positions where reference = Mtb: {len(conserved)}; BLOSUM62-negative: {len(radical)}; length {len(p)}")
        print("    ", ", ".join(f"{x}{pos}{y}" + ("*" if mtbmap.get(pos) == x else "") for pos, x, y in d if pos) + ("  [* = reference residue identical in Mtb]"))
        print(f"     identity to reference {100*(1-len([1 for pos,_,_ in d if pos])/len(ref[lt])):.1f}%")
        # does the persicum-like allele agree with Mtb where the reference does not?
        toMtb = [(pos, x, y) for pos, x, y in d if pos and mtbmap.get(pos) == y]
        print(f"     substitutions that match the Mtb residue: {len(toMtb)}")
    A[f"{name}_id"].to_csv(OUT / f"{name}_allele_per_isolate.csv")
    # public genomes via tblastn with the reference allele and the main persicum-like allele
    main = next((p for p, a in ids.items() if a == "alleleP1"), None)
    if main:
        for sp in ("persicum", "kansasii"):
            fa = OUT / f"pub_{sp}.fa"
            if not fa.exists():
                with open(fa, "w") as f:
                    for gf in sorted((PUB / sp).glob("*.fna")):
                        for r in SeqIO.parse(gf, "fasta"): f.write(f">{gf.stem}|{r.id}\n{str(r.seq)}\n")
                subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(OUT / f"pub_{sp}")], check=True, stdout=subprocess.DEVNULL)
            res: dict[str, dict[str, tuple[float, float, int]]] = {}
            for tag, q in (("ref", ref[lt]), ("alleleP1", main)):
                (OUT / "q.faa").write_text(f">q\n{q}\n")
                out = subprocess.run(["tblastn", "-query", str(OUT / "q.faa"), "-db", str(OUT / f"pub_{sp}"), "-outfmt", "6 sseqid pident length bitscore", "-evalue", "1e-30", "-max_target_seqs", "200"],
                                     check=True, capture_output=True, text=True).stdout.strip().split("\n")
                best: dict[str, tuple[float, float, int]] = {}
                for l in out:
                    if not l: continue
                    sid, pid, ln, bs = l.split("\t"); g = sid.split("|")[0]
                    if g not in best or float(bs) > best[g][1]: best[g] = (float(pid), float(bs), int(ln))
                res[tag] = best
            gen = sorted(set(res["ref"]) | set(res["alleleP1"]))
            assign: Counter[str] = Counter()
            for g in gen:
                sr, sp1 = res["ref"].get(g, (0, 0, 0)), res["alleleP1"].get(g, (0, 0, 0))
                assign["alleleP1-like" if sp1[1] > sr[1] else "REF-like" if sr[1] > sp1[1] else "equal"] += 1
            print(f"  public {sp} genomes ({len(gen)}): best match = {dict(assign)}")
