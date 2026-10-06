"""Variants carried by an MIC-outlier isolate but by no other genome of the same GTDB species (and, as a looser
tier, by at most one other). Annotated against the Bakta annotation of the reference; multi-base calls split per
base and merged per codon. Background = all same-species genomes (not only the ones with MIC), minus contaminated cultures."""
from __future__ import annotations

from typing import Any
import re
import pandas as pd
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas"
RES = ROOT / "runs/mkan329/assembly/results"
REFDIR = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation"
CONTIG = {"contig_1": "NC_022663.1", "contig_2": "NC_022654.1"}
EXCL = {35, 113, 106, 10, 21}
FOCAL = {"Mkan329-054": "persicum", "Mkan329-055": "persicum", "Mkan329-047": "kansasii", "Mkan329-016": "pseudokansasii",
         "Mkan329-050": "kansasii", "Mkan329-009": "persicum", "Mkan329-008": "innocens"}
genome = {CONTIG[r.id]: str(r.seq) for r in SeqIO.parse(REFDIR / "GCF_000157895.3_query.fna", "fasta")}

feats: list[tuple[str, int, int, str, str, str, str, str]] = []
for l in open(REFDIR / "GCF_000157895.3_query.gff3"):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("CDS", "rRNA", "tRNA", "ncRNA"): continue
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    feats.append((CONTIG[f[0]], int(f[3]), int(f[4]), f[6], f[2], attrs.get("gene", ""), attrs.get("product", ""), attrs.get("locus_tag", "")))
byc: defaultdict[str, list[tuple[str, int, int, str, str, str, str, str]]] = defaultdict(list)
for x in feats: byc[x[0]].append(x)
for c in byc: byc[c].sort(key=lambda x: x[1])

def hit(chrom: str, pos: int) -> tuple[str, int, int, str, str, str, str, str] | None:
    for x in byc[chrom]:
        if x[1] <= pos <= x[2]: return x
        if x[1] > pos: break
    return None

meta = pd.read_csv(ROOT / "data/imm/screening_map_results.csv", usecols=["PROBENNUMMER", "NR", "species"]).set_index("PROBENNUMMER")
sets: dict[str, set[tuple[str, int, str, str, str]]] = {}
for tab in sorted(RES.glob("Mkan329-*/5_typing/kansasii_snippy/snippy_out/snps.tab")):
    s = tab.parts[-5]
    if s not in meta.index or pd.isna(meta.loc[s, "species"]) or meta.loc[s, "NR"] in EXCL: continue
    t = pd.read_csv(tab, sep="\t", dtype=str, usecols=["CHROM", "POS", "TYPE", "REF", "ALT"])
    sets[s] = set(zip(t.CHROM, t.POS.astype(int), t.TYPE, t.REF, t.ALT))
species = {s: meta.loc[s, "species"] for s in sets}
print("genomes:", len(sets), pd.Series(species).value_counts().to_dict())

def describe(chrom: str, pos: int, typ: str, ref: str, alt: str) -> tuple[str, str, str, str]:
    x = hit(chrom, pos)
    if x is None: return "intergenic", "", "", ""
    _, st, en, strand, ft, gene, prod, lt = x
    if ft != "CDS": return f"{ft}", gene, prod, lt
    if typ in ("ins", "del"): return ("frameshift" if abs(len(alt) - len(ref)) % 3 else "inframe_indel"), gene, prod, lt
    if len(ref) != len(alt): return "complex_indel", gene, prod, lt
    cds = genome[chrom][st - 1:en]; cds = cds if strand == "+" else str(Seq(cds).reverse_complement())
    subs: dict[int, list[tuple[int, str]]] = {}
    for i, (r, a) in enumerate(zip(ref, alt)):
        if r == a or not st <= pos + i <= en: continue
        off = pos + i - st if strand == "+" else en - (pos + i)
        a2 = a if strand == "+" else str(Seq(a).reverse_complement())
        subs.setdefault(off // 3, []).append((off % 3, a2))
    out: list[str] = []
    for ci, lst in subs.items():
        c0 = cds[ci * 3:ci * 3 + 3]; c1 = list(c0)
        for k, a in lst: c1[k] = a
        a0, a1 = str(Seq(c0).translate(table=11)), str(Seq("".join(c1)).translate(table=11))
        out.append(f"{a0}{ci + 1}{a1}" + (" (syn)" if a0 == a1 else ""))
    eff = ",".join(out)
    if "*" in eff.split(str(ci + 1))[-1:][0] or re.search(r"\d+\*", eff): eff = eff + " STOP_GAINED"
    return eff, gene, prod, lt

rows: list[dict[str, Any]] = []
for foc, sp in FOCAL.items():
    others = [s for s in sets if species[s] == sp and s != foc]
    cnt: defaultdict[tuple[str, int, str, str, str], int] = defaultdict(int)
    for s in others:
        for v in sets[s]: cnt[v] += 1
    for v in sets[foc]:
        n = cnt.get(v, 0)
        if n <= 1:
            eff, gene, prod, lt = describe(*v)
            rows.append(dict(focal=foc, species=sp, n_species_bg=len(others), n_others_with=n, chrom=v[0], pos=v[1], type=v[2],
                             ref=v[3], alt=v[4], effect=eff, gene=gene, locus=lt, product=prod))
df = pd.DataFrame(rows); df.to_csv(OUT / "outlier_private_variants.csv", index=False)
df["impact"] = ~df.effect.isin(["intergenic", "rRNA", "tRNA", "ncRNA"]) & ~df.effect.str.contains("syn", na=False)
priv = df[df.n_others_with == 0]
print("\nvariants private to the focal genome (0 others) / carried by <=1 other / of those, protein-altering:")
for foc in FOCAL:
    a, b = priv[priv.focal == foc], df[df.focal == foc]
    print(f"  {foc} ({FOCAL[foc]}, bg {int(b.n_species_bg.iloc[0]) if len(b) else 0}): private={len(a)}  <=1 other={len(b)}  protein-altering private={int(a.impact.sum())}")
