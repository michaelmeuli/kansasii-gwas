"""Extract variants in candidate resistance genes per isolate from snippy snps.tab (unannotated),
using the Bakta annotation of the reference (identical sequence; contig_1==NC_022663.1, contig_2==NC_022654.1)."""
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii")
REFDIR = ROOT / "runs/kansasii_complex_gtdb_representatives/assembly/results/GCF_000157895.3_query/2_annotation"
RES = ROOT / "runs/mkan329/assembly/results"
OUT = ROOT / "output/gwas"
CONTIG = {"contig_1": "NC_022663.1", "contig_2": "NC_022654.1"}
CANDIDATES = {"gyrA", "gyrB", "rpoB", "rpsL", "embA", "embB", "embC", "katG", "inhA", "ethA", "ethR",
              "whiB7", "eis", "folP", "dfrA", "rplC", "rrl", "rrs", "erm"}

genome = {CONTIG[r.id]: str(r.seq) for r in SeqIO.parse(REFDIR / "GCF_000157895.3_query.fna", "fasta")}

genes = []
for line in open(REFDIR / "GCF_000157895.3_query.gff3"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("CDS", "rRNA"):
        continue
    attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    name = attrs.get("gene", "")
    if re.sub(r"\d+$", "", name) in CANDIDATES or name in CANDIDATES:
        genes.append(dict(gene=name, locus=attrs.get("locus_tag"), chrom=CONTIG[f[0]], start=int(f[3]),
                          end=int(f[4]), strand=f[6], type=f[2], product=attrs.get("product", "")))
genes = pd.DataFrame(genes)
genes.to_csv(OUT / "candidate_gene_coords.csv", index=False)
print("candidate loci found:\n", genes[["gene", "chrom", "start", "end", "strand", "type"]].to_string(index=False))

def codon_change(g, pos, ref, alt):
    """Amino-acid change for a single-base substitution inside CDS g."""
    seq = genome[g.chrom]
    cds = seq[g.start - 1:g.end]
    off = pos - g.start if g.strand == "+" else g.end - pos
    cds = cds if g.strand == "+" else str(Seq(cds).reverse_complement())
    r, a = (ref, alt) if g.strand == "+" else (str(Seq(ref).reverse_complement()), str(Seq(alt).reverse_complement()))
    ci = off // 3
    codon = cds[ci * 3:ci * 3 + 3]
    new = codon[:off % 3] + a + codon[off % 3 + 1:]
    aa0, aa1 = str(Seq(codon).translate(table=11)), str(Seq(new).translate(table=11))
    return ci + 1, aa0, aa1

rows = []
for tab in sorted(RES.glob("Mkan329-*/5_typing/kansasii_snippy/snippy_out/snps.tab")):
    sample = tab.parts[-5]
    t = pd.read_csv(tab, sep="\t", dtype=str)
    t["POS"] = t["POS"].astype(int)
    for g in genes.itertuples():
        sub = t[(t.CHROM == g.chrom) & (t.POS >= g.start) & (t.POS <= g.end)]
        for v in sub.itertuples():
            eff = ""
            if g.type == "CDS":
                if v.TYPE == "snp":
                    ci, a0, a1 = codon_change(g, v.POS, v.REF, v.ALT)
                    eff = f"{a0}{ci}{a1}" + ("" if a0 != a1 else " (syn)")
                elif v.TYPE in ("ins", "del"):
                    d = abs(len(v.ALT) - len(v.REF))
                    eff = "frameshift" if d % 3 else "inframe_indel"
                else:
                    eff = "complex/mnp"
            rows.append(dict(sample=sample, gene=g.gene, locus=g.locus, pos=v.POS, type=v.TYPE,
                             ref=v.REF, alt=v.ALT, effect=eff))
out = pd.DataFrame(rows)
out.to_csv(OUT / "candidate_variants_per_isolate.csv", index=False)
print("\nvariants:", len(out), "| isolates scanned:", len(list(RES.glob('Mkan329-*/5_typing/kansasii_snippy/snippy_out/snps.tab'))))
print(out[out.effect.ne("") & ~out.effect.str.contains("syn", na=False)].groupby(["gene", "effect"])["sample"].nunique()
      .sort_values(ascending=False).head(40).to_string())
