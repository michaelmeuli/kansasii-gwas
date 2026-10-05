"""BLAST the three haplotypes of the 26 kb region (NC_022663.1:1,944,065-1,970,050) against the 72 public M. kansasii-complex genomes.
Queries: H880 (Mkan329-101), H211 (Mkan329-005), reference-like (reference sequence). For each query and subject genome: query bases covered by a
non-overlapping chain of HSPs and the identity over those bases."""
import subprocess
import pandas as pd
from pathlib import Path
from Bio import SeqIO

ROOT = Path("/shares/sander.imm.uzh/MM/kansasii"); OUT = ROOT / "output/gwas/blast"; OUT.mkdir(exist_ok=True)
RES = ROOT / "runs/mkan329/assembly/results"
PUB = ROOT / "data/gtdb_genomes/Mycobacteriaceae/mlsa-kansasii"
LO, HI = 1944065, 1970050
ref = str(next(SeqIO.parse(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa", "fasta")).seq).upper()
refseg = ref[LO - 1:HI]
def hap_from_assembly(sample):
    """Extract the region from the isolate's own assembly: BLAST the reference segment, take the contig with most aligned bases and span its hits."""
    fa = RES / sample / f"1_unicycler/{sample}.fasta"
    tmp = OUT / f"tmp_{sample}"; tmp.mkdir(exist_ok=True)
    (tmp / "q.fa").write_text(f">q\n{refseg}\n")
    subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "nucl", "-out", str(tmp / "db")], check=True, stdout=subprocess.DEVNULL)
    out = subprocess.run(["blastn", "-query", str(tmp / "q.fa"), "-db", str(tmp / "db"), "-evalue", "1e-10", "-outfmt", "6 sseqid sstart send length pident"],
                         check=True, capture_output=True, text=True).stdout
    h = pd.DataFrame([l.split("\t") for l in out.strip().split("\n")], columns=["c", "ss", "se", "len", "pid"]).astype({"ss": int, "se": int, "len": int, "pid": float})
    best = h.groupby("c").len.sum().idxmax(); d = h[h.c == best]
    a, b = min(d.ss.min(), d.se.min()), max(d.ss.max(), d.se.max())
    contigs = {r.id: str(r.seq).upper() for r in SeqIO.parse(fa, "fasta")}
    seq = contigs[best][a - 1:b]
    if d.iloc[0].se < d.iloc[0].ss:   # minus strand
        from Bio.Seq import Seq
        seq = str(Seq(seq).reverse_complement())
    print(f"  {sample}: contig {best}, extracted {len(seq)} bp (reference segment {len(refseg)} bp), hits {len(d)}, other contigs with hits: {h.c.nunique() - 1}")
    return seq
q = {"H880_Mkan329-101": hap_from_assembly("Mkan329-101"), "H211_Mkan329-005": hap_from_assembly("Mkan329-005"), "refpattern_ATCC12478": refseg}
with open(OUT / "queries.fasta", "w") as f:
    for k, v in q.items():
        print(k, "length", len(v), "non-ACGT:", sum(c not in "ACGT" for c in v)); f.write(f">{k}\n{v}\n")
# subject database with species-tagged contig ids
with open(OUT / "subjects.fasta", "w") as f:
    n = 0
    for sp in sorted(p.name for p in PUB.iterdir()):
        for fa in sorted((PUB / sp).glob("*.fna")):
            for r in SeqIO.parse(fa, "fasta"):
                f.write(f">{sp}|{fa.stem}|{r.id}\n{str(r.seq)}\n"); n += 1
print("subject contigs:", n)
subprocess.run(["makeblastdb", "-in", str(OUT / "subjects.fasta"), "-dbtype", "nucl", "-out", str(OUT / "subjects")], check=True, stdout=subprocess.DEVNULL)
subprocess.run(["blastn", "-query", str(OUT / "queries.fasta"), "-db", str(OUT / "subjects"), "-evalue", "1e-20", "-max_target_seqs", "2000", "-max_hsps", "50",
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore", "-out", str(OUT / "hits.tsv"), "-num_threads", "4"], check=True)
h = pd.read_csv(OUT / "hits.tsv", sep="\t", names="q s pident length qs qe ss se bits".split())
h[["species", "genome", "contig"]] = h.s.str.split("|", expand=True)
rows = []
for (qq, sp, g), d in h.groupby(["q", "species", "genome"]):
    cov = []; ident = 0.0; covered = 0
    for r in d.sort_values("bits", ascending=False).itertuples():
        a, b = min(r.qs, r.qe), max(r.qs, r.qe)
        if any(not (b < x or a > y) for x, y in cov): continue
        cov.append((a, b)); covered += b - a + 1; ident += r.pident * (b - a + 1)
    rows.append(dict(query=qq, species=sp, genome=g, covered_bp=covered, cov_frac=covered / len(q[qq]), identity=ident / covered if covered else float("nan"), n_hsp=len(cov)))
R = pd.DataFrame(rows); R.to_csv(OUT / "per_genome_summary.csv", index=False)
pd.set_option("display.width", 200)
for qq in q:
    x = R[R["query"] == qq]
    print(f"\n=== {qq}: genomes with >=80% of the region covered: {(x.cov_frac >= .8).sum()} of {x.genome.nunique()}")
    print(x[x.cov_frac >= .8].groupby("species").agg(n=("genome", "size"), identity_mean=("identity", "mean"), identity_min=("identity", "min"), identity_max=("identity", "max")).round(2).to_string())
