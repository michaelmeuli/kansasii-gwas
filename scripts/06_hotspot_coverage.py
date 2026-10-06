"""Base state and coverage at hotspot codons straight from snippy snps.aligned.fa ('-' = no coverage, 'N' = low coverage)."""
from __future__ import annotations

import os
import json
import pandas as pd
from pathlib import Path
from collections import Counter

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii")); OUT = ROOT / "output/gwas"
RES = ROOT / "runs/mkan329/assembly/results"
sites = json.load(open(OUT / "hotspot_sites.json", encoding="utf-8"))
ref: dict[str, list[str]] = {}
cur = ""
for l in open(RES / "Mkan329-001/5_typing/kansasii_snippy/snippy_out/ref.fa", encoding="utf-8"):
    if l.startswith(">"):
        cur = l[1:].split()[0]; ref[cur] = []
    else:
        ref[cur].append(l.strip())
chrom = "".join(ref["NC_022663.1"])
need = {p for v in sites.values() for p in v}

rows = []
for d in sorted(RES.glob("Mkan329-*")):
    f = d / "5_typing/kansasii_snippy/snippy_out/snps.aligned.fa"
    if not f.exists():
        continue
    seq, on = [], False
    for l in open(f, encoding="utf-8"):
        if l.startswith(">"):
            if on: break
            on = l[1:].split()[0] == "NC_022663.1"
        elif on:
            seq.append(l.strip())
    s = "".join(seq)
    for site, ps in sites.items():
        bases = "".join(s[p - 1] for p in ps)
        refb = "".join(chrom[p - 1] for p in ps)
        state = "no_coverage" if set(bases) <= {"-"} else "low_cov_N" if "N" in bases or "-" in bases else "ref" if bases.upper() == refb.upper() else "DIFFERS"
        rows.append(dict(sample=d.name, site=site, ref=refb, isolate=bases, state=state))
df = pd.DataFrame(rows); df.to_csv(OUT / "hotspot_coverage_states.csv", index=False)
print("genomes read:", df["sample"].nunique())
print(pd.crosstab(df.site, df.state).to_string())
print("\nDIFFERS:\n", df[df.state == "DIFFERS"].groupby(["site", "ref", "isolate"]).size().to_string())
