#!/usr/bin/env bash
# 00_fetch_references.sh - download the M. tuberculosis H37Rv proteins (UniProt) and the
# E. coli rrnA region (NCBI) that scripts 04-07 and 09 use as numbering references.
#
# Usage: ./00_fetch_references.sh [OUTDIR]     (default $KANSASII_ROOT/output/gwas/ref)
# Existing files are kept. Output is identical to the files used for the 2026-10 analysis.
set -euo pipefail

ROOT=${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}
OUT=${1:-$ROOT/output/gwas/ref}
mkdir -p "$OUT"

uniprot () {  # $1 = UniProt accession, $2 = output name
  [[ -s "$OUT/$2" ]] && { echo "have $2"; return; }
  echo "fetching $1 -> $2"
  curl -fsSL "https://rest.uniprot.org/uniprotkb/$1.fasta" -o "$OUT/$2"
}

uniprot P9WGY9 Mtb_rpoB.faa   # RpoB   Rv0667
uniprot P9WG45 Mtb_gyrB.faa   # GyrB   Rv0005
uniprot P9WG47 Mtb_gyrA.faa   # GyrA   Rv0006
uniprot P9WIE5 Mtb_katG.faa   # KatG   Rv1908c
uniprot P9WNL7 Mtb_embB.faa   # EmbB   Rv3795
uniprot P9WH63 Mtb_rpsL.faa   # RpsL   Rv0682
uniprot P9WH87 Mtb_rplC.faa   # RplC   Rv0701
uniprot P9WNF9 Mtb_ethA.faa   # EthA   Rv3854c
uniprot P9WMC1 Mtb_ethR.faa   # EthR   Rv3855
uniprot P9WH07 Mtb_ksgA.faa   # RsmA   Rv1010
uniprot P9WGW7 Mtb_rsmI.faa   # RsmI   Rv3579c

if [[ ! -s "$OUT/Ecoli_rrnA_region.fa" ]]; then
  echo "fetching E. coli K-12 MG1655 rrnA region -> Ecoli_rrnA_region.fa"
  curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=U00096.3&rettype=fasta&retmode=text&seq_start=4035000&seq_stop=4041000" \
    -o "$OUT/Ecoli_rrnA_region.fa"
fi

for f in "$OUT"/*.fa*; do printf '  %-24s %s\n' "$(basename "$f")" "$(head -1 "$f" | cut -c1-70)"; done
