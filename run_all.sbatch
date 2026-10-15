#!/usr/bin/env bash
#SBATCH --job-name=kansasii-gwas
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log
# Whole genotype-phenotype analysis: scripts/00 .. 26 in numeric order (that order is a valid
# dependency order), with the Gubbins run between 16 and 17. Stops at the first failure.
#
#   sbatch run_all.sbatch            # everything (about 2-3 h, gubbins ~40 min)
#   FROM=13 sbatch run_all.sbatch    # restart at script 13 (earlier outputs must exist)
#
# Inputs  : data/imm/screening_map_results.csv, output/mic/{mhk,mgit}/, runs/mkan329/assembly/results
# Outputs : $KANSASII_ROOT/output/gwas/  (KANSASII_ROOT defaults to the cluster root)
set -euo pipefail
module load apptainer
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa
# blastn/makeblastdb/tblastn (scripts 17, 21, 24) are not in env_mlsa; take them from env_arm
if ! command -v blastn >/dev/null; then
  export PATH="$PATH:/home/mimeul/data/conda/envs/env_arm/bin"
fi
command -v blastn >/dev/null || { echo "blastn not found (env_arm)" >&2; exit 1; }

HERE=/shares/sander.imm.uzh/MM/kansasii/repos/kansasii-gwas
ROOT=${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}
GUBBINS_SIF=/shares/sander.imm.uzh/singu/gubbins_3.3.1--py310pl5321h83093d7_0.sif
FROM=${FROM:-0}
cd "$HERE"
mkdir -p "$ROOT/output/gwas"/{mgit,lineage,gubbins,blast,boundaries}

for s in scripts/[0-9][0-9]_*; do
  n=$((10#$(basename "$s" | cut -c1-2)))
  (( n < FROM )) && continue
  echo "=== $s ($(date +%T))"
  case "$s" in
    *.sh) bash "$s" ;;
    *.py) python "$s" ;;
  esac
  if (( n == 16 )); then
    echo "=== gubbins ($(date +%T))"
    cd "$ROOT/output/gwas/gubbins"
    # FastTree for both tree stages (approximate ML; supports GTRGAMMA, not GTR)
    apptainer exec "$GUBBINS_SIF" run_gubbins.py --prefix kansasii_gubbins --threads 8 \
      --first-tree-builder fasttree --tree-builder fasttree --model GTRGAMMA \
      --filter-percentage 50 --min-snps 3 kansasii_chr_aln.fasta
    [[ -s kansasii_gubbins.final_tree.tre ]] || { echo "gubbins produced no final tree" >&2; exit 1; }
    cd "$HERE"
  fi
done
echo "done $(date +%T)"
