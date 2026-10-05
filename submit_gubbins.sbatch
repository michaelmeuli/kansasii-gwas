#!/usr/bin/env bash
#SBATCH --job-name=kansasii-gubbins
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
module load apptainer
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa
OUT=/shares/sander.imm.uzh/MM/kansasii/output/gwas/gubbins
python /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-gwas/scripts/16_build_gubbins_alignment.py
cd "$OUT"
# Gubbins 3.3.1 from the group shared container cache; FastTree for both tree stages (fast; approximate ML; FastTree supports GTRGAMMA, not GTR)
apptainer exec /shares/sander.imm.uzh/singu/gubbins_3.3.1--py310pl5321h83093d7_0.sif \
  run_gubbins.py --prefix kansasii_gubbins --threads 8 \
  --first-tree-builder fasttree --tree-builder fasttree --model GTRGAMMA \
  --filter-percentage 50 --min-snps 3 kansasii_chr_aln.fasta
