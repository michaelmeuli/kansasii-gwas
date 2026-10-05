#!/usr/bin/env bash
#SBATCH --job-name=kansasii-mgit-gene
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa
python /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-gwas/scripts/11_mgit_whole_gene.py
