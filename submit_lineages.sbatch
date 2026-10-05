#!/usr/bin/env bash
#SBATCH --job-name=kansasii-lineages
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mlsa
python /shares/sander.imm.uzh/MM/kansasii/repos/kansasii-gwas/scripts/13_kansasii_lineages.py
