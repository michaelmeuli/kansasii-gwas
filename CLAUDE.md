# CLAUDE.md

See README.md; it is superseded by `kansasii-amr-marker-screen` and kept for the record. Easy to get wrong:

- **Run:** `sbatch run_all.sbatch` (scripts 00-26 in numeric order, Gubbins after 16). Not on the login node.
- **Environment:** `env_mlsa` (see `environment.yml`); blast tools come from `env_arm` via PATH, set in the sbatch file.
- **Inputs must be rerun first:** `data/imm/screening_map_results.csv` and `output/mic/{mhk,mgit}/` (run order in `immensekansasii/WORKFLOW.md`).
- **Outputs** go to `output/gwas/`, not the repo. The earlier run is archived in `output-old/gwas_2026-10-06_pre-jc/`.
- Usable N is small (see "Cohort that can be analysed"); do not report any lead as a finding.
- The sbatch file writes `output/gwas/MANIFEST.txt` at the end.
- Commit directly on `master`; push only when asked.
