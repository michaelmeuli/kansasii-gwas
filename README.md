# kansasii-gwas: genotype vs drug susceptibility in the M. kansasii complex

Status as of 2026-10-05. Exploratory analysis; nothing here is a validated resistance finding.
Run scripts in `env_mlsa` (`module load miniforge3; source activate env_mlsa`). Outputs go to `output/gwas/`.

## Cohort that can be analysed
- 44 isolates have broth-microdilution MIC (`output/mic/mhk/mic_parsed.csv`, 13 drugs); 30 also have a genome and are not excluded
  (21 kansasii, 5 persicum, 2 pseudokansasii, 1 attenuatum, 1 innocens).
- Excluded as mixed or contaminated cultures: NR 35, 113, 106, 10, 21. Control strains are dropped.
- Almost all isolates are susceptible by CLSI category, so continuous log2 MIC is used. Censored values (`<=`, `>`) are treated as exact.
- Persicum has only 5 isolates with MIC. Per-species statistics are not possible for it.

## Scripts (`scripts/`)
| Script | Purpose | Main output |
|---|---|---|
| `01_build_table.py` | join MIC with species and QC flags, apply exclusions | `geno_pheno_table_long.csv` |
| `02_candidate_variants.py` | find candidate resistance loci in the Bakta reference annotation; list variants per isolate | `candidate_gene_coords.csv`, `candidate_variants_per_isolate.csv` |
| `03_within_species_mic.py` | carrier vs non-carrier MIC tests per species (Mann-Whitney) | `within_species_mic_tests.csv`, `candidate_features_per_sample.csv` |
| `04_hotspots.py` | known hotspot codons mapped to M. tuberculosis numbering, verified by residue identity | `hotspot_numbering_check.csv`, `hotspot_genotypes_per_isolate.csv` |
| `05_rrna_hotspots.py` | rrl A2058/A2059 and rrs A1408 in E. coli numbering | (printed) |
| `06_hotspot_coverage.py` | base state and coverage at hotspots from snippy `snps.aligned.fa` | `hotspot_coverage_states.csv` |
| `07_gyrA_from_assemblies.py` | gyrA QRDR residues read from each genome's own annotation | `gyrA_QRDR_from_assemblies.csv` |
| `08_outlier_private_variants.py` | variants carried by an MIC-outlier isolate and by no other genome of the species | `outlier_private_variants.csv` |
| `09_ethA_054.py`, `10_ethA_promoter_and_paralog.py` | follow-up on Mkan329-054's ethA calls | `ethionamide_ethA_upstream_3145420.csv` |
| `11_mgit_whole_gene.py` (+ `submit_mgit_whole_gene.sbatch`) | MGIT whole-gene association (gene-level binary calls, Fisher within species) | `mgit/gene_matrix.csv`, `mgit/mgit_gene_tests.csv`, `mgit/mgit_tests_run.csv` |
| `12_mgit_lineage_followup.py` | characterise the top MGIT signal; ANI vs amikacin in the MIC cohort | `mgit/kansasii_mgit_block_ani_amikacin.csv`, `mgit/mic_cohort_ani_vs_mic.csv` |

Reference files downloaded for numbering (UniProt H37Rv proteins, E. coli rrnA region) are in `output/gwas/ref/`.

## Technical pitfalls found (these affect any rerun)
1. Snippy `snps.tab` has empty gene/effect columns because the snippy `ref.gff` is 0 bytes. Genes are mapped with the Bakta GFF3 of
   GCF_000157895.3 (same sequence; `contig_1` = NC_022663.1, `contig_2` = plasmid NC_022654.1).
2. The reference strain carries rpsL R43 while nearly every other isolate has K43 (the M. tuberculosis wild-type residue). "Differs
   from reference" is therefore not "mutant" at some sites.
3. kansasii gyrA is 1,258 aa because it contains an intein; other species have a normal ~837 aa gene named "DNA topoisomerase
   (ATP-hydrolyzing) subunit A". QRDR numbering still matches M. tuberculosis (A90, D94). Search both product names.
4. Snippy masks gyrA in persicum, pseudokansasii and attenuatum (N calls), so read gyrA from the assemblies for those.
5. rpoB clinical numbering (S450, H445, D435) equals the UniProt index minus 6. The Bakta kansasii rpoB lacks about 21 N-terminal
   residues, so S450 is residue 429. Always map by alignment.
6. Four loci in the reference are named ethA. Only `FCNCLM_02843` (87% identical to M. tuberculosis EthA) is the real ortholog.
   `FCNCLM_03714` (49%), `FCNCLM_02718` (52%) and `FCNCLM_04416` (74%) are other FAD monooxygenases.

## Results
### Candidate-gene screen (script 03)
- 1,282 features tested. Within kansasii (21 isolates, 507 tests) the smallest FDR-adjusted q is 0.75. Within persicum (5 isolates) the
  smallest q is 0.61. Nothing is significant.
- The pooled analysis with species-centred MIC produced small q-values (about 0.01) that are an artifact: rifabutin barely varies, and
  variants fixed in a species line up with species. Do not use it.
- Many features are carried by exactly the same isolates, so the real number of independent tests is much smaller than the count.

### Known hotspots (scripts 04-07)
- No resistance-associated change at gyrA A90/D94 (and 88, 89, 91), gyrB D461/N499/E501, rpoB D435/H445/S450/L452, embB M306/G406/Q497,
  katG S315, rpsL K88, rplC C154, rrl A2058/A2059 or rrs A1408 in any genome that has coverage.
- The only differences found are a synonymous rpoB S450 codon change (TCG to TCC) in one isolate and one ambiguous base at rpsL K88.
- gyrA was read from the assembly annotation for 127 of 131 genomes: wild-type at all QRDR residues. The 4 without an annotated gyrA are
  NR 10, 21, 113 (excluded) and NR 55 (not checked).
- rRNA sites have no read coverage in 22 genomes (rrl) and 10 (rrs); these are not called wild-type.

### MIC outlier tails (script 08)
- 27 outlier calls, 17 in genome-having isolates (7 isolates).
- Mkan329-054 (persicum) is above the cohort median for all 13 drugs and by more than 1.5 log2 steps for 7 (ciprofloxacin, linezolid,
  clarithromycin +3; moxifloxacin, ethionamide +5; isoniazid, rifampicin +2). Rifabutin is unchanged. A uniform rise across drug classes
  suggests a general mechanism (efflux, permeability, regulation) or a testing artifact; the data cannot separate them.
- Mkan329-055 (persicum) is low for 4 drugs. It has only 97 private variants and none that disrupt a gene.
- Persicum isolates carry long private stretches (054: about 900, 009: about 1,100) because they sit on separate branches, so counts and
  gene-class keyword searches do not discriminate.
- None of the outlier isolates has MGIT data.

### Mkan329-054 and ethA (scripts 09-10)
- **Retracted:** an earlier interim note said 054 has private ethA changes (E492D; A497T/V496I/P495A). Those are in `FCNCLM_03714`, a
  paralog, not ethA. 054's assembled protein for that locus is identical to every other persicum, and the calls sit in a low-depth
  cluster (11 reads vs about 20). They are treated as a mapping artifact.
- The real EthA protein (`FCNCLM_02843`) is identical in all 22 persicum genomes, including 054. No coding change in ethA or ethR is specific to 054.
- One candidate remains: genome position 3145420 (A to C; T to G on the ethA strand), 4 bp upstream of the ethA start codon.
  Read evidence is clean (054: 20 of 20 reads; 009: 19 of 19). Among the 28 isolates with ethionamide MIC, only 054 (MIC 10 mg/L) and
  009 (0.6 mg/L) carry it; they are the two highest of the 5 persicum with MIC. Another change at the same position, A to T, occurs in one
  innocens isolate with a low MIC.
  Caveats: the carriers are two isolates, both persicum; the unadjusted p = 0.001 is not meaningful because the site was chosen after
  looking at the outliers and is confounded with species; within persicum 2 vs 3 isolates cannot reach p below 0.1; the two carriers
  differ 16-fold in MIC, so the change would not explain 054's extreme value on its own. It does not explain 054's resistance to the
  other drugs.

### MGIT whole-gene analysis (scripts 11-12)
- 66 isolates have MGIT results; 47 have a genome and are not excluded (35 kansasii, 8 persicum, 4 pseudokansasii). Only 1 overlaps
  with the MIC cohort (a persicum isolate), so MGIT is an independent set. Unreadable results (`U`) are dropped.
- Resistance is rare. Only one drug and concentration has enough spread for a test (3 or more R and 3 or more non-R in a species):
  amikacin 1 mg/L in kansasii, 11 R vs 6 S/I. Isoniazid 0.1 mg/L is R in all 13 tested isolates (no variation); isoniazid 1 mg/L has 2 R;
  every other pair has 0-1 R, apart from clofazimine 0.25 mg/L (7 isolates, 4 R, all persicum). Persicum amikacin and the R-vs-S-only
  contrast did not reach the group-size minimum.
- Features: 8,321 gene-level binary calls (protein-altering variant, rRNA variant, or any change in the 150 bp upstream region). Fisher
  exact test, BH correction over distinct carrier patterns (271 patterns).
- Result: nothing significant (smallest q = 0.15). The top pattern is a block of 36 genes (9 across all 47 isolates) carried by 10 of 11
  resistant and 0 of 6 non-resistant isolates.
- That block is a divergence-from-reference signal, not a gene: carriers have ANI 98.9-99.7% to the reference, non-carriers 99.1-100%;
  carriers differ from each other by about 26,000 sites, as much as from non-carriers, so they are not one clone. ANI alone separates
  amikacin R from S/I (Mann-Whitney p = 0.05, uncorrected). The kansasii here is heterogeneous (82 genomes, 54 clusters at 1,000
  differing sites), consistent with several M. kansasii subtypes; the reference represents only one.
- The known amikacin/aminoglycoside loci show nothing: rrs (2 R carriers, 0 non-R, p = 0.51; the hotspot A1408 itself is wild-type), and
  no gene named eis or whiB7 is annotated in the reference.
- Independent check in the MIC cohort (21 kansasii, no overlap): amikacin MIC median 2.0 log2 for ANI below 99.67 (n = 13) vs 1.25 for ANI
  of 99.67 or more (n = 8), Mann-Whitney p = 0.054; Spearman across ANI r = -0.23, p = 0.32. Same direction, not confirmed. The threshold
  came from the MGIT data and only amikacin was tested for this. Ciprofloxacin also shows a weak trend with ANI (rho = -0.44, p = 0.047
  unadjusted across 13 drugs), which does not survive correction.
- Caveat on the phenotype: 1 mg/L is a low MGIT concentration, near the MIC of wild-type isolates (MIC cohort amikacin 1-8 mg/L), so
  "R" at 1 mg/L reflects a slightly higher MIC distribution, not clinical resistance.

## Interpretation so far
No known resistance mutation explains any phenotype in either cohort. The two reproducible patterns are lineage-level: kansasii isolates
more divergent from the reference tend toward higher amikacin MIC and growth at 1 mg/L amikacin, and persicum isolate 054 is a
multi-drug high outlier. Both need subtype assignment (e.g. hsp65/ITS or a within-species tree), more isolates, and repeat testing.

## Open items
1. Retest 054 and 055 (all drugs). A repeat decides whether the multi-drug shift is real or a testing artifact.
2. Test whether position 3145420 changes ethA expression or ethionamide MIC (e.g. in other persicum isolates or by reporter assay).
3. Look at permeability and regulator genes in 054 beyond the filtered variant lists. Read-mapping calls in divergent genes need
   verification against the assembly.
4. Assign M. kansasii subtype (hsp65/ITS or a relaxed-core within-species tree) and test subtype, not ANI, against amikacin; then redo the
   gene-level analysis with subtype as a stratum. MGIT cannot validate the MIC outliers (none have MGIT data).
5. Persicum needs more MIC-tested isolates; per-species analysis is currently impossible there.
