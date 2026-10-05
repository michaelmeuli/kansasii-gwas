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
| `13_kansasii_lineages.py` (+ `submit_lineages.sbatch`) | core-SNP distances, NJ tree and clades within kansasii (82 genomes) from snippy alignments | `lineage/core_snp_distance_matrix.csv`, `lineage/kansasii_nj.nwk`, `lineage/kansasii_core_snps.fasta` |
| `14_stratified_analysis.py` | amikacin and gene-level tests stratified by clade | `lineage/kansasii_clades.csv`, `lineage/cladeB_gene_tests.csv`, `lineage/mic_cohort_clade_vs_mic.csv` |
| `15_region_haplotype.py` | variant counts and haplotype classes in the 26 kb region FCNCLM_01713-01736 | `lineage/kansasii_region_haplotype.csv` |
| `16_build_gubbins_alignment.py`, `submit_gubbins.sbatch` | full-length chromosome alignment (82 kansasii + reference) and the Gubbins 3.3.1 run (shared container `/shares/sander.imm.uzh/singu/gubbins_3.3.1--py310pl5321h83093d7_0.sif`; FastTree, GTRGAMMA; about 41 min) | `gubbins/kansasii_gubbins.*` |
| `17_blast_region.py` | BLAST of the three region haplotypes against 72 public complex genomes | `blast/per_genome_summary.csv` |
| `18_gubbins_region_check.py`, `19_clonal_frame_pairs.py` | does Gubbins flag the region; clonal-frame nearest-relative test of haplotype vs amikacin | `gubbins/region_block_coverage_by_isolate.csv`, `gubbins/amikacin_pairs_table.csv` |
| `20_gubbins_stratified.py` | tree-aware re-analysis using the Gubbins recombination-filtered tree (replaces the raw-distance clades of 13-14) | `gubbins/filtered_tree_distance_matrix.csv`, `gubbins/mantel_mic_vs_filtered_tree.csv`, `gubbins/main_cluster_gene_tests.csv` |
| `21_segment_boundaries.py` | window-by-window identity of each carrier's own assembly to the reference vs public persicum genomes; boundaries in reference coordinates | `boundaries/window_calls.csv` |
| `22_reclassify_segment.py`, `23_segment_vs_amikacin.py` | reclassify all 82 kansasii genomes on the true extent (diagnostic sites); test classes against amikacin (MGIT, MIC, tree-aware) | `boundaries/segment_classes_v2.csv`, `boundaries/segment_genes.csv` |
| `24_rsmA_rsmI_alleles.py`, `25_segment_protein_divergence.py`, `26_rsm_alleles_vs_amikacin.py` | RsmA/RsmI alleles across cohort and public genomes; divergence ranking within the segment; alleles vs amikacin | `boundaries/RsmA_allele_per_isolate.csv`, `boundaries/RsmI_allele_per_isolate.csv`, `boundaries/segment_protein_identity.csv` |
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
  differing sites). (Correction: an earlier draft said this was consistent with several M. kansasii subtypes. It is not: the GTDB
  species 'kansasii' is the former subtype I only; the other subtypes are separate species, 88-94% ANI away.)
- The known amikacin/aminoglycoside loci show nothing: rrs (2 R carriers, 0 non-R, p = 0.51; the hotspot A1408 itself is wild-type), and
  no gene named eis or whiB7 is annotated in the reference.
- Independent check in the MIC cohort (21 kansasii, no overlap): amikacin MIC median 2.0 log2 for ANI below 99.67 (n = 13) vs 1.25 for ANI
  of 99.67 or more (n = 8), Mann-Whitney p = 0.054; Spearman across ANI r = -0.23, p = 0.32. Same direction, not confirmed. The threshold
  came from the MGIT data and only amikacin was tested for this. Ciprofloxacin also shows a weak trend with ANI (rho = -0.44, p = 0.047
  unadjusted across 13 drugs), which does not survive correction.
- Caveat on the phenotype: 1 mg/L is a low MGIT concentration, near the MIC of wild-type isolates (MIC cohort amikacin 1-8 mg/L), so
  "R" at 1 mg/L reflects a slightly higher MIC distribution, not clinical resistance.

### Lineages within kansasii and stratified analysis (scripts 13-15) - clades A/B superseded, see the Gubbins-tree section below
- All 85 GTDB-kansasii isolates are former subtype I, so there are no named subtypes to assign. Structure was derived instead from
  4.19 Mb of core sequence (sites with an unambiguous base in all 82 non-excluded kansasii genomes; 218,694 variable sites). Median
  pairwise distance is 22,588 core SNPs (max 54,123). The tree is neighbour-joining/average-linkage on SNP counts, without recombination
  filtering, so it is an approximation.
- Two main groups: clade A (n = 20, ANI to reference 99.14-99.84%) and clade B (n = 60, up to 100%), plus two single long-branch genomes
  (Mkan329-071, Mkan329-124). The 3-way "deep split" in the merge heights is just those two singletons.
- Amikacin 1 mg/L by clade (R / S+I): A 3/0, B 7/5, singletons 1/1 (A vs B Fisher p = 0.51). Clade does not explain the phenotype.
  MIC cohort amikacin: clade A median 2.0 vs B 1.5 log2 (p = 0.077; q = 0.88 over 13 drugs); no drug is significant.
- Within clade B (7 R, 5 S+I) a block of 32 gene features splits the groups (7/0 vs 0/5; q = 0.13 over 106 patterns; not significant). ANI
  within clade B also differs (R mean 99.44 vs S+I 99.85; p = 0.010 uncorrected).
- The block is mostly one 26 kb region, NC_022663.1:1,944,065-1,970,050 (FCNCLM_01713-01736: sigJ, an aminodeoxychorismate synthase,
  rsmI (16S rRNA C1402 2'-O-methyltransferase), pMT1 mannosyltransferase, arcA, a TetR regulator, aurF, a PPE protein and cell-wall
  glycosyltransferases), plus devS elsewhere. In the 7 R isolates the region has about 880 variants against the reference (0 in the 5
  non-R isolates; neighbouring windows have 0.2-440). They share 817 of them (pairwise Jaccard 0.97). The raw core-SNP distance between these isolates (about
  19,500) is recombination-inflated; see the Gubbins section, which shows the carriers are not one clade but are modestly closer than average.
- Across all 82 kansasii the region has several states: none (30 genomes), about 211 variants (25), about 400 (5, apparently combinations),
  about 880 (22). The same haplotype states appear to be polarised against the reference; the reference-like state is the minority outside
  clade B's non-resistant group.
- Haplotype ~880 vs amikacin 1 mg/L in all kansasii MGIT isolates: R in 8 of 9 carriers vs 3 of 8 others (Fisher p = 0.0498). A coarser
  "any divergence" class gave 10/11 vs 1/6 (p = 0.005) but mixes several haplotypes, so it is not interpretable as one event.
- Independent MIC cohort: haplotype ~880 (n = 4) median amikacin 2.0 log2 vs 1.5 in others (n = 17), Mann-Whitney p = 0.135. Same
  direction, underpowered, not confirmed.
- Statistical status: the region was chosen after searching 106 patterns within clade B, so its p-values are optimistic. It is a lead, not
  a finding. A functional link is speculative: rsmI and the cell-wall genes are plausible only in a general sense.

### Source of the 26 kb haplotype and Gubbins (scripts 16-19)
- BLAST (script 17) against 72 public complex genomes (kansasii 42, persicum 12, pseudokansasii 6, attenuatum 4, innocens 3, ostraviense 3,
  gastri 2). Query haplotypes were extracted from isolate assemblies (the snippy alignment is 17% uncovered in the divergent haplotype):
  - The ~880-variant haplotype (Mkan329-101, 25,529 bp) is 99.89% identical to M. persicum over its whole length (all 12 persicum genomes,
    99.81-99.93%). Against kansasii it is 94.1% on average, against pseudokansasii 90.2%, gastri 93.0%, innocens 93.5%, ostraviense 93.6%.
  - The reference-like haplotype is 93.5% identical to persicum, so this segment is persicum-derived DNA inside kansasii (an interspecies
    transfer), not a kansasii allele. The ~211-variant haplotype is a kansasii allele (98.8% within kansasii, 94.6% to persicum).
  - Only 2 of 42 public kansasii genomes carry the persicum-like segment (>= 99.5% identity); 27 of 38 carry the reference-like one. In this cohort
    22 of 82 kansasii genomes carry it.
- Gubbins final run (82 kansasii + reference; 1,275 recombination blocks, median 4,659 bp; 15.2 Mb summed over branches). Recombination is
  extensive, so raw SNP distances, and the raw-distance clades in scripts 13-14, are recombination-inflated.
  - The region is flagged in nearly every isolate (29 overlapping blocks, 10 on branches shared by 20 or more isolates). Fraction of the region
    covered: reference-like 0.93, ~211 haplotype 0.96, ~880 haplotype 0.88. So Gubbins marks it as a general recombination hotspot with several
    circulating haplotypes; it does not single out the carriers.
  - **Retraction:** at iteration 4 (interim) Gubbins predicted a 64.8 kb block shared by 20 carriers (monophyletic clade, 21/22 carriers vs 0/60
    others covered). That does not hold in the final run and was an artifact of recombinant SNPs still in the tree. Use only final files.
  - Final tree: the 22 carriers' common ancestor spans all 82 tips, so they are not one clade. They are modestly closer to each other than
    non-carriers are (median 79 vs 96 SNPs; range 19-189).
- Clonal-frame test (script 19, 17 phenotyped kansasii): resistant isolates are closer to each other (84 SNPs) than to S/I isolates (114), so the
  phenotype has lineage signal. Only 3 nearest-relative pairs differ in haplotype (2 distinct carriers, 81-104 SNPs apart); all 3 have the carrier
  resistant, sign test p = 0.25. This cannot separate the haplotype from the genetic background. Inconclusive.
- All 5 persicum isolates with a readable amikacin 1 mg/L result grew at 1 mg/L (persicum carries the donor haplotype by definition). This is
  consistent but cannot be counted as support because persicum differs from kansasii across the whole genome.

### Stratification with the Gubbins tree (script 20)
Replaces the clade A/B stratification, which used raw SNP distances inflated by recombination (15.2 Mb recombined, summed over branches).
- The recombination-filtered tree is shallow and star-like: median pairwise distance 97 SNPs (quartiles 76-121), median nearest-neighbour distance
  34. There are no deep clades: 6 clusters at 100 SNPs (one holds 66 of 81 isolates), 30 at 60 SNPs, 51 at 40 SNPs. Clusters are average-linkage on
  tree distances.
- Mkan329-068 is excluded from these analyses. Its terminal branch is 28,959 SNPs (median 18), but its raw nearest-neighbour distance is normal (9,407),
  its missing data is normal, its ANI is 99.68% and Gubbins found no private blocks in it. This is judged a tree artifact; the cause is not resolved.
- Phenotype has lineage signal (Mantel-type test, 5,000 permutations): amikacin at 1 mg/L in MGIT rho = 0.31, p = 0.013 (n = 17). In the MIC cohort
  ciprofloxacin (p = 0.021), doxycycline (0.022) and streptomycin (0.027) are nominal; none survives BH correction over 13 drugs (q = 0.118).
- The ~880-variant haplotype occurs on many backgrounds: the 22 carriers (of 81) fall in 13 different clusters at 40 SNPs and 9 at 60 SNPs; the
  largest 60-SNP cluster holds only 8 of them (with 15 non-carriers) and 14 carriers lie outside it. At fine scale carriers still cluster: 91% have a
  carrier as nearest relative (random expectation about 27%). This could reflect spread within lineages or residual recombinant signal in the tree.
- Stratified test of haplotype vs amikacin 1 mg/L: with strata = clusters of close relatives, there are no informative strata at 40 SNPs and one at
  60 and 100 SNPs. Exact Fisher p = 0.20 (4/0 vs 0/1) and 0.070 (7/0 vs 3/3). Not significant. The asymptotic Mantel-Haenszel p (0.04-0.05) is
  unreliable with a single small stratum and should not be quoted.
- Gene-level within the main 66-isolate cluster (13 phenotyped: 10 R, 3 S/I): nothing significant (147 patterns, min q = 0.51). The top pattern
  (10/0 vs 0/3) is a block of genes next to the introgressed segment (FCNCLM_01704-01709: rpfB, metG, tatD), meaning the segment likely extends
  beyond the region as I defined it.
- Net effect: the earlier conclusion that clade membership does not explain amikacin stands, but for a different reason. The clonal frame has little
  structure to stratify on, and with only 17 phenotyped isolates the data contain too few closely related pairs that differ in haplotype to separate
  it from background.

### Redefined extent of the persicum-derived segment (scripts 21-23)
This supersedes the 26 kb region (FCNCLM_01713-01736) used in scripts 15 and 17-20, and the numbers quoted for it (8/9 vs 3/8, p = 0.0498).
- Each carrier's own assembly was scanned over 1.80-2.10 Mb (1 kb windows, 500 bp step) against the reference and 12 public persicum genomes. A window
  is persicum-like when identity to persicum is at least 99.3% and at least 2 points above identity to the reference. Controls: two reference-like
  isolates (0 persicum-like windows) and two ~211-variant isolates.
- Full segment: NC_022663.1:1,929,793-1,994,501 (64.7 kb), seen with both boundaries inside the contig in Mkan329-001 and -116, and with the left
  boundary in five more (Mkan329-003, -013, -020, -066, -101; contig ends at 1,978,371) and the right boundary in two (-012, -022). It matches the
  64.75 kb block of the interim Gubbins iteration 4 (1,929,757-1,994,509) to within 40 bp. It contains 59 genes, FCNCLM_01702-01761 (16 in the left
  half, 43 in the right), including two 16S rRNA methyltransferases (rsmA, rsmI), a TetR regulator and cell-wall genes. Whether any matters for
  amikacin is unknown.
- The "~211-variant" haplotype is the left part only (about 1,929-1,950 kb, 21 kb): persicum-like to about 1,950,330 and kansasii-like to the right. The
  old 26 kb region straddled that internal boundary, which produced the 0 / ~211 / ~880 variant classes. Two control isolates also carry other
  persicum-like patches elsewhere (1.82-1.86 Mb in Mkan329-072, 2.03-2.05 Mb in Mkan329-005); not followed up.
- Reclassification of all 82 kansasii genomes with 2,926 diagnostic sites (840 left, 2,086 right); persicum-like threshold 0.75, none 0.10:
  full 21, left-only 27, none 26, partial/mosaic 8. Persicum-like DNA in the left part is therefore common (48 of 82). Nineteen left-only isolates share
  an identical left fraction of 0.86, a fixed shorter variant of the left piece.
- Amikacin 1 mg/L, MGIT (17 kansasii, 11 R). Full vs rest: R 7/8 vs 4/9, p = 0.13. Any persicum-like piece (full or left-only) vs none/partial: R 9/10
  vs 2/7, Fisher p = 0.0345 (nominal; four class definitions were compared here and an earlier one before, so this does not survive correction).
  Left-only vs rest: p = 0.51. By class: full 7R/1S-I, left-only 2R/0, none 0R/4, partial 2R/1.
- Independent MIC cohort (21 kansasii, amikacin log2 MIC): full median 2.0 (n = 4) vs others, p = 0.135; full or left-only vs others, p = 0.58 (left-only
  median 1.5). Not replicated for the broader class.
- Tree-aware (Gubbins filtered tree): nearest-relative pairs discordant for the full segment: 3, of which 1 supports and 1 goes the other way (sign
  test p = 1.0); for any persicum-like piece, 4 pairs, 1 vs 1. One informative stratum at each cut: full segment p = 0.20 (60 SNPs) and 0.19 (100
  SNPs); full or left-only p = 0.20 and 0.035 (the latter 8/0 vs 2/3, single stratum, exploratory).
- Net: the amikacin link is not supported by the tree-aware or the independent-cohort checks. It is a hypothesis with a nominal association in one
  cohort and an unresolved dependence on how carriers are defined.

### The two 16S rRNA methyltransferases in the segment: rsmA and rsmI (scripts 24-26)
RsmA (FCNCLM_01703, 317 aa; Mtb KsgA/Rv1010, 84% identical) and RsmI (FCNCLM_01716, 285 aa; Mtb Rv2372c, 84% identical) both lie in the left part of
the segment. Alleles were taken from each genome's own Bakta proteins (82 cohort kansasii) and checked with tblastn in 12 public persicum and 42 public
kansasii genomes.
- Persicum-like alleles keep the same length (317 and 285 aa); there is no truncation or frameshift in any allele.
- RsmA: main persicum-like allele (32 genomes) has 10 differences from the kansasii reference (96.8% identity): S121G, R123H, D129V, V188M, H191Y,
  Q293R, T295S, A299S, A308T, S315P. Seven of ten are in the C-terminal tail (293-315) and only two (R123H, S315P) are at positions where the reference
  equals Mtb. Rarer alleles carry more changes (up to 21, with a 3-residue shorter protein in 3 genomes).
- RsmI: main persicum-like allele (47 genomes) has 5 differences (98.2%): Q76R, T231S, Q235R, H243N, E280A; one (E280A) is at a position where the reference
  equals Mtb.
- Both alleles are persicum-like in all 12 public persicum genomes; among public kansasii genomes 7 of 42 (RsmA) and 11 of 42 (RsmI) carry them.
- Allele by segment class: RsmI is persicum-like in all 21 full carriers and 24 of 27 left-only isolates, whereas RsmA is persicum-like in all 21 full but
  only 9 of 27 left-only isolates (the other left-only isolates have a shorter left piece that starts after rsmA).
- Not unusual divergence: across the 59 segment proteins the median identity between the kansasii reference and a full carrier is 94.2% (quartiles
  89.8-96.5%). RsmA (96.8%) is more conserved than 81% of them and RsmI (98.2%) more than 93% of them. The differences look like ordinary kansasii-persicum
  divergence, not a damaged enzyme. A core-motif check was attempted but the search pattern matched neither protein, so nothing is reported on motifs;
  whether any difference affects catalysis is unknown.
- Gene-specific amikacin test (persicum-like vs reference-like allele): MGIT 1 mg/L: RsmA R 9/11 vs 2/6 (p = 0.11); RsmI R 10/11 vs 1/6 (Fisher p = 0.0054). MIC
  cohort: RsmA median 1.5 vs 2.0 log2 (p = 0.23); RsmI 1.75 vs 1.5 (p = 0.48). The RsmI MGIT result is essentially the "any persicum-like piece" association
  again (RsmI is persicum-like in nearly all left-only and full carriers), so it is not independent evidence. It was chosen after seeing the segment result,
  and about seven related tests have now been run on this hypothesis (old region, four class definitions, two alleles), so p = 0.0054 is nominal. It is
  not replicated in the MIC cohort, and the tree-aware check was not repeated for the allele variable.
- Net: no sequence-level reason to expect a functional change in either enzyme, and the allele-level association adds no independent support.

## Interpretation so far
No known resistance mutation explains any phenotype in either cohort. The patterns that recur are not single mutations: persicum isolate 054 is a
multi-drug high outlier, and amikacin growth at 1 mg/L in kansasii tracks a persicum-derived segment (64.7 kb, FCNCLM_01702-01761, about 99.9% identical to M. persicum) that sits in
a recombination hotspot and is carried by several backgrounds. The tree-aware checks do not support it (1 supporting and 1 contradicting nearest-relative pair; MIC cohort p = 0.135 for the full segment). Both need more isolates and repeat testing; neither is statistically confirmed.

## Open items
1. Retest 054 and 055 (all drugs). A repeat decides whether the multi-drug shift is real or a testing artifact.
2. Test whether position 3145420 changes ethA expression or ethionamide MIC (e.g. in other persicum isolates or by reporter assay).
3. Look at permeability and regulator genes in 054 beyond the filtered variant lists. Read-mapping calls in divergent genes need
   verification against the assembly.
4. Test the 26 kb persicum-derived haplotype in new isolates. Redo the clade stratification (scripts 13-14) with the Gubbins clonal-frame tree
   instead of raw SNP distances. Add more closely related isolate pairs that differ in haplotype. MGIT cannot validate the MIC outliers
   (none have MGIT data).
5. Calibrate MGIT 1 mg/L against MIC: almost all MIC-cohort isolates have amikacin MIC of 2 mg/L or more, yet 5 MGIT isolates are
   S/I at 1 mg/L, so the two methods may not be on a comparable scale.
6. Persicum needs more MIC-tested isolates; per-species analysis is currently impossible there.
