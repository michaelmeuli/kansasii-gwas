# kansasii-gwas: genotype vs drug susceptibility in the M. kansasii complex

> **Status: superseded by `kansasii-amr-marker-screen`; no positive findings.** This was an exploratory genotype-phenotype analysis, not a true
> GWAS. With 80 usable isolates (59 kansasii), almost no CLSI resistance, extensive recombination and strong species/lineage structure,
> the cohort is too small for genome-wide association. No hotspot resistance mutations were found, and no association survived multiple
> testing or replication in the independent cohort. The persicum-segment/amikacin and ethA-upstream leads are unconfirmed hypotheses. Use
> `kansasii-amr-marker-screen` (targeted AMR markers vs MIC) for resistance work; this repo is kept for the record, the numbering/pitfall notes below
> and the negative hotspot result.

Status as of 2026-10-16: rerun from scratch with `run_all.sbatch` on the full cohort (186 genomes, 83 isolates with a
broth-microdilution MIC, 66 with MGIT). Exploratory analysis; nothing here is a validated resistance finding. The earlier run
(2026-10-05/06: 131 genomes, 44 MIC isolates) is archived in `output-old/gwas_2026-10-06_pre-jc/`; several of its leads did not
hold on the larger cohort (see "What changed").

Outputs go to `output/gwas/`.

## Run

```bash
sbatch run_all.sbatch            # scripts/00 .. 26 in numeric order, Gubbins after 16 (about 2-3 h)
FROM=13 sbatch run_all.sbatch    # restart at script 13 (earlier outputs must exist)
```

The numeric order is a valid dependency order. `run_all.sbatch` uses `env_mlsa` (python, biopython,
statsmodels) and takes `blastn`/`makeblastdb`/`tblastn` (scripts 17, 21, 24) from `env_arm`. It needs
`data/imm/screening_map_results.csv` and `output/mic/{mhk,mgit}/` (rerun those first, see
`immensekansasii/WORKFLOW.md`, "Run order"). Script 00 downloads the numbering references to
`output/gwas/ref/`; the scripts are otherwise unchanged by it.

## Cohort that can be analysed
- 83 isolates have broth-microdilution MIC (`output/mic/mhk/mic_parsed.csv`, 13 drugs, includes the jc MIC table); 80 of them also have a genome
  and are not excluded: 59 kansasii, 12 persicum, 6 pseudokansasii, 1 attenuatum, 1 innocens and 1 avium (NR 184, not a complex member, see
  "Open items"). 66 isolates have MGIT results and a genome; the two cohorts overlap in 15 kansasii isolates.
- Excluded in `01_build_table.py`: isolates with `contamination_flag == "contaminated"` (CheckM > 5 %: NR 10, 21, 35, 106, 113; one rule for all scripts via `scripts/exclusions.py`) and control strains. NR 10 and 21 (7.7 %) were excluded by hand before 2026-10-20, included 2026-10-20 (flag at 10 %), and are excluded again by the 5 % flag.
- Almost all isolates are susceptible by CLSI category, so continuous log2 MIC is used. Censored values (`<=`, `>`) are treated as exact.
- Persicum has 12 isolates with MIC, the other non-kansasii species 1-6 each; per-species statistics are only meaningful for kansasii.

## Scripts (`scripts/`)
| Script | Purpose | Main output |
|---|---|---|
| `00_fetch_references.sh` | UniProt H37Rv proteins and E. coli rrnA region used for numbering | `ref/*.faa`, `ref/Ecoli_rrnA_region.fa` |
| `01_build_table.py` | join MIC with species and QC flags, apply exclusions | `geno_pheno_table_long.csv` |
| `02_candidate_variants.py` | find candidate resistance loci in the Bakta reference annotation; list variants per isolate | `candidate_gene_coords.csv`, `candidate_variants_per_isolate.csv` |
| `03_within_species_mic.py` | carrier vs non-carrier MIC tests per species (Mann-Whitney) | `within_species_mic_tests.csv`, `candidate_features_per_sample.csv` |
| `04_hotspots.py` | known hotspot codons mapped to M. tuberculosis numbering, verified by residue identity | `hotspot_numbering_check.csv`, `hotspot_genotypes_per_isolate.csv`, `hotspot_sites.json` |
| `05_rrna_hotspots.py` | rrl A2058/A2059 and rrs A1408 in E. coli numbering; adds the rRNA sites to `hotspot_sites.json` | (printed) |
| `06_hotspot_coverage.py` | base state and coverage at hotspots from snippy `snps.aligned.fa` | `hotspot_coverage_states.csv` |
| `07_gyrA_from_assemblies.py` | gyrA QRDR residues read from each genome's own annotation | `gyrA_QRDR_from_assemblies.csv` |
| `08_outlier_private_variants.py` | variants carried by an MIC-outlier isolate and by no other genome of the species | `outlier_private_variants.csv` |
| `09_ethA_054.py`, `10_ethA_promoter_and_paralog.py` | follow-up on Mkan329-054's ethA calls | `ethionamide_ethA_upstream_3145420.csv` |
| `11_mgit_whole_gene.py` | MGIT whole-gene association (gene-level binary calls, Fisher within species) | `mgit/gene_matrix.csv`, `mgit/mgit_gene_tests.csv`, `mgit/mgit_tests_run.csv` |
| `13_kansasii_lineages.py` | core-SNP distances, NJ tree and clades within kansasii (82 genomes) from snippy alignments | `lineage/core_snp_distance_matrix.csv`, `lineage/kansasii_nj.nwk`, `lineage/kansasii_core_snps.fasta` |
| `14_stratified_analysis.py` | amikacin and gene-level tests stratified by clade | `lineage/kansasii_clades.csv`, `lineage/cladeB_gene_tests.csv`, `lineage/mic_cohort_clade_vs_mic.csv` |
| `15_region_haplotype.py` | variant counts and haplotype classes (`haplotype` divergent/ref-like; `hap` none / H~211 / H~400(mix) / H~880, cut points 50/250/600 variants) in the 26 kb region FCNCLM_01713-01736 | `lineage/kansasii_region_haplotype.csv` |
| `16_build_gubbins_alignment.py` + Gubbins step of `run_all.sbatch` | full-length chromosome alignment (82 kansasii + reference) and the Gubbins 3.3.1 run (shared container `/shares/sander.imm.uzh/singu/gubbins_3.3.1--py310pl5321h83093d7_0.sif`; FastTree, GTRGAMMA; about 41 min) | `gubbins/kansasii_gubbins.*` |
| `17_blast_region.py` | BLAST of the three region haplotypes against 72 public complex genomes | `blast/per_genome_summary.csv` |
| `18_gubbins_region_check.py`, `19_clonal_frame_pairs.py` | does Gubbins flag the region; clonal-frame nearest-relative test of haplotype vs amikacin | `gubbins/region_block_coverage_by_isolate.csv`, `gubbins/amikacin_pairs_table.csv` |
| `20_gubbins_stratified.py` | tree-aware re-analysis using the Gubbins recombination-filtered tree (replaces the raw-distance clades of 13-14) | `gubbins/filtered_tree_distance_matrix.csv` (patristic distances of the Gubbins final tree, written here), `gubbins/mantel_mic_vs_filtered_tree.csv`, `gubbins/main_cluster_gene_tests.csv` |
| `21_segment_boundaries.py` | window-by-window identity of each carrier's own assembly to the reference vs public persicum genomes; boundaries in reference coordinates | `boundaries/window_calls.csv` |
| `22_reclassify_segment.py`, `23_segment_vs_amikacin.py` | reclassify all 82 kansasii genomes on the true extent (diagnostic sites); test classes against amikacin (MGIT, MIC, tree-aware) | `boundaries/segment_classes_v2.csv`, `boundaries/segment_genes.csv` |
| `24_rsmA_rsmI_alleles.py`, `25_segment_protein_divergence.py`, `26_rsm_alleles_vs_amikacin.py` | RsmA/RsmI alleles across cohort and public genomes; divergence ranking within the segment; alleles vs amikacin | `boundaries/RsmA_allele_per_isolate.csv`, `boundaries/RsmI_allele_per_isolate.csv`, `boundaries/segment_protein_identity.csv` |
| `12_mgit_lineage_followup.py` | characterise the top MGIT signal; ANI vs amikacin in the MIC cohort | `mgit/kansasii_mgit_block_ani_amikacin.csv`, `mgit/mic_cohort_ani_vs_mic.csv` |

Reference files for numbering (UniProt H37Rv proteins, E. coli rrnA region) are in `output/gwas/ref/` (made by `00_fetch_references.sh`).

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

## What changed with the larger cohort
| | 2026-10-06 run | 2026-10-16 run |
|---|---|---|
| genomes / usable with MIC | 131 / 30 | 186 / 80 |
| kansasii genomes in lineage analyses | 82 | 122 |
| amikacin 1 mg/L MGIT, kansasii (R / S+I) | 11 / 6 | 16 / 6 |
| MGIT top block (36 genes), R carriers vs non-R carriers | 10/11 vs 0/6, q = 0.15 | 12/16 vs 0/6, p = 0.003, q = 0.57 |
| ANI vs amikacin MIC, MIC cohort | p = 0.054 | p = 0.84 (n = 59) |
| persicum-like haplotype vs amikacin MIC | ~880: p = 0.135 (n = 4 carriers) | divergent class: identical medians, p = 0.84 (n = 59) |

The ANI lead and the independent-cohort support for the segment both disappear on the larger MIC cohort. Treat the persicum-segment
association as unconfirmed.

## Results
### Candidate-gene screen (script 03)
- 1,375 distinct features, 5,837 tests (kansasii 1,235; persicum 325; pooled species-centred 4,277). Within kansasii (59 isolates) the smallest
  FDR-adjusted q is 0.10: six folP variants carried by the same 2 isolates, rifabutin MIC lower by about 1 log2 (a folate-pathway gene is not
  a rifabutin target, so this is most likely chance with n = 2). Gene-drug pairs with a known mechanism: rpsL R43K vs streptomycin p = 0.014,
  q = 0.38; rplC T67S vs linezolid q = 0.46. Within persicum (12 isolates) the smallest q is 1.0. Nothing is significant.
- The pooled analysis with species-centred MIC again produces small q-values (about 0.0015) that are an artifact: rifabutin barely varies, and
  variants fixed in a species line up with species. Do not use it.
- Many features are carried by exactly the same isolates, so the real number of independent tests is much smaller than the count.

### Known hotspots (scripts 04-07)
- 50,664 variants were called in the candidate genes over 186 isolates. rpsL R43K is in 178 isolates (the reference carries R43).
- No resistance-associated change at gyrA A90/D94 (and 88, 89, 91), gyrB N499, rpoB H445/S450, embB M306/G406/Q497, katG S315, rpsL K88,
  rrl A2058/A2059/G2057 or rrs A1408/A1409/G1491 in any genome with coverage. The only differences at hotspots are synonymous codon changes
  (gyrB D461 GAT-GAC in 3, E501 in 2, rpoB D435, L452, Q432, S450 in 1-3 each; rplC C154 TGT-TGC in 35 plus 3 ambiguous; rpsL K88 in 2, one of them ambiguous).
- gyrA was read from the assembly annotation for 182 of 186 genomes: wild-type at all QRDR residues in every species. Not annotated:
  NR 10, 21, 113 (excluded) and 55.
- Read coverage at the hotspots (`hotspot_coverage_states.csv`): gyrA QRDR positions are low-coverage or masked in about 50 genomes (the
  snippy mask in persicum, pseudokansasii and attenuatum, see pitfall 4); rrl sites have no coverage in 34 genomes and rrs in 24. These are
  not called wild-type.

### MIC outlier tails (script 08)
- Outlier isolates analysed: Mkan329-054 (persicum), -055 (persicum), -047 and -050 (kansasii), -016 (pseudokansasii), -009 (persicum) and -008
  (innocens, only one innocens genome, so everything is "private").
- Mkan329-054 remains the multi-drug high outlier (above the cohort median for all 13 drugs; 773 private variants, 270 protein-altering).
  A uniform rise across drug classes suggests a general mechanism (efflux, permeability, regulation) or a testing artifact; the data cannot
  separate them. Mkan329-055 (persicum) is low for 4 drugs, with 62 private variants (23 protein-altering).
- Persicum isolates carry long private stretches (054: 773, 009: 932) because they sit on separate branches, so counts and gene-class
  keyword searches do not discriminate. None of the outlier isolates has MGIT data.

### Mkan329-054 and ethA (scripts 09-10)
- **Retracted earlier:** the apparent private ethA changes in 054 (E492D; A497T/V496I/P495A) are in `FCNCLM_03714`, a paralog, and are a
  mapping artifact (054's assembled protein for that locus is identical to every other persicum).
- The real EthA (`FCNCLM_02843`, 87% identical to Mtb EthA) has the same 11 differences from the kansasii reference in all 29 persicum genomes,
  including 054; no coding change in ethA or ethR is specific to 054. Its 44 calls against the reference are persicum-vs-kansasii divergence.
- One candidate remains: genome position 3145420 (A to C), 4 bp upstream of the ethA start. Among the isolates with ethionamide MIC only 054
  (10 mg/L) and 009 (0.6 mg/L) carry it (unadjusted p = 0.002, all species); an A to T change at the same position occurs in one innocens
  isolate with a low MIC. Caveats: two isolates, both persicum; the site was chosen after looking at the outliers; the carriers differ
  16-fold in MIC; it does not explain 054's other drugs. A hypothesis only.

### MGIT whole-gene analysis (scripts 11-12)
- 66 isolates have MGIT results and a genome. Unreadable results (`U`) are dropped.
- Resistance is rare. Tests with enough spread (3 or more R and 3 or more non-R in a species): amikacin 1 mg/L in kansasii (R 16, S+I 6) and
  amikacin 4 mg/L in kansasii (R 3, S+I 47 or S 46). 8,347 gene-level binary calls (protein-altering variant, rRNA variant, or any change in
  the 150 bp upstream region); Fisher exact test, BH correction over distinct carrier patterns (406 patterns at 1 mg/L, 1,115 at 4 mg/L).
- Result: nothing significant. Amikacin 1 mg/L: smallest q = 0.57; the top patterns are FCNCLM_02264 (0/16 R vs 4/6 S+I carriers, p = 0.002)
  and the block of 36 genes (carried by 12 of 16 resistant and 0 of 6 non-resistant, p = 0.003). Amikacin 4 mg/L: smallest q = 0.54, led by
  FCNCLM_04073 (3/3 R carriers vs 3 of 47 others, p = 0.001). rrs: no association (p = 0.53 at 1 mg/L, 0.37-0.38 at 4 mg/L).
- The 36-gene block is a divergence-from-reference signal, not a gene: block carriers have ANI 98.9-99.7% to the reference, non-carriers
  99.1-100%; carriers differ from each other by about 25,800 sites, as much as from non-carriers, so they are not one clone. ANI separates amikacin
  R from S+I only weakly (Mann-Whitney p = 0.060). The kansasii here is heterogeneous (122 genomes, 73 clusters at 1,000 differing sites).
  (The GTDB species 'kansasii' is the former subtype I only; the other subtypes are separate species.)
- Independent check in the MIC cohort (59 kansasii, 15 also in MGIT): amikacin median 1.0 log2 on both sides of the ANI threshold (p = 0.84);
  Spearman ANI vs MIC is positive for streptomycin (rho = 0.30, p = 0.019, unadjusted over 13 drugs) and not significant for amikacin
  (rho = 0.09). The earlier ANI lead (p = 0.054) is not reproduced.
- Caveat on the phenotype: 1 mg/L is a low MGIT concentration near the MIC of wild-type isolates, so "R" at 1 mg/L reflects a slightly higher
  MIC distribution, not clinical resistance.

### Lineages within kansasii (scripts 13-15) - raw-distance clades superseded
- 122 kansasii genomes; 3.95 Mb of core sequence (sites with an unambiguous base in all), 222,644 variable sites. Median pairwise distance
  20,948 core SNPs (max 56,386); median nearest-neighbour distance 32. The tree is neighbour-joining/average linkage on SNP counts without
  recombination filtering, so it is an approximation.
- The "clade A / clade B" split of script 14 was tuned to the 82-genome cohort; on 122 genomes it degenerates (A = 2, B = 118, 2 singletons)
  and its stratified tests (amikacin by clade, MIC by clade) are not meaningful. Script 14 is kept for the record; use the Gubbins-tree
  section below instead. Within "clade B" (n = 118) a block of genes around FCNCLM_01713-01717 (sigJ, rsmI, pMT1) is carried by 12 of 15 R
  and 0 of 5 S+I (p = 0.004, q = 0.88).
- Script 15 classes the 26 kb region FCNCLM_01713-01736 (NC_022663.1:1,944,065-1,970,050) by variant count against the reference. Counts fall in
  clean groups: none (0-30 variants), H~211 (179-215), H~400 mix (292-405) and H~880 (820-917). Over the 122 genomes the classes are 46
  none, 36 H~211, 6 mix, 34 H~880 (one isolate with 292 variants was placed in the mix class by the cut point). Pairwise Jaccard among
  all divergent haplotype isolates is 0.53 (the divergent class combines H~211, mix and H~880).
- Amikacin 1 mg/L MGIT (22 kansasii: 16 R): divergent haplotype R 13/14 vs ref-like R 3/8 (Fisher p = 0.011). In the independent MIC cohort
  (59 kansasii) amikacin MIC is identical in both groups (median 1.0 log2, p = 0.84), and no drug differs (smallest p = 0.075,
  sulfamethoxazole/trimethoprim).

### Source of the haplotype and Gubbins (scripts 16-19)
- BLAST against 72 public complex genomes (script 17): the H~880 haplotype (Mkan329-101, 25,529 bp) is 99.89% identical to M. persicum over its
  whole length (all 12 genomes, 99.81-99.93%) and 94.1% to kansasii on average. The reference-like haplotype is 93.5% to persicum, so the segment
  is persicum-derived DNA inside kansasii (interspecies transfer). The H~211 query (Mkan329-005) is a kansasii allele (98.8% within kansasii,
  94.6% to persicum).
- Gubbins final run (122 kansasii + reference, FastTree/GTRGAMMA): 1,543 recombination blocks (median 4,497 bp; 19.0 Mb summed over branches).
  Recombination is extensive, so raw SNP distances and the raw-distance clades of scripts 13-14 are recombination-inflated.
  - Fraction of the region covered by blocks: H~880 1.00 (0.96-1.00), H~400 mix 0.41, H~211 0.23, none 0.02 (0-0.28). 26 blocks overlap the
    region, 2 on branches shared by 20 or more isolates, 17 in single isolates. The region is flagged in nearly all H~880 isolates and in few
    others, so Gubbins now marks the carriers, not just a general hotspot (the 82-genome run flagged it in nearly every isolate).
  - The MRCA of the 34 H~880 isolates spans 120 tips (86 non-carriers inside), so they are not one clade. Carrier pairs are closer on
    the filtered tree (median 85 SNPs, range 11-192) than non-carrier pairs (114).
- Clonal-frame test (script 19, 22 phenotyped kansasii, 10 carriers): R-R tree distance 95 SNPs, S/I-S/I 128, R-S/I 118, so the phenotype has
  lineage signal. 17 nearest-relative pairs, 5 discordant for H~880: 3 support (carrier R, relative S/I), 2 do not; sign test p = 1.0.

### Stratification with the Gubbins tree (script 20)
- The recombination-filtered tree is shallow and star-like. Excluding Mkan329-068 (a 28,959-SNP terminal branch judged a tree artifact), 121
  isolates: 79 clusters at 40 SNPs, 54 at 60 SNPs, 20 at 100 SNPs (largest 70). The distance matrix is computed from the final tree by the script.
- Lineage signal in phenotype (permutation Mantel-type test): MGIT amikacin 1 mg/L rho = 0.27, p = 0.029 (n = 22). MIC cohort (n = 59):
  streptomycin p = 0.006 (q = 0.075), doxycycline p = 0.030, moxifloxacin p = 0.043 (negative rho), none significant after BH over 13 drugs;
  amikacin p = 0.42.
- H~880 carriers are spread over many backgrounds (34 carriers in 23 clusters at 40 SNPs, 14 at 60) but at fine scale carriers cluster:
  94% have a carrier as nearest relative (random expectation about 28%). This could reflect spread within lineages or residual recombinant
  signal in the tree.
- Stratified haplotype vs amikacin 1 mg/L: no informative stratum at 40 or 60 SNPs; at 100 SNPs two (exact Fisher p = 1.0 and 0.46; pooled OR
  3.0, MH p = 0.42). Not significant.
- Gene-level within the main 70-isolate cluster (14 phenotyped: 12 R, 2 S/I): nothing significant (166 patterns, min q = 0.68).

### Extent of the persicum-derived segment (scripts 21-23)
- Boundaries (script 21, ten H~880 and four control assemblies scanned in 1 kb windows over 1.80-2.10 Mb) are unchanged by the larger cohort:
  full segment NC_022663.1:1,929,793-1,994,501 (64.7 kb, 59 genes FCNCLM_01702-01761, 16 in the left half and 43 in the right, including rsmA and
  rsmI, a TetR regulator and cell-wall genes). The H~211 "haplotype" is the left part only (about 1,929-1,950 kb); the old 26 kb region straddled
  that internal boundary.
- Reclassification of the 122 kansasii genomes with 2,926 diagnostic sites (840 left, 2,086 right): full 31, left-only 39, none 42,
  partial/mosaic 10. Persicum-like DNA in the left part is common (70 of 122). 29 H~211 isolates are "ambiguous" because they share a fixed left
  fraction of 0.86 (a shorter left variant); 7 more H~211 isolates are left-only.
- Amikacin 1 mg/L, MGIT (22 kansasii, 16 R): full vs rest R 8/9 vs 8/13 (p = 0.33); any persicum-like piece (full or left-only) vs none/partial
  R 12/13 vs 4/9 (Fisher p = 0.023, nominal; several class definitions were compared, so it does not survive correction); left-only vs rest
  p = 0.54. By class: full 8R/1S-I, left-only 4R/0, none 2R/4, partial 2R/1.
- MIC cohort (59 kansasii): amikacin median 1.0 in the full, left-only and none classes; full vs others p = 0.71, full or left-only vs others
  p = 1.0. Not replicated.
- Tree-aware: five nearest-relative pairs discordant for the full segment (1 supports, 1 goes the other way, sign test p = 1.0); seven for any
  piece (2 vs 1, p = 1.0). At the 100 SNP cut two informative strata: full segment exact p = 1.0 and 0.46; full-or-left-only p = 1.0 and 0.066.
- Net: the amikacin association is nominal in the MGIT data only and is not supported by the independent MIC cohort or the tree-aware
  checks.

### rsmA and rsmI in the segment (scripts 24-26)
RsmA (FCNCLM_01703, 317 aa; 84% identical to Mtb KsgA) and RsmI (FCNCLM_01716, 285 aa) lie in the left part of the segment. Alleles were taken from
each genome's own Bakta proteins and checked with tblastn in 12 public persicum and 42 public kansasii genomes.
- Persicum-like alleles keep the length (317 and 285 aa); no truncation or frameshift in the main alleles (RsmA alleleP2 and P5 are 3 aa shorter).
- RsmA main persicum-like allele (44 cohort genomes): 10 differences from the kansasii reference (96.8% identity; 7 of 10 in the C-terminal tail).
  RsmI main persicum-like allele (71 genomes): 5 differences (98.2%). Both are persicum-like in all 12 public persicum genomes and in 7 (RsmA) and
  11 (RsmI) of 42 public kansasii genomes.
- RsmI is persicum-like in all 31 full and 36 of 39 left-only isolates; RsmA in all full but only 11 of 39 left-only isolates.
- Divergence is ordinary for the segment: across its 59 proteins the median identity is 94.2% (quartiles 89.8-96.5%); RsmA is more conserved
  than 19% of them and RsmI than 7%. No sequence-level reason to expect a functional change.
- Amikacin 1 mg/L MGIT (22 kansasii): persicum-like RsmA R 11/13 vs 5/9 (p = 0.18); RsmI R 13/14 vs 3/8 (Fisher p = 0.011). MIC cohort
  (59 kansasii): RsmA median 1.0 vs 1.0 (p = 0.35); RsmI 1.0 vs 1.0 (p = 0.90). The RsmI MGIT result is essentially the "any persicum-like
  piece" association again, was chosen after seeing it and is not replicated in the MIC cohort.

## Interpretation so far
No known resistance mutation explains any phenotype in either cohort. Recurrent patterns are not single mutations: persicum isolate 054 is a
multi-drug high outlier, and amikacin growth at 1 mg/L in kansasii tracks a persicum-derived segment (64.7 kb, FCNCLM_01702-01761, about
99.9% identical to M. persicum) in a recombination hotspot carried on several backgrounds. With 122 kansasii genomes the association is
nominal in MGIT (p = 0.01-0.02) but absent in the independent MIC cohort (59 isolates) and not supported by tree-aware tests (all p >= 0.066).
Both need more isolates and repeat testing; neither is confirmed.

## Open items
1. `01_build_table.py` still treats NR 184 (M. avium) as usable; exclude non-complex species and rerun (a one-isolate change, but it enters every
   pooled analysis). Also hard-coded: Mkan329-068 in scripts 20 and 23, the carrier and control
   isolates in script 21, and the clade labels of script 14, which no longer separate anything on 122 genomes.
2. Retest 054 and 055 (all drugs): a repeat decides whether the multi-drug shift is real.
3. Test whether position 3145420 changes ethA expression or ethionamide MIC (other persicum isolates, reporter assay).
4. Look at permeability and regulator genes in 054 beyond the filtered variant lists.
5. Calibrate MGIT 1 mg/L against MIC: the MIC-cohort amikacin median is 2 mg/L (log2 1.0), yet 6 of 22 MGIT kansasii isolates are
   S/I at 1 mg/L, so the two methods may not be on a comparable scale. MGIT cannot validate the MIC outliers (none has MGIT data).
6. Persicum needs more MIC-tested isolates; per-species analysis is limited to kansasii.
7. The three SXT-resistant MIC isolates (NR 168, 179, 182) belong to the jc batch (NR 134-186), which now has genomes; they have not been
   examined genotypically here.

## Windows checkout

- `.gitattributes` forces LF line endings, so a Windows checkout (even with
  `core.autocrlf=true`) keeps scripts runnable. Recommended: `git config core.autocrlf false`
  and `git config core.longpaths true`.
- Scripts default to the cluster root `/shares/sander.imm.uzh/MM/kansasii`; set the
  `KANSASII_ROOT` environment variable to point elsewhere (e.g. a mapped drive).
- Nextflow, Singularity and sbatch steps only run on the cluster (or WSL).


## Provenance

The sbatch files write `MANIFEST.txt` (commit, clean/dirty flag, script, conda env hash, inputs) into the output
directory when they finish; see `kansasii-docs/docs/PROVENANCE.md`. Project index: `kansasii-docs/README.md`.
