# BeversLitovchenko2018
Source code to support the paper: "Extensive  mitochondrial  population  structure  and  haplotype-specific  variation  in  metabolic  phenotypes  in  the  Drosophila  Genetic  Reference  Panel"

| File  | Description |
| ------------- | ------------- |
| 1_QC.sh   | step 1: perform qc check of the raw files |
| 2_trim.sh   | step 2: perform trimming of the adapters on the raw files |
| 3_dm3_map_dedupl.sh   | step 3: perform mapping to dm3 (for genotyping)  |
| 4_dm3_getMap_stats.sh   | step 4: check the quality of the maping on dm3  |
| 5_dm3_HaplotypeCaller.sh   | step 5: perform calling of variants with haplotype caller, create gvcf  |
| 6_dm3_GenotypeGVCF.sh   | step 6: perform merging of raw gvcf with genotypegvcfs  |
| 7_dm3_GenotypeWithR.sh   | step 7: perform genotyping with use of custom r script  |
| 7a_dm3_GenotypeGVCF.R   |  step 7: custom r script to perform genotyping |
| 8_dm6_map_dedupl_dm6.sh   |  step 8: perform mapping to dm6 & simultaniously correct genotypes |
| 9_dm6_getMap_stats.sh   |  step 9: check the quality of the maping on dm6 |
| 10_dm6_merge_BeversData.sh   | step 10: merge mitochandrial gvcfs into one  |
| 11_dm6_5950_5975_repeats.sh   |  step 11: extract reads mapping to the intergenic repeat region of the mitochondrial genome |
| 11A_Heteroplasmy.py   | step 12: perform calculations of heteroplasmy |
| 12_Heteroplasmy.sh   | step 12: perform calculations of heteroplasmy   |
| 13_parse_nuclVars_for_GRD.sh   | step 13: move dgrp2 vcf from dm3 to dm6, select variants by maf  |
| 14_GRDs.sh   | step 14: detect GRDs  |
| 14a_GRDs.R   | step 14: custom r script to detect GRDs  |
| 15_annotate_GRD.sh   |  step 15: perform annotation of GRDs with use of snpeff |
| 16_functions.R   | source file for the functions used in r scripts  |
| 17_baseAnalysis.R   |  base analysis of the variation in mitochondrial genomes of dgrps |
| 18_compToRich_richGenome.R   | comparison of the variants from this study to the previously published one  |
| 19_heteroplasmy.R   |  post-processing of the heteroplasmy detection results |
| 20_GRD_postProc.R   |  post-processing of the GRDs |
| 21_parse_PhenoDB_MitoVcf_forGWAS.R   |  parsing of the data base of the phenotypes for further use in GWAS |
| 22_GWAS_Haplogroups.R   | GWAS on mitochondrial haplogroups  |
| 23_GRD_phenotypes.R   | association of GRDs with phenotypes  |
| Comparison to Richardson   |  folder containing scripts used for the comparison to the Richardson study |
| MitoRepeats   | java source code used to process reads falling into repetitive region of mitochondrial genome |
| NUMTS   | folder containing scripts used for the detection and quantification of NUMTs |
