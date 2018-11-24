#!/bin/bash

#------------------------------------------------------------------------------
# STEP 15: PERFORM ANNOTATION OF GRDs WITH USE OF SNPEFF
#------------------------------------------------------------------------------

snpEffJar=/home/litovche/bin/snpEff/snpEff.jar
baseDir=/home/litovche/Desktop/MitoSeq_RPJB_ML_Aug2017/
inputVCF=$baseDir"dgrp2_dm6_BeversLinesOnly_MAF005.vcf.gz"

outputDir=/home/litovche/Desktop/MitoSeq_RPJB_ML_Aug2017/results/GRD_DGRP2/
cd $outputDir

GRDsNuclPositions=$outputDir"GRDs_DGRP2_NuclVars.txt"

vcftools --gzvcf $inputVCF --positions-overlap $GRDsNuclPositions \
         --recode --recode-INFO-all --out $outputDir"GRDs_DGRP2_NuclVars"
mv $outputDir"GRDs_DGRP2_NuclVars.recode.vcf" $outputDir"GRDs_DGRP2_NuclVars.vcf"
java -jar $snpEffJar BDGP6.82 \
         -upDownStreamLen 2000 \
         $outputDir"GRDs_DGRP2_NuclVars.vcf" > $outputDir"GRDs_DGRP2_NuclVars.annot.vcf"

rm $outputDir"GRDs_DGRP2_NuclVars.vcf"
rm "*.log"
rm "snpEff_genes.txt" "snpEff_summary.html"
rm $outputDir"GRDs_DGRP2_NuclVars.txt"
echo "Finished annotation"


