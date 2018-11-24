#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 32_GRD_DGRP2.%I.err
#BSUB -o 32_GRD_DGRP2.%I.out
#BSUB -J 32_GRD_DGRP2[1-72]
#BSUB -M 32000000
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 14: DETECT GRDs
#------------------------------------------------------------------------------

export PATH=/software/bin:$PATH;
module use /software/module/;

chroms=( chr2L chr3L chr2R chr3R chrX chr4 )
# there is 33 mito variants which pass cutoff of 150 genotyped lines and MAF > 0.05
chroms=($(for x in "${chroms[@]}"; do printf "$x%.0s " {1..12}; done))
mitoVars=($(seq 12))
mitoVars=( "${mitoVars[@]}" "${mitoVars[@]}" "${mitoVars[@]}" "${mitoVars[@]}" "${mitoVars[@]}" "${mitoVars[@]}" )

zeroArr=( zero )
chroms=("${zeroArr[@]}" "${chroms[@]}")
mitoVars=("${zeroArr[@]}" "${mitoVars[@]}")

chr=${chroms[${LSB_JOBINDEX}]};
mito=${mitoVars[${LSB_JOBINDEX}]};

baseDir=/scratch/el/monthly/mlitovch/GRD_DGRP/
mitoVCF=$baseDir"Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf"
chromVCF=$baseDir"dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"

Rscript --vanilla 15a_GRDs.R $mitoVCF $mito $chromVCF $baseDir"GRD_DGRP2_"$chr"_"$mito".Rds"

exit 0;
