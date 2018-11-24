#!/bin/bash

#------------------------------------------------------------------------------
# STEP 7: PERFORM GENOTYPING WITH USE OF CUSTOM R SCRIPT
#------------------------------------------------------------------------------

#===  FUNCTION  ===============================================================
# NAME:  waitall
# DESCRIPTION: pauses the system until all PIDs are finished
# PARAMETER  1: PIDs
# EXAMPLE: for sample in ${samples[@]}
#	   do
# 	     for fileR1 in $fastqR1
#            do
#	     done
#	     waitall $pids
#	    done
#==============================================================================
waitall() { # PID...
  ## Wait for children to exit and indicate whether all exited with 0 status.
  local errors=0
  while :; do
    #debug "Processes remaining: $*"
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
        #echo "$pid is still alive."
        set -- "$@" "$pid"
      elif wait "$pid"; then
        echo "$pid exited with zero exit status."
      else
        echo "$pid exited with non-zero exit status."
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    # TODO: how to interrupt this sleep when a child terminates?
    sleep ${WAITALL_DELAY:-1}
   done
  ((errors == 0))
}

#===  FUNCTION  ===============================================================
# NAME:  createGenotypeVCF
# DESCRIPTION: creates vcf files used further for genotyping from DGRP2 and 
#	       deplancke vcf
# PARAMETER  1: generic path to deplancke variants,
#               I'll need *.vcf.gz and *.GT.FORMAT
# PARAMETER  2: column (sample effectively) to use, ranges from >= 3
# PARAMETER  3: path to dgrp2 snp vcf
# PARAMETER  4: path to the final destination
# EXAMPLE: 
#==============================================================================
createGenotypeVCF() {
      echo $1
      echo $2
      sampleName=`head -1 $1".GT.FORMAT" | cut -f $2`;
      echo "Working with "$sampleName
      # get genotyped sites for this sample from deplancke vcf
      cut -f 1,2,$2 $1".GT.FORMAT" | awk '$3 !="./." && $1 !="Wolbachia" && $1 !="CHROM" {print $1 "\t" $2 }' > $sampleName".pos";

      # extract those sites from deplancke vcf
      deplCommand="vcftools --gzvcf "$1".vcf.gz --positions "$sampleName".pos --indv "$sampleName" --recode --recode-INFO-all --out "$sampleName"_GT_in_Depl";
      eval $deplCommand;
      bgzip -c $sampleName"_GT_in_Depl.recode.vcf" > $sampleName"_GT_in_Depl.vcf.gz";
      tabix -p vcf $sampleName"_GT_in_Depl.vcf.gz";
      rm $sampleName"_GT_in_Depl.recode.vcf";

      # extract matching sites from dgrp2
      dgrpCommand="vcftools --gzvcf "$3" --positions "$sampleName".pos --recode --recode-INFO-all --out "$sampleName"_GT_in_DGRP2";
      eval $dgrpCommand;
      bgzip -c $sampleName"_GT_in_DGRP2.recode.vcf" > $sampleName"_GT_in_DGRP2.vcf.gz";
      tabix -p vcf $sampleName"_GT_in_DGRP2.vcf.gz";
      rm $sampleName"_GT_in_DGRP2.recode.vcf";
           
      # merge them together. Although I do it here, I'm not sure that original alleles will be preserved
      vcf-merge $sampleName"_GT_in_Depl.vcf.gz" $sampleName"_GT_in_DGRP2.vcf.gz" | bgzip -c > $sampleName"_Depl_DGRP2.vcf.gz";

      # move to the result folder
      mv $sampleName"_GT_in_DGRP2.vcf.gz" $sampleName"_GT_in_Depl.vcf.gz" $sampleName"_Depl_DGRP2.vcf.gz" $4; 
}


#===  INPUTS  =================================================================
deplVars=~/Desktop/MitoSeq_RPJB_ML_Aug2017/MitoSeq_AllRuns_dm3.gvcf.snps.fltr
dgrpSNPsVCF=~/Desktop/dgrp2.SNPs.vcf.gz
numbOfProc=5
outputDir=~/Desktop/MitoSeq_RPJB_ML_Aug2017/results/genotyping/vcfForGenotyping/

# step 1: get GT information from deplancke variants
vcftools --gzvcf $deplVars".vcf.gz" --extract-FORMAT-info GT --out $deplVars

# step 2: get number of samples (actual # of samples is this - 2, but I need 
# index)
numbOfSampl=`awk '{print NF}' $deplVars".GT.FORMAT" | sort -nu | tail -n 1`

# step 3: get vcf
indices=($(seq 3 $numbOfProc $numbOfSampl))
printf '%s\n' "${indices[@]}"

for i in "${indices[@]}"
do
   pids=""
   numbOfProc1=$((numbOfProc - 1))
   indToAdd=($(seq 0 1 $numbOfProc1))
   for j in "${indToAdd[@]}"
   do
      createGenotypeVCF $deplVars $((i + j)) $dgrpSNPsVCF $outputDir &
   done
   pids="$pids $!"
   waitall $pids
done
