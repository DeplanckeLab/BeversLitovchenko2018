#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 1_trim.%I.err
#BSUB -o 1_trim.%I.out
#BSUB -J 1_trim[1-243]
#BSUB -M 8000000
#BSUB -R rusage[mem=8000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Quality_control/fastqc/0.11.2;
module add UHTS/Analysis/samtools/1.2;
module add UHTS/Quality_control/cutadapt/1.8;

#------------------------------------------------------------------------------
# STEP 2: PERFORM TRIMMING OF THE ADAPTERS ON THE RAW FILES, 
# only paired-end runs need it (run1 and run2), but I'll trim run3 as well, 
# to be consistent
#------------------------------------------------------------------------------
fastqFolder=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/0_rawData/

# make new directory
putTrimmedIn=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/1_trimmed/
putFastqcIn=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/1_trimmed-QC/
mkdir $putTrimmedIn
mkdir $putFastqcIn
# path to trim galore
trimGalore=/scratch/el/monthly/mlitovch/tools/trim_galore

# list all R1s
zeroArr=( zero )
R1Arr=($(ls -d -1 $fastqFolder*/* | grep gz$ | grep _R1_))
firstInPair=("${zeroArr[@]}" "${R1Arr[@]}")

# list all R2s
zeroArr=( zero )
R2Arr=($(ls -d -1 $fastqFolder*/* | grep gz$ | grep _R2_))
secondInPair=("${zeroArr[@]}" "${R2Arr[@]}")

R1=${firstInPair[${LSB_JOBINDEX}]};
R2=${secondInPair[${LSB_JOBINDEX}]};

# check, whatever or not sequencing was pair-end and trim
if [ "$R2" == "" ]; then 
  $trimGalore $R1 --nextera -o $putTrimmedIn --fastqc
else
  $trimGalore --paired $R1 $R2 --nextera -o $putTrimmedIn --fastqc
fi

exit 0;
