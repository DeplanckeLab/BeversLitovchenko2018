#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 0_QC.%I.err
#BSUB -o 0_QC.%I.out
#BSUB -J 0_QC[1-405]
#BSUB -M 8000000
#BSUB -R rusage[mem=8000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Quality_control/fastqc/0.11.2;
module add UHTS/Analysis/samtools/1.2;

#------------------------------------------------------------------------------
# STEP 1: PERFORM QC CHECK OF THE RAW FILES
#------------------------------------------------------------------------------
# list all the fastqs
fastqFolder=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/0_rawData/
zeroArr=( zero )
samplesArr=($(ls -d -1 $fastqFolder*/* | grep gz$))
samples=("${zeroArr[@]}" "${samplesArr[@]}")

# do QC
# make new directory
mkdir /scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/0_rawData_QC
cd /scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/0_rawData_QC
sample=${samples[${LSB_JOBINDEX}]};
fastqc $sample

exit 0;
