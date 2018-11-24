#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 14_heteroplasmy.%I.err
#BSUB -o 14_heteroplasmy.%I.out
#BSUB -J 14_heteroplasm[1-179]
#BSUB -M 16000000
#BSUB -R rusage[mem=16000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 12: PERFORM CALCULATIONS OF HETEROPLASMY
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;

baseDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/

zeroArr=( zero )
allBams=($(ls -d -1 $baseDir"8_chrM_dm6/"* | grep bam$ ))
bams=("${zeroArr[@]}" "${allBams[@]}")

bam=${bams[${LSB_JOBINDEX}]};
bai=$bam".bai"
ID=`echo $bam | sed 's@.*dm6/@@g' | sed 's@_BWA_dm6.align.srt.dedupl.chrM.bam@@g'`
echo $ID
mkdir $baseDir$ID"/"
cp $bam $baseDir$ID"/"
cp $bai $baseDir$ID"/"
cp $baseDir"8_chrM_dm6/13_Heteroplasmy.py" $baseDir$ID"/"
cd $baseDir$ID"/"
python $baseDir$ID"/13A_Heteroplasmy.py"
mv *.het $baseDir

exit 0;

