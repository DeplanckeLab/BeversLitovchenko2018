#!/bin/bash

#------------------------------------------------------------------------------
# STEP 11: EXTRACT READS MAPPING TO THE INTERGENIC REPEAT REGION OF THE
# MITOCHONDRIAL GENOME
#------------------------------------------------------------------------------

module add UHTS/Analysis/BEDTools/2.26.0;
echo -e "chrM\t5959\t5983" > chrM_5959_5983.bed
regionDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/7_chrM_5959_5983_bams
mkdir $regionDir
allFiles=($(ls -d -1 7_dedupl_dm6/* | grep align.srt.dedupl.bam$))

for ((i=0; i<=${#allFiles[@]}; i++)); do
    sample=${allFiles[i]}
    outputTo=`echo $sample | sed 's@7_dedupl_dm6/@@g' | sed 's/dm6/dm6_chrM_5959_5983/g'` 
    bedtools intersect -F 1 -abam $sample -b chrM_5959_5983.bed > $regionDir/$outputTo
    samtools index $regionDir/$outputTo
    txtOutput=`echo $outputTo | sed 's/bam/txt/g'`
    samtools view  $regionDir/$outputTo | cut -f 10 > $regionDir/$txtOutput
    finalOutput=`echo $txtOutput | sed 's/.txt/_stats.txt/g'`
    grep -o CTA.*GGG $regionDir/$txtOutput | sort | uniq -c | sort -nr | sed 's/  */ /g'  | sed 's/^ //g' > $regionDir/$finalOutput
done

rm chrM_5959_5983.bed
