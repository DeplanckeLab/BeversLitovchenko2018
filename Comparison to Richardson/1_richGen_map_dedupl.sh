#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 13_BWA_rich.%I.err
#BSUB -o 13_BWA_rich.%I.out
#BSUB -J 13_BWA_rich[1-243]
#BSUB -M 32000000
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 14: PERFORM MAPPING TO THE GENOME USED IN RICHARDSON PAPER
# (chrU:5288528–5305749), Wolbachia (GenBank ID: AE017196); and an
# equivalently-sized (1.2 Mb) nuclear random region 3L (GenBank ID: NT037436;
# positions 10000000–11200000))
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Aligner/bwa/0.7.13;
module add UHTS/Analysis/picard-tools/2.2.1;

#===  FUNCTION  ===============================================================
# NAME:  mapWithBWA_SE
# DESCRIPTION: maps single-end (SE) reads with BWA
# PARAMETER  1: path to first read pair file
# PARAMETER  2: path to reference genome
# PARAMETER  3: ID of the sequencing, i.e. DGRP_100
# PARAMETER  4: sample info, i.e. A1
# PARAMETER  5: library name, i.e lib1 or run1
# PARAMETER  6: folder to put files
#==============================================================================
mapWithBWA_SE() {
  groupInfo="@RG\tID:"$3"\tSM:"$3"_"$4"\tPL:illumina\tLB:"$5"\tPU:unit1"
  alignFileName=$6$3"_"$4"_"$5".align.sam"
  alignBAMname=$6$3"_"$4"_"$5".align.bam"
  bwa mem -M -R $groupInfo $2 $1 > $alignFileName
  samtools view -bS $alignFileName > $alignBAMname
  samtools sort $alignBAMname -o $6$3"_"$4"_"$5".align.srt.bam"
  samtools index $6$3"_"$4"_"$5".align.srt.bam"
  samtools flagstat $6$3"_"$4"_"$5".align.srt.bam" > $6$3"_"$4"_"$5".align.srt.flagstat"
  #rm $alignFileName
}

#===  FUNCTION  ===============================================================
# NAME:  mapWithBWA_PE
# DESCRIPTION: maps paired-end (PE) reads with BWA
# PARAMETER  1: path to first read pair file
# PARAMETER  2: path to second read pair file
# PARAMETER  3: path to reference genome
# PARAMETER  4: ID of the sequencing
# PARAMETER  5: sample info
# PARAMETER  6: library name
# PARAMETER  7: folder to put files
#==============================================================================
mapWithBWA_PE() {    
  groupInfo="@RG\tID:"$4"\tSM:"$4"_"$5"\tPL:illumina\tLB:"$6"\tPU:unit1"
  alignFileName=$7$4"_"$5"_"$6".align.sam"
  alignBAMname=$7$4"_"$5"_"$6".align.bam"
  bwa mem -M -R $groupInfo $3 $1 $2 > $alignFileName
  samtools view -bS $alignFileName > $alignBAMname
  samtools sort $alignBAMname -o $7$4"_"$5"_"$6".align.srt.bam"
  samtools index $7$4"_"$5"_"$6".align.srt.bam"
  samtools flagstat $7$4"_"$5"_"$6".align.srt.bam" > $7$4"_"$5"_"$6".align.srt.flagstat"
  #rm $alignFileName
}


# path to trimmed files which will be aligned
trimDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/1_trimmed/
bwaDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/2_mapped_richardson/
mkdir $bwaDir
deduplDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/3_dedupl_richardson/
mkdir $deduplDir
refGen=/scratch/el/monthly/mlitovch/RefGen/richardson/richardson.fa

zeroArr=( zero )

# list all R1s
R1Arr=($(ls -d -1 $trimDir/* | grep gz$ | grep _R1_))
firstInPair=("${zeroArr[@]}" "${R1Arr[@]}")

# IDs
IDsArr=($(ls -d -1 $trimDir/* | grep gz$ | grep _R1_ | sed 's@.*/@@' | cut -f1 -d'_'))
IDs=("${zeroArr[@]}" "${IDsArr[@]}")

# info
infoArr=($(ls -d -1 $trimDir/* | grep gz$ | grep _R1_ | sed 's@.*/@@' | cut -f2 -d'_'))
infos=("${zeroArr[@]}" "${infoArr[@]}")

R1=${firstInPair[${LSB_JOBINDEX}]};
# matching R2, if any
R2=`ls -d -1 $trimDir/* | grep $(echo $R1 | sed 's@.*/@@' | cut -f1,2 -d "_") | grep _R2_`;
ID=${IDs[${LSB_JOBINDEX}]};
info=${infos[${LSB_JOBINDEX}]};
# I use genome as library
lib=richardson

# check, whatever or not sequencing was pair-end and align accordingly
if [ "$R2" == "" ]; then 
  mapWithBWA_SE $R1 $refGen $ID $info $lib $bwaDir
else
  mapWithBWA_PE $R1 $R2 $refGen $ID $info $lib $bwaDir
fi

java -jar /software/UHTS/Analysis/picard-tools/2.2.1/bin/picard.jar \
          MarkDuplicates I=$bwaDir$ID"_"$info"_"$lib".align.srt.bam" \
          O=$deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" \
          CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
          M=$deduplDir$ID"_"$info"_"$lib".metrix" \
          REMOVE_DUPLICATES=TRUE
samtools flagstat $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" > $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.flagstat"
samtools idxstats $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" | cut -f 1,3 > $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.idxstat"

exit 0;
