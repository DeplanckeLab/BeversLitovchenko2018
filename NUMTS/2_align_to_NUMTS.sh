#BSUB -L /bin/bash
#BSUB -e align_to_NUMTS.%I.err
#BSUB -o align_to_NUMTS.%I.out
#BSUB -J align_to_NUMTS[1-8]
#BSUB -M 36000000
#BSUB -R rusage[mem=36000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 2: ALIGN DATA TO THE PREVIOUSLY CREATED NUMTs GENOME
#------------------------------------------------------------------------------

export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Aligner/bwa/0.7.13;
module add UHTS/Analysis/BEDTools/2.26.0;

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
  groupInfo="@RG\tID:"$3"\tSM:"$4"\tPL:illumina\tLB:"$5"\tPU:unit1"
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
  groupInfo="@RG\tID:"$4"\tSM:"$5"\tPL:illumina\tLB:"$6"\tPU:unit1"
  alignFileName=$7$4"_"$5"_"$6".align.sam"
  alignBAMname=$7$4"_"$5"_"$6".align.bam"
  bwa mem -M -R $groupInfo $3 $1 $2 > $alignFileName
  samtools view -bS $alignFileName > $alignBAMname
  samtools sort $alignBAMname -o $7$4"_"$5"_"$6".align.srt.bam"
  samtools index $7$4"_"$5"_"$6".align.srt.bam"
  samtools flagstat $7$4"_"$5"_"$6".align.srt.bam" > $7$4"_"$5"_"$6".align.srt.flagstat"
  #rm $alignFileName
}

# path to raw files which will be aligned
trimmedDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/trimmed/
mergeDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/merged_2018/
mkdir $mergeDir
bwaDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/mapped_numts_2018/
mkdir $bwaDir
deduplDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/mapped_numts_dedupl_2018/
mkdir $deduplDir
clearMatchDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/clearMatch_2018/
mkdir $clearMatchDir
clearMatchFastqDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/clearMatchFastq_2018/
mkdir $clearMatchFastqDir
clearMatchMitoDir=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/clearMatchMito_2018/
mkdir $clearMatchMitoDir

mouseNumts=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/mm10_numts.fa
flyNumts=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/dm6_numts.fa
mouseMT=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/mm10_chrMT.fa
flyMT=/scratch/el/monthly/mlitovch/MitoSeq_Ref_DMel_MMuc/dm6_chrM.fa

IDs=( zero 
BRSX4
BRX1
BRX2
BRX3
MA1E
MA2F
MB1G
MB2H
)

species=( zero
fly
fly
fly
fly
mouse
mouse
mouse
mouse
)

echo "STEP 1: map to numts genome"
# current file to work with
ID=${IDs[${LSB_JOBINDEX}]};
specie=${species[${LSB_JOBINDEX}]};
R1s=($(ls -d -1 $trimmedDir/* | grep $ID"_" | grep _R1_));
for ((i=0; i < ${#R1s[@]}; i++))
do
  info=$specie
  R1=${R1s[$i]}
  lib=`echo $R1 | sed 's@.*_00@00@g' | sed 's@.trimmed.fastq@@g'`
  # matching R2
  R2=`ls -d -1 $trimmedDir/* | grep $ID"_" | grep _R2_ | grep $lib`
 
  echo $R1
  echo $R2
  echo $lib
  echo $info 
 
  if [ "$specie" == "mouse" ]; then 
    echo "mouse"
   # mapWithBWA_PE $R1 $R2 $mouseNumts $ID $info $lib $bwaDir
  else
    echo "fly"
   # mapWithBWA_PE $R1 $R2 $flyNumts $ID $info $lib $bwaDir
  fi
  #rm $bwaDir$ID"_"$info"_"$lib".align.sam"
  #rm $bwaDir$ID"_"$info"_"$lib".align.bam"
done

echo "STEP 2: merge mapped files"
lib=numts
toMerge=`ls -d -1 $bwaDir/* | grep $ID"_" | grep srt.bam$ | xargs`
#samtools merge $bwaDir$ID"_"$info"_"$lib".align.srt.bam" $toMerge

echo "STEP 3: mark duplicates"
java -jar /scratch/el/monthly/mlitovch/tools/picard/picard.jar \
          MarkDuplicates I=$bwaDir$ID"_"$info"_"$lib".align.srt.bam" \
          O=$deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" \
          CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
          M=$deduplDir$ID"_"$info"_"$lib".metrix" \
          REMOVE_DUPLICATES=TRUE

echo "STEP 4: extract reads with perfect match to NUMTS"
#  Reads with perfect match to Numts were extracted out 
samtools view -H $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" > $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.sam"
samtools view -h $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" | awk '$6 ~ /^[0-9]+M$/' >> $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.sam"
samtools view -bS $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.sam" > $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.bam"
rm $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.sam"
bedtools bamtofastq -i $clearMatchDir$ID"_"$info"_"$lib".align.srt.dedupl.clearMatch.bam" \
                    -fq $clearMatchFastqDir$ID"_"$info"_"$lib".mapped.clearMatch.fastq" 

# and aligned to the mitochondrial genome
R1=$clearMatchFastqDir$ID"_"$info"_"$lib".mapped.clearMatch.fastq"
lib=numts_to_mito
if [ "$specie" == "mouse" ]; then
  echo "mouse"
  mapWithBWA_SE $R1 $mouseMT $ID $info $lib $clearMatchMitoDir
else
  echo "fly"
  mapWithBWA_SE $R1 $flyMT $ID $info $lib $clearMatchMitoDir
fi
rm $clearMatchMitoDir$ID"_"$info"_"$lib".align.sam"
rm $clearMatchMitoDir$ID"_"$info"_"$lib".align.bam"

# After filtering out reads with a perfect match to the mitochondrial genome, 
# the remaining reads were considered as contaminations from Numts.
fromMito=`samtools view -h $clearMatchMitoDir$ID"_"$info"_"$lib".align.srt.bam" | awk '$6 ~ /^[0-9]+M$/' | wc -l`
total=`samtools view -h $clearMatchMitoDir$ID"_"$info"_"$lib".align.srt.bam" | wc -l`
echo $ID" "$fromMito" "$total > $ID".numts.perc.txt"

exit 0;
