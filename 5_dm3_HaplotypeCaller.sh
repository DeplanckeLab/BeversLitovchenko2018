#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 4_HC.%I.err
#BSUB -o 4_HC.%I.out
#BSUB -J 4_HC[1-243]
#BSUB -M 25000000
#BSUB -R rusage[mem=25000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 5: PERFORM CALLING OF VARIANTS WITH HAPLOTYPE CALLER, CREATE GVCF
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add Development/java_jdk/1.8.0_66;
module add UHTS/Analysis/picard-tools/2.2.1;

GATKjar=/scratch/el/monthly/mlitovch/tools/GenomeAnalysisTK.jar
ref=/scratch/el/monthly/mlitovch/RefGen/dm3/dm3.Wolb.fa

deduplDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/3_dedupl_dm3/
targetDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/4_dm3_HC_a_target/
mkdir $targetDir
replaceRGdir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/4_dm3_HC_b_rgReplaced/
mkdir $replaceRGdir
realignDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/4_dm3_HC_c_realign/
mkdir $realignDir
vcfDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/4_dm3_HC_d_rawGVCF/
mkdir $vcfDir

# array of samples
zeroArr=( zero )
IDsArr=($(ls -d -1 $deduplDir/* | grep bam$ | sed 's@.*/@@' | cut -f1 -d'.'))
samples=("${zeroArr[@]}" "${IDsArr[@]}")

# sample
sample=${samples[${LSB_JOBINDEX}]};

# Local realignment around indels
java -jar $GATKjar -T RealignerTargetCreator \
                   -R $ref -I $deduplDir$sample.align.srt.dedupl.bam \
                   -o $targetDir$sample.intervals

# Replace read group
java -jar /software/UHTS/Analysis/picard-tools/2.2.1/bin/picard.jar \
          AddOrReplaceReadGroups I=$deduplDir$sample.align.srt.dedupl.bam \
          O=$replaceRGdir$sample.RGrepl.dedupl.bam SORT_ORDER=coordinate \
          RGID=$sample RGLB=$sample PL=illumina RGPU=illumina RGSM=Sample1 \
          CREATE_INDEX=True

# Perform realignment 
java -jar $GATKjar -T IndelRealigner -R $ref \
          -I $replaceRGdir$sample.RGrepl.dedupl.bam \
          -targetIntervals $targetDir$sample.intervals \
          -o $realignDir$sample.realign.dedupl.bam

# HaplotypeCaller - create GVCF
java -jar $GATKjar -T HaplotypeCaller -nct 8 -R $ref \
     -I $realignDir$sample.realign.dedupl.bam \
     -stand_call_conf 30 -stand_emit_conf 10 \
     --emitRefConfidence GVCF -o $vcfDir$sample.gatk.var.raw.g.vcf

# add sample name, important for merging
toReplace=Sample1
sed -i "s@$toReplace@$sample@g" $vcfDir$sample.gatk.var.raw.g.vcf

# Select no indels
java -jar $GATKjar -T SelectVariants -R $ref \
     --variant $vcfDir$sample.gatk.var.raw.g.vcf \
     -o $vcfDir$sample.gatk.var.raw.noIndelMNP.g.vcf \
     --selectTypeToExclude INDEL --selectTypeToExclude MNP -select "DP > 5"

exit 0;
