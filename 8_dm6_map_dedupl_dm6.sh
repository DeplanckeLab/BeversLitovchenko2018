#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 9_BWA_dm6.%I.err
#BSUB -o 9_BWA_dm6.%I.out
#BSUB -J 9_BWA_dm6[1-179]
#BSUB -M 32000000
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 8: PERFORM MAPPING TO DM6 & SIMULTANIOUSLY CORRECT GENOTYPES
# ALSO CALL VARIANTS
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Aligner/bwa/0.7.13;
module add UHTS/Analysis/picard-tools/2.2.1;
module add UHTS/Analysis/vcftools/0.1.14;

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


# path to trimmed files which will be aligned
trimDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/1_trimmed/
bwaDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/6_mapped_dm6/
mkdir $bwaDir
deduplDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/7_dedupl_dm6/
mkdir $deduplDir
chrMDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/8_chrM_dm6/
mkdir $chrMDir
targetDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/9_dm6_HC_a_target/
mkdir $targetDir
replaceRGdir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/9_dm6_HC_b_rgReplaced/
mkdir $replaceRGdir
realignDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/9_dm6_HC_c_realign/
mkdir $realignDir
vcfDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/9_dm6_HC_d_rawVCF/
mkdir $vcfDir

GATKjar=/scratch/el/monthly/mlitovch/tools/GenomeAnalysisTK.jar
refGen=/scratch/el/monthly/mlitovch/RefGen/dm6/dm6.Wolb.fa

# list all files names
initialNames=( zero
Berk1_S71
Berk2_S98
21E9_L5
26D9_L5
28A4_L8
31D7_L5
32F9_L5
508C9_L5
41_S94
45_S57
48B12_L8
57B_RB
101_N712
69_S59
73E8_L5
75A5_L8
83B5_L8
85B9_L5
88_S72
91A10_L8
93A6_L8
100_S96
59_S91
105_S89
859_N712
129A2_L8
136B8_L8
138_S109
142_S56
149A8_L8
153A12_L5
362A3_L5
161B6_L8
181G9_L5
189E10_L8
195C6_L8
208_S121
217B1_L8
223B11_L8
228A9_L8
229C10_L8
235F1_L5
237D10_L8
239G6_L5
280_S120
287C3_L8
301E4_L5
303C4_L8
306A3_L8
307D1_L8
309F2_L8
313A7_L8
315A1_L8
317C1_L5
318D8_L8
319B5_L5
321B7_L8
324A5_L5
332_S108
336_S50
338E11_L8
737C5_L8
348F10_L8
352G7_L5
354_Ad2
355C9_L8
356F8_L5
357G12_L8
358_S44
359B3_L8
361E2_L8
362E4_L8
365F2_L5
367_Ad2
819C10_L5
371C3_L5
373G2_L8
375G8_L8
377G5_L8
379A7_L5
383G10_L8
227_N712
382A9_L5
383_S48
385G2_L5
386E10_L5
390F12_L8
392B4_L5
395C4_L5
397F3_L8
399F11_L8
405_S46
406E3_L8
409A1_L5
426G1_L8
427G3_L8
437G6_L8
440E7_L5
441B2_L5
443E5_L8
486B10_L5
491_S45
492H2_L8
502E5_L5
508_S81
509D10_L5
513B2_L8
517_S122
42C1_L8
530E1_L5
535F3_L5
551E6_L5
555E12_L8
559F4_L8
563_S74
566E6_L8
584E7_L8
596_S110
627E9_L8
630F1_L8
634D12_L8
639G11_L5
235_N712
703G4_L5
707B10_L8
712B4_L8
714_S117
716C12_L8
721B1_L5
730_S70
732C11_L8
340D9_L8
738D5_L8
360_S115
761_S83
765G10_L5
774E12_L5
776C2_L5
783E8_L8
786G5_L5
787A6_L5
790D12_L5
799_S100
801C7_L8
802A11_L5
808F5_L5
810D3_L5
812_S76
370B8_L5
820_S68
821B11_L5
822F4_L5
832_S65
837G7_L8
843_S53
849_S93
853H6_L8
855F8_L8
857_S77
109_S103
879F7_L5
882G8_L5
884_S52
887_S101
890H7_L8
892_S107
894_S80
897C8_L5
900G1_L5
907H4_L8
908H5_L8
911H11_L8
913H12_L8
ore_S64
ore2F12_L5
ore3D8_L5
w1E2_L5
w1118_N712
w2A2_L5
)

# list all new, real GT, names
realGTs=( zero
Berk1
Berk2
DGRP-021
DGRP-026
DGRP-028
DGRP-031
DGRP-032
DGRP-038
DGRP-041
DGRP-045
DGRP-048
DGRP-057
DGRP-059
DGRP-069
DGRP-073
DGRP-075
DGRP-083
DGRP-085
DGRP-088
DGRP-091
DGRP-093
DGRP-100
DGRP-101
DGRP-105
DGRP-109
DGRP-129
DGRP-136
DGRP-138
DGRP-142
DGRP-149
DGRP-153
DGRP-158
DGRP-161
DGRP-181
DGRP-189
DGRP-195
DGRP-208
DGRP-217
DGRP-223
DGRP-228
DGRP-229
DGRP-235
DGRP-237
DGRP-239
DGRP-280
DGRP-287
DGRP-301
DGRP-303
DGRP-306
DGRP-307
DGRP-309
DGRP-313
DGRP-315
DGRP-317
DGRP-318
DGRP-319
DGRP-321
DGRP-324
DGRP-332
DGRP-336
DGRP-338
DGRP-340
DGRP-348
DGRP-352
DGRP-354
DGRP-355
DGRP-356
DGRP-357
DGRP-358
DGRP-359
DGRP-361
DGRP-362
DGRP-365
DGRP-367
DGRP-370
DGRP-371
DGRP-373
DGRP-375
DGRP-377
DGRP-379
DGRP-380
DGRP-381
DGRP-382
DGRP-383
DGRP-385
DGRP-386
DGRP-390
DGRP-392
DGRP-395
DGRP-397
DGRP-399
DGRP-405
DGRP-406
DGRP-409
DGRP-426
DGRP-427
DGRP-437
DGRP-440
DGRP-441
DGRP-443
DGRP-486
DGRP-491
DGRP-492
DGRP-502
DGRP-508
DGRP-509
DGRP-513
DGRP-517
DGRP-528
DGRP-530
DGRP-535
DGRP-551
DGRP-555
DGRP-559
DGRP-563
DGRP-566
DGRP-584
DGRP-596
DGRP-627
DGRP-630
DGRP-634
DGRP-639
DGRP-646
DGRP-703
DGRP-707
DGRP-712
DGRP-714
DGRP-716
DGRP-721
DGRP-730
DGRP-732
DGRP-737
DGRP-738
DGRP-748
DGRP-761
DGRP-765
DGRP-774
DGRP-776
DGRP-783
DGRP-786
DGRP-787
DGRP-790
DGRP-799
DGRP-801
DGRP-802
DGRP-808
DGRP-810
DGRP-812
DGRP-819
DGRP-820
DGRP-821
DGRP-822
DGRP-832
DGRP-837
DGRP-843
DGRP-850
DGRP-853
DGRP-855
DGRP-857
DGRP-859
DGRP-879
DGRP-882
DGRP-884
DGRP-887
DGRP-890
DGRP-892
DGRP-894
DGRP-897
DGRP-900
DGRP-907
DGRP-908
DGRP-911
DGRP-913
ore
ore2
ore3
w1
w1118
w2
)

#initialName=${initialNames[${LSB_JOBINDEX}]};
#R1=`ls -d -1 $trimDir/* | grep $initialName | grep _R1_`;
# matching R2, if any
#R2=`ls -d -1 $trimDir/* | grep $(echo $R1 | sed 's@.*/@@' | cut -f1,2 -d "_") | grep _R2_`;
# real GT
ID=${realGTs[${LSB_JOBINDEX}]};
echo $ID
# I use mapper as info
info=BWA
# I use genome as library
lib=dm6

# check, whatever or not sequencing was pair-end and align accordingly
#if [ "$R2" == "" ]; then 
#  mapWithBWA_SE $R1 $refGen $ID $info $lib $bwaDir
#else
#  mapWithBWA_PE $R1 $R2 $refGen $ID $info $lib $bwaDir
#fi

#java -jar /software/UHTS/Analysis/picard-tools/2.2.1/bin/picard.jar \
#          MarkDuplicates I=$bwaDir$ID"_"$info"_"$lib".align.srt.bam" \
#          O=$deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" \
#          CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
#          M=$deduplDir$ID"_"$info"_"$lib".metrix" \
#          REMOVE_DUPLICATES=TRUE
#samtools flagstat $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" > $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.flagstat"
#samtools idxstats $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" | cut -f 1,3 > $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.idxstat"
#samtools view -b $deduplDir$ID"_"$info"_"$lib".align.srt.dedupl.bam" chrM > $chrMDir$ID"_"$info"_"$lib".align.srt.dedupl.chrM.bam"
#samtools index $chrMDir$ID"_"$info"_"$lib".align.srt.dedupl.chrM.bam"

# DO NOT TRY TO RUN GENOTYPE GVCF ON DM6! IT PRODUCES VERY BIG FILES WHICH ARE
# IMPOSSIBLE TO PROCESS!
# BUT DO RUN GENOTYPE GVCF ON chrM OF DM6!

# Local realignment around indels
java -jar $GATKjar -T RealignerTargetCreator \
                   -R $refGen -I $chrMDir$ID"_"$info"_"$lib".align.srt.dedupl.chrM.bam" \
                   -o $targetDir$ID.intervals

# Replace read group
java -jar /software/UHTS/Analysis/picard-tools/2.2.1/bin/picard.jar \
          AddOrReplaceReadGroups I=$chrMDir$ID"_"$info"_"$lib".align.srt.dedupl.chrM.bam" \
          O=$replaceRGdir$ID.RGrepl.dedupl.bam SORT_ORDER=coordinate \
          RGID=$ID RGLB=$ID PL=illumina RGPU=illumina RGSM=Sample1 \
          CREATE_INDEX=True

# Perform realignment 
java -jar $GATKjar -T IndelRealigner -R $refGen \
          -I $replaceRGdir$ID.RGrepl.dedupl.bam \
          -targetIntervals $targetDir$ID.intervals \
          -o $realignDir$ID.realign.dedupl.bam

# HaplotypeCaller - create GVCF
java -jar $GATKjar -T HaplotypeCaller -nct 8 -R $refGen \
     -I $realignDir$ID.realign.dedupl.bam \
     -stand_call_conf 30 -stand_emit_conf 10 \
     --emitRefConfidence GVCF \
     -o $vcfDir$ID"_chrM.dm6.gatk.var.raw.g.vcf"

# Select variants only for chrM
java -jar $GATKjar -T SelectVariants -R $ref \
   -V $vcfDir$ID"_chrM.dm6.gatk.var.raw.g.vcf" \
   -o $vcfDir$ID"_chrM.only.dm6.gatk.var.raw.g.vcf" \
   -L chrM

mv $vcfDir$ID"_chrM.only.dm6.gatk.var.raw.g.vcf" $vcfDir$ID"_chrM.dm6.gatk.var.raw.g.vcf"

# add sample name, important for merging
toReplace=Sample1
sed -i "s@$toReplace@$ID@g" $vcfDir$ID"_chrM.dm6.gatk.var.raw.g.vcf"

exit 0;
