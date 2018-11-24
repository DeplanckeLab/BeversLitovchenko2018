#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 14_mpileup_rich.err
#BSUB -o 14_mpileup_rich.out
#BSUB -J 14_mpileup_rich
#BSUB -M 64000000
#BSUB -R rusage[mem=64000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 3: VARIANT CALLING IN MITOCHONDRIAL GENOMES ACCORDING TO THE
# RICHARDON PAPER
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Analysis/vcftools/0.1.14;

GATKjar=/scratch/el/monthly/mlitovch/tools/GenomeAnalysisTK.jar
ref=/scratch/el/monthly/mlitovch/RefGen/richardson/richardson.fa
deduplDir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/3_dedupl_richardson/

# list samples to merge and remove bad samples
toMerge=`ls -d -1 $deduplDir* | grep align.srt.dedupl.bam$ | grep -v 177A12_L8 | grep -v 304A10_L5 | xargs`

#samtools mpileup -ugf $ref $toMerge | bcftools call -vmO z -o MitoSeq_AllRuns_richardson.vcf.gz  
#tabix -p vcf MitoSeq_AllRuns_richardson.vcf.gz
#bcftools filter -O z -o MitoSeq_AllRuns_richardson.flt.vcf.gz -s LOWQUAL -i'%QUAL>10' MitoSeq_AllRuns_richardson.vcf.gz

# extract chrU only
#vcftools --gzvcf MitoSeq_AllRuns_richardson.flt.vcf.gz \
#         --chr chrU --recode \
#         --out MitoSeq_AllRuns_richardson_chrU.flt
mv MitoSeq_AllRuns_richardson_chrU.flt.recode.vcf MitoSeq_AllRuns_richardson_chrU.flt.vcf

# select samples
java -jar $GATKjar -T SelectVariants -R $ref \
     -V MitoSeq_AllRuns_richardson_chrU.flt.vcf \
     -o MitoSeq_AllRuns_richardson_chrU_rigthGT.fltr.vcf -sn Berk1_S71 -sn Berk2_S98 -sn 21E9_L5 -sn 26D9_L5 -sn 28A4_L8 -sn 31D7_L5 -sn 32F9_L5 -sn 508C9_L5 -sn 41_S94 -sn 45_S57 -sn 48B12_L8 -sn 57B_RB -sn 101_N712 -sn 69_S59 -sn 73E8_L5 -sn 75A5_L8 -sn 83B5_L8 -sn 85B9_L5 -sn 88_S72 -sn 91A10_L8 -sn 93A6_L8 -sn 100_S96 -sn 59_S91 -sn 105_S89 -sn 859_N712 -sn 129A2_L8 -sn 136B8_L8 -sn 138_S109 -sn 142_S56 -sn 149A8_L8 -sn 153A12_L5 -sn 362A3_L5 -sn 161B6_L8 -sn 181G9_L5 -sn 189E10_L8 -sn 195C6_L8 -sn 208_S121 -sn 217B1_L8 -sn 223B11_L8 -sn 228A9_L8 -sn 229C10_L8 -sn 235F1_L5 -sn 237D10_L8 -sn 239G6_L5 -sn 280_S120 -sn 287C3_L8 -sn 301E4_L5 -sn 303C4_L8 -sn 306A3_L8 -sn 307D1_L8 -sn 309F2_L8 -sn 313A7_L8 -sn 315A1_L8 -sn 317C1_L5 -sn 318D8_L8 -sn 319B5_L5 -sn 321B7_L8 -sn 324A5_L5 -sn 332_S108 -sn 336_S50 -sn 338E11_L8 -sn 737C5_L8 -sn 348F10_L8 -sn 352G7_L5 -sn 354_Ad2 -sn 355C9_L8 -sn 356F8_L5 -sn 357G12_L8 -sn 358_S44 -sn 359B3_L8 -sn 361E2_L8 -sn 362E4_L8 -sn 365F2_L5 -sn 367_Ad2 -sn 819C10_L5 -sn 371C3_L5 -sn 373G2_L8 -sn 375G8_L8 -sn 377G5_L8 -sn 379A7_L5 -sn 383G10_L8 -sn 227_N712 -sn 382A9_L5 -sn 383_S48 -sn 385G2_L5 -sn 386E10_L5 -sn 390F12_L8 -sn 392B4_L5 -sn 395C4_L5 -sn 397F3_L8 -sn 399F11_L8 -sn 405_S46 -sn 406E3_L8 -sn 409A1_L5 -sn 426G1_L8 -sn 427G3_L8 -sn 437G6_L8 -sn 440E7_L5 -sn 441B2_L5 -sn 443E5_L8 -sn 486B10_L5 -sn 491_S45 -sn 492H2_L8 -sn 502E5_L5 -sn 508_S81 -sn 509D10_L5 -sn 513B2_L8 -sn 517_S122 -sn 42C1_L8 -sn 530E1_L5 -sn 535F3_L5 -sn 551E6_L5 -sn 555E12_L8 -sn 559F4_L8 -sn 563_S74 -sn 566E6_L8 -sn 584E7_L8 -sn 596_S110 -sn 627E9_L8 -sn 630F1_L8 -sn 634D12_L8 -sn 639G11_L5 -sn 235_N712 -sn 703G4_L5 -sn 707B10_L8 -sn 712B4_L8 -sn 714_S117 -sn 716C12_L8 -sn 721B1_L5 -sn 730_S70 -sn 732C11_L8 -sn 340D9_L8 -sn 738D5_L8 -sn 360_S115 -sn 761_S83 -sn 765G10_L5 -sn 774E12_L5 -sn 776C2_L5 -sn 783E8_L8 -sn 786G5_L5 -sn 787A6_L5 -sn 790D12_L5 -sn 799_S100 -sn 801C7_L8 -sn 802A11_L5 -sn 808F5_L5 -sn 810D3_L5 -sn 812_S76 -sn 370B8_L5 -sn 820_S68 -sn 821B11_L5 -sn 822F4_L5 -sn 832_S65 -sn 837G7_L8 -sn 843_S53 -sn 849_S93 -sn 853H6_L8 -sn 855F8_L8 -sn 857_S77 -sn 109_S103 -sn 879F7_L5 -sn 882G8_L5 -sn 884_S52 -sn 887_S101 -sn 890H7_L8 -sn 892_S107 -sn 894_S80 -sn 897C8_L5 -sn 900G1_L5 -sn 907H4_L8 -sn 908H5_L8 -sn 911H11_L8 -sn 913H12_L8 -sn ore_S64 -sn ore2F12_L5 -sn ore3D8_L5 -sn w1E2_L5 -sn w1118_N712 -sn w2A2_L5

# replace manually names of lines! check the order first
# like below is checked
# DGRP100	DGRP059	DGRP105 DGRP859 DGRP129 DGRP136 DGRP138 DGRP142 DGRP149 DGRP153 DGRP161 DGRP181 DGRP189 DGRP195 DGRP208 DGRP217 DGRP021 DGRP223 DGRP381 DGRP228 DGRP229 DGRP235 DGRP646 DGRP237 DGRP239 DGRP026 DGRP280 DGRP287 DGRP028 DGRP301 DGRP303 DGRP306 DGRP307 DGRP309 DGRP313 DGRP315 DGRP317 DGRP318 DGRP319 DGRP031 DGRP321 DGRP324 DGRP032 DGRP332 DGRP336 DGRP338 DGRP737 DGRP348 DGRP352 DGRP354 DGRP355 DGRP356 DGRP357 DGRP358 DGRP359 DGRP748 DGRP361 DGRP158 DGRP362 DGRP365 DGRP367 DGRP819 DGRP371 DGRP373 DGRP375 DGRP377 DGRP379 DGRP382 DGRP383 DGRP380 DGRP385 DGRP386 DGRP390 DGRP392 DGRP395 DGRP397 DGRP399 DGRP405 DGRP406 DGRP409 DGRP041 DGRP426 DGRP427 DGRP528 DGRP437 DGRP440 DGRP441 DGRP443 DGRP045 DGRP486 DGRP048 DGRP491 DGRP492 DGRP502 DGRP038 DGRP508 DGRP509 DGRP513 DGRP517 DGRP530 DGRP535 DGRP551 DGRP555 DGRP559 DGRP563 DGRP566 DGRP057 DGRP584 DGRP596 DGRP101 DGRP627 DGRP630 DGRP634 DGRP639 DGRP069 DGRP703 DGRP707 DGRP712 DGRP714 DGRP716 DGRP721 DGRP730 DGRP732 DGRP340 DGRP738 DGRP073 DGRP075 DGRP761 DGRP765 DGRP774 DGRP776 DGRP783 DGRP786 DGRP787 DGRP790 DGRP799 DGRP801 DGRP802 DGRP808 DGRP810 DGRP812 DGRP370 DGRP820 DGRP821 DGRP822 DGRP832 DGRP837 DGRP083 DGRP843 DGRP850 DGRP853 DGRP855 DGRP857 DGRP109 DGRP085 DGRP879 DGRP882 DGRP884 DGRP887 DGRP088 DGRP890 DGRP892 DGRP894 DGRP897 DGRP900 DGRP907 DGRP908 DGRP911 DGRP913 DGRP091 DGRP093 Berk1   Berk2   ore2    ore3    ore     w1118   w1      w2
exit 0;
