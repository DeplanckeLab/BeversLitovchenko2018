#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 11_annot.err
#BSUB -o 11_annot.out
#BSUB -J 11_annot
#BSUB -M 64000000
#BSUB -R rusage[mem=64000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 10: MERGE MITOCHANDRIAL GVCFs INTO ONE
#------------------------------------------------------------------------------
export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Analysis/vcftools/0.1.14;

GATKjar=/scratch/el/monthly/mlitovch/tools/GenomeAnalysisTK.jar
ref=/scratch/el/monthly/mlitovch/RefGen/dm6/dm6.Wolb.fa
outputdir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017

toMerge=`ls -d -1 9_dm6_HC_d_rawVCF/* | grep _chrM.dm6.gatk.var.raw.g.vcf$ | xargs | sed "s@9_dm6_HC@-V 9_dm6_HC@g"`

# do GenotypeGVCFs, it improves genotype calls
java -Xms64g -Xmx64g -jar $GATKjar -T GenotypeGVCFs \
     -stand_call_conf 30 -stand_emit_conf 10 \
     -R $ref -nt 4 $toMerge \
     -o $outputdir/Bevers_dm6_chrM_GenotypeGVCF.vcf
echo "Finished GenotypeGVCFs"

# filter out bad quality
# SEE https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set
# FOR DETAILS

# This filter will mark all variants with low quality and will set all
# genotypes of variants which do not pass to no call

# DP depth of coverage

# FS = Phred-scaled p-value using Fisher’s Exact Test to detect strand bias 
# (the variation being seen on only the forward or only the reverse 
# strand) in the reads. More bias is indicative of false positive calls.

# QD = QualByDepth (QD) 2.0 This is the variant confidence (from the QUAL 
# field) divided by the unfiltered depth of non-reference samples.
# THIS FILTER ISN'T POSSIBLE, BECAUSE VARIANT SCORE RECALIBRATION ISN'T 
# POSSIBLE TO DO DUE TO THE ABSENCE OF KNOWN FILE FOR DROSOPHILA

java -Xms64g -Xmx64g -jar $GATKjar -T VariantFiltration \
     --genotypeFilterExpression "DP < 10 || FS > 60.0 || MQ < 10" \
     --genotypeFilterName "Low QUAL" \
     --setFilteredGtToNocall \
     -R $ref \
     -V $outputdir/Bevers_dm6_chrM_GenotypeGVCF.vcf \
     -o $outputdir/Bevers_dm6_chrM_GenotypeGVCF.flrt.vcf
echo "Finished VariantFiltration"

java -jar /scratch/el/monthly/mlitovch/tools/snpEff/snpEff.jar \
     -no-downstream -no-upstream -no-utr dm6.mitoMANUALCHECK \
     $outputdir/Bevers_dm6_chrM_GenotypeGVCF.flrt.vcf > $outputdir/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf
echo "Finished annotation"


exit 0;
