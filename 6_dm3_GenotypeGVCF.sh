#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e 5_GVCFs.err
#BSUB -o 5_GVCFs.out
#BSUB -J 5_GVCFs
#BSUB -M 100000000
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add Development/java_jdk/1.8.0_66;
module add UHTS/Analysis/vcftools/0.1.12b;
module add R/3.3.2;

#------------------------------------------------------------------------------
# STEP 6: PERFORM MERGING OF RAW GVCF WITH GenotypeGVCFs
#------------------------------------------------------------------------------
outputdir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017;
GATKjar=/scratch/el/monthly/mlitovch/tools/GenomeAnalysisTK.jar
ref=/scratch/el/monthly/mlitovch/RefGen/dm3/dm3.Wolb.fa

rawGVCFdir=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/4_dm3_HC_d_rawGVCF

# list samples to merge and remove bad samples
toMerge=`ls -d -1 $rawGVCFdir/* | grep noIndelMNP.g.vcf$ | grep -v 177A12_L8 | grep -v 304A10_L5 | xargs | sed "s@/scratch@-V /scratch@g"`

# try first use mergeGVCF to make it faster

# do GenotypeGVCFs, it improves genotype calls
java -Xms64g -Xmx64g -jar $GATKjar -T GenotypeGVCFs \
     -allSites -stand_call_conf 30 -stand_emit_conf 10 \
     -R $ref -nt 4 --variant $toMerge \
     -o $outputdir/MitoSeq_AllRuns_dm3.gvcf.vcf
echo "Finished GenotypeGVCFs"

# select only SNPs
java -Xms64g -Xmx64g -jar $GATKjar -T SelectVariants \
     -R $ref -V $outputdir/MitoSeq_AllRuns_dm3.gvcf.vcf \
     -selectType SNP --selectTypeToExclude MNP -select "DP > 5" \
     -o $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.vcf
echo "Finished SelectVariants"

# filter out bad quality
java -Xms64g -Xmx64g -jar $GATKjar -T VariantFiltration \
     -R $ref -V $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.vcf \
     -filterName FS -filter "FS > 30.0" -filterName QD \
     -filter "QD < 2.0" \
     -o $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr.vcf
echo "Finished VariantFiltration"

# extract only genotype info
vcftools --vcf $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr.vcf \
         --extract-FORMAT-info GT \
         --out $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr
echo "Finished GT extraction"

# create table suitable for R
sed -i 's@\.\/\.@-@g' $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr.GT.FORMAT
sed -i 's@0\/0@0@g' $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr.GT.FORMAT
sed -i 's@1\/1@2@g' $outputdir/MitoSeq_AllRuns_dm3.gvcf.snps.fltr.GT.FORMAT
echo "Finished preparation of file for R"

exit 0;
