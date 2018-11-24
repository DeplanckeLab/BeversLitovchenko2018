#BSUB -L /bin/bash
#BSUB -e NUMTS_genome.err
#BSUB -o NUMTS_genome.out
#BSUB -J NUMTS_genome
#BSUB -M 36000000
#BSUB -R rusage[mem=36000]
#BSUB -n 4
#BSUB -u maria.litovchenko@epfl.ch

#------------------------------------------------------------------------------
# STEP 1: CREATES NUMTS GENOME TO ALIGN THE DATA TO. BOTH MICE AND FLY.
#------------------------------------------------------------------------------

export PATH=/software/bin:$PATH;
module use /software/module/;
module add UHTS/Analysis/samtools/1.3;
module add UHTS/Aligner/bwa/0.7.13;
module add Blast/ncbi-blast/2.2.28+;
module add UHTS/Analysis/BEDTools/2.26.0;

# extract chrM for mouse
samtools faidx mm10.fa
samtools faidx mm10.fa chrMT > mm10_chrMT.fa
# extract chrM for fly
samtools faidx /scratch/el/monthly/mlitovch/RefGen/dm6/dm6.Wolb.fa chrM > dm6_chrM.fa

# blastn mouse genome to get the location of numts
# chrM, mm10 was aligned to the mouse nuclear genome (mm10) by BLASTn 
# (NCBI Blast-2.2.26) program using the default parameters. All hits with 
# length over 50 bp were kept as Numts fragments
blastn -query mm10_chrMT.fa -subject mm10.fa -outfmt 6 | cut -f2,4,9,10 | grep -v chrMT | awk -v OFS='\t' '$2 > 50  {print $1, $3, $4}' > mm10_numts.bed
blastn -query dm6_chrM.fa -subject /scratch/el/monthly/mlitovch/RefGen/dm6/dm6.Wolb.fa -outfmt 6 | cut -f2,4,9,10 | grep -v chrM | awk -v OFS='\t' '$2 > 50  {print $1, $3, $4}' > dm6_numts.bed

# create NUMTs genome
bedtools getfasta -fi mm10.fa -bed mm10_numts.bed -fo > mm10_numts.fa 
samtools faidx mm10_numts.fa
bwa index mm10_numts.fa
java -jar /scratch/el/monthly/mlitovch/tools/picard/picard.jar CreateSequenceDictionary \
    REFERENCE=mm10_numts.fa \ 
    OUTPUT=mm10_numts.dict 

bedtools getfasta -fi /scratch/el/monthly/mlitovch/RefGen/dm6/dm6.Wolb.fa -bed dm6_numts.bed -fo > dm6_numts.fa 
samtools faidx dm6_numts.fa
bwa index dm6_numts.fa
java -jar /scratch/el/monthly/mlitovch/tools/picard/picard.jar CreateSequenceDictionary \
    REFERENCE=dm6_numts.fa \
    OUTPUT=dm6_numts.dict

# index mito genomes
samtools faidx mm10_chrMT.fa
bwa index mm10_chrMT.fa
java -jar /scratch/el/monthly/mlitovch/tools/picard/picard.jar CreateSequenceDictionary \
    REFERENCE=mm10_chrMT.fa \
    OUTPUT=mm10_chrMT.dict
samtools faidx dm6_chrM.fa
bwa index dm6_chrM.fa
java -jar /scratch/el/monthly/mlitovch/tools/picard/picard.jar CreateSequenceDictionary \
    REFERENCE=dm6_chrM.fa \
    OUTPUT=dm6_chrM.dict

exit 0;
