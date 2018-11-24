#!/bin/bash

# Based on Genetic Incompatibilities are Widespread Within Species
# Putative variant be supported by a minimum of five reads
# No MNPs
# Excluded sites wherein fewer than 150 individuals have a supported genotype
# Where the minor allele was present in fewer than 10 individuals
# Where more than 15% percent of individuals with data had heterozygous genotypes

#------------------------------------------------------------------------------
# STEP 13: MOVE DGRP2 VCF FROM DM3 TO DM6, SELECT VARIANTS BY MAF
#------------------------------------------------------------------------------

# Download DGRP2
wget dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf

# Move it to dm6 
# stupid crossmap doesn't work if I modify chain file by removing chr
# so I have to modify everything else
sed -i 's/^/chr/g' dgrp2.vcf
sed -i 's/chr#/#/g' dgrp2.vcf
wget http://hgdownload.cse.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
refGen=/home/litovche/Documents/RefGen/dm6/dm6.Wolb.fa
# do cross map
CrossMap.py vcf dm3ToDm6.over.chain.gz dgrp2.vcf $refGen dgrp2_dm6.vcf
sed -i "s/line_/DGRP-/g" dgrp2_dm6.vcf
bgzip -c dgrp2_dm6.vcf > dgrp2_dm6.vcf.gz
tabix -p vcf dgrp2_dm6.vcf.gz
rm dgrp2.vcf
rm dgrp2_dm6.vcf

# Restrict it to only our samples
vcftools --gzvcf /home/litovche/Documents/RefGen/dm6/dgrp2_dm6.vcf.gz \
         --indv DGRP-21 --indv DGRP-26 --indv DGRP-28 --indv DGRP-31 \
         --indv DGRP-32 --indv DGRP-38 --indv DGRP-41 --indv DGRP-45 \
         --indv DGRP-48 --indv DGRP-57 --indv DGRP-59 --indv DGRP-69 \
	 --indv DGRP-73 --indv DGRP-75 --indv DGRP-83 --indv DGRP-85 \
	 --indv DGRP-88 --indv DGRP-91 --indv DGRP-93 --indv DGRP-100 \
	 --indv DGRP-101 --indv DGRP-105 --indv DGRP-109 --indv DGRP-129 \
	 --indv DGRP-136 --indv DGRP-138 --indv DGRP-142 --indv DGRP-149 \
	 --indv DGRP-153 --indv DGRP-158 --indv DGRP-161 --indv DGRP-181 \
	 --indv DGRP-189 --indv DGRP-195 --indv DGRP-208 --indv DGRP-217 \
	 --indv DGRP-223 --indv DGRP-228 --indv DGRP-229 --indv DGRP-235 \
	 --indv DGRP-237 --indv DGRP-239 --indv DGRP-280 --indv DGRP-287 \
	 --indv DGRP-301 --indv DGRP-303 --indv DGRP-306 --indv DGRP-307 \
	 --indv DGRP-309 --indv DGRP-313 --indv DGRP-315 --indv DGRP-317 \
	 --indv DGRP-318 --indv DGRP-319 --indv DGRP-321 --indv DGRP-324 \
     --indv DGRP-332 --indv DGRP-336 --indv DGRP-340 \ #--indv DGRP-338
	 --indv DGRP-348 --indv DGRP-352 --indv DGRP-354 --indv DGRP-355 \
     --indv DGRP-357 --indv DGRP-358 --indv DGRP-359 \ # --indv DGRP-356
	 --indv DGRP-361 --indv DGRP-362 --indv DGRP-365 --indv DGRP-367 \
	 --indv DGRP-370 --indv DGRP-371 --indv DGRP-373 --indv DGRP-375 \
	 --indv DGRP-377 --indv DGRP-379 --indv DGRP-380 --indv DGRP-381 \
	 --indv DGRP-382 --indv DGRP-383 --indv DGRP-385 --indv DGRP-386 \
	 --indv DGRP-390 --indv DGRP-392 --indv DGRP-395 --indv DGRP-397 \
	 --indv DGRP-399 --indv DGRP-405 --indv DGRP-406 --indv DGRP-409 \
	 --indv DGRP-426 --indv DGRP-427 --indv DGRP-437 --indv DGRP-440 \
	 --indv DGRP-441 --indv DGRP-443 --indv DGRP-486 --indv DGRP-491 \
	 --indv DGRP-492 --indv DGRP-502 --indv DGRP-508 --indv DGRP-509 \
     --indv DGRP-513 --indv DGRP-517 --indv DGRP-528 --indv DGRP-530 \
	 --indv DGRP-535 --indv DGRP-551 --indv DGRP-555 --indv DGRP-559 \
	 --indv DGRP-563 --indv DGRP-566 --indv DGRP-584 --indv DGRP-596 \
	 --indv DGRP-627 --indv DGRP-630 --indv DGRP-634 --indv DGRP-639 \
	 --indv DGRP-646 --indv DGRP-703 --indv DGRP-707 --indv DGRP-712 \
	 --indv DGRP-714 --indv DGRP-716 --indv DGRP-721 --indv DGRP-730 \
	 --indv DGRP-732 --indv DGRP-737 --indv DGRP-738 --indv DGRP-748 \
	 --indv DGRP-761 --indv DGRP-765 --indv DGRP-774 --indv DGRP-776 \
	 --indv DGRP-783 --indv DGRP-786 --indv DGRP-787 --indv DGRP-790 \
	 --indv DGRP-799 --indv DGRP-801 --indv DGRP-802 --indv DGRP-808 \
	 --indv DGRP-810 --indv DGRP-812 --indv DGRP-819 --indv DGRP-820 \
	 --indv DGRP-821 --indv DGRP-822 --indv DGRP-832 --indv DGRP-837 \
	 --indv DGRP-843 --indv DGRP-850 --indv DGRP-853 --indv DGRP-855 \
	 --indv DGRP-857 --indv DGRP-859 --indv DGRP-879 --indv DGRP-882 \
	 --indv DGRP-884 --indv DGRP-887 --indv DGRP-890 --indv DGRP-892 \
	 --indv DGRP-894 --indv DGRP-897 --indv DGRP-900 --indv DGRP-907 \
	 --indv DGRP-908 --indv DGRP-911 --indv DGRP-913 \
         --recode --recode-INFO-all \
         --out dgrp2_dm6_BeversLinesOnly

# Remove variants with less than 150 lines genotyped and >5% MAF only
vcftools --vcf dgrp2_dm6_BeversLinesOnly.recode.vcf \
	     --max-missing 0.88 \
         --maf 0.05 \
         --min-alleles 2 --max-alleles 2 \
         --recode --recode-INFO-all \
         --out dgrp2_dm6_BeversLinesOnly_MAF005
bgzip -c dgrp2_dm6_BeversLinesOnly_MAF005.recode.vcf > dgrp2_dm6_BeversLinesOnly_MAF005.vcf.gz
tabix -p vcf dgrp2_dm6_BeversLinesOnly_MAF005.vcf.gz

# divide in chromosomes

chroms=( chr2L chr3L chr2R chr3R chrX chr4 )

for chr in "${chroms[@]}"
do
    vcftools --gzvcf dgrp2_dm6_BeversLinesOnly_MAF005.vcf.gz \
            --chr $chr \
            --recode --recode-INFO-all \
            --out "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr
    bgzip -c "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".recode.vcf" > "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".recode.vcf.gz"
    tabix -p vcf "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".recode.vcf.gz"
    rm "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".recode.vcf"
    vcftools --gzvcf "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".recode.vcf.gz" \
            --extract-FORMAT-info GT \
            --out "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr
    sed -i "s@0/0@1@g" "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"
    sed -i "s@1/1@3@g" "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"
    sed -i "s@0/1@2@g" "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"
    sed -i "s@1/0@2@g" "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"
    sed -i "s@./.@0@g" "dgrp2_dm6_BeversLinesOnly_MAF005_"$chr".GT.FORMAT"
done
