#!/usr/bin/env Rscript
# FILE: 14a_GRDs -------------------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: Performs chi-tests between nuclear and mitochondrial variants
#              to detect genome ratio disturbance (GRD). Suitable for vital it.
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: Package variant annotation is depricated on vital it and can not read
#       our mitovcf, so just same mitoVcfGeno as Rds and pass it to the script
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  26.01.2018
# REVISION: 26.01.2018

baseDir <- '~/Desktop/MitoSeq_RPJB_ML_Aug2017/'
#baseDir <- '/scratch/el/monthly/mlitovch/GRD_DGRP/'
setwd(baseDir)
source('17_functions.R')

# PATHS TO INPUTS -------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
beversDm6Path <- args[1]
mitoVarIndex <- as.integer(args[2])
dgrp2Dm6Path <- args[3]
saveTo <- args[4]

# bevers mitochronrial variants
beversDm6Path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
# for Pool
# beversDm6Path <- 'results/Pool_Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
# where coding part of mtDNA ends
codingEnds <- 14917
# cut-off on MAF
MAFcutoff <- 0.05
# This number ensures that we took all the mt variants from our dataset, but
# nothing extra
#MACcutoff <- 8.5 #  for Pool. 
# how many genotyped samples are required?
sampsWithGenoCut <- 150
# minimum frequncy of the alele
minAlleleFreqNuclMito <- 0.05

# Read in BEVERS data ---------------------------------------------------------
mitoVcfdm6 <- readInMitoVCF(beversDm6Path, 'dm6', codingEnds, 
                            removeRefLines = T, interGen = "merge")
mitoVcfGeno <- mitoVcfdm6[[2]]
# remove MNPs
mitoVcfGeno <- mitoVcfGeno[getVariantStructType(mitoVcfdm6[[1]]) != 'MNP', ]

# there are only 33 variants which pass MAF > 5% cut off. However, many of them
# have the same "configuration" in the sequenced DGRP lines aka clusters. I 
# will scan only 1 variant per cluster, which results into 12 variants
mitoClust <- readRDS('results/mitoClust.Rds')
mitoVcfGeno <- mitoVcfGeno[mitoClust[tested == 1]$MitoVar, ]

# I run GRD calculation with heteroplasmy in mito variants and that gave me
# inflated numbers of GRDs. I exclude heteroplasmic variants in individuals
# I do not though the whole variant away!
# Although the variants should not be heteroplasmic...
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/0", ".", x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/1", ".", x))
# select variants with number of genotyped samples > cut off
mitoVcfGeno <- mitoVcfGeno[apply(mitoVcfGeno, 1, 
                                 function(x) sum(x != './.' & 
                                                 x != '.') >= sampsWithGenoCut), ]
# select variants with MAF > MAFcutoff. All variants should pass.
mitoVcfGeno <- mitoVcfGeno[apply(mitoVcfGeno, 1,
                                 function(x) calcMacMaf(x)['MAF'] >= MAFcutoff),]
# for Pool
#mitoVcfGeno <- mitoVcfGeno[apply(mitoVcfGeno, 1,
#                              function(x) calcMacMaf(x)['MAC'] >= MACcutoff),]

mitoVcfdm6 <- NULL # clean up memory
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/0", 1, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/1", 3, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/0", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/1", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("\\.", 0, x))
mitoVcfGeno <- as.matrix(mitoVcfGeno)

# Read in DGRP2 DM6 vcf GENO --------------------------------------------------
dgrp2dm6Geno <- fread(dgrp2Dm6Path)
addZero <- grepl('DGRP-..$', colnames(dgrp2dm6Geno))
setnames(dgrp2dm6Geno, colnames(dgrp2dm6Geno)[addZero],
         gsub('DGRP-', 'DGRP-0', colnames(dgrp2dm6Geno)[addZero]))
dgrp2dm6Geno <- cbind(ID = paste0(dgrp2dm6Geno$CHROM, "_", dgrp2dm6Geno$POS),
                      dgrp2dm6Geno)
dgrp2dm6Geno[, CHROM := NULL]
dgrp2dm6Geno[, POS := NULL]
# some positions will be repetetive in the dgrp2dm6Geno - they are MNPs
# I remove them
nuclMNPs <- dgrp2dm6Geno$ID[duplicated(dgrp2dm6Geno$ID)]
dgrp2dm6Geno <- dgrp2dm6Geno[!dgrp2dm6Geno$ID %in% nuclMNPs, ]

# translate to matrix
dgrp2dm6Geno <- as.matrix(dgrp2dm6Geno)
rownames(dgrp2dm6Geno) <- dgrp2dm6Geno[, 1]
dgrp2dm6Geno <- dgrp2dm6Geno[, -1]
mitoVcfGeno <- mitoVcfGeno[, colnames(dgrp2dm6Geno)]

# Perform calculations --------------------------------------------------------
print(paste('Started at', Sys.time()))
start.time <- Sys.time()
chiStats <- apply(dgrp2dm6Geno, 1, autoChiSq, mitoVcfGeno[mitoVarIndex, ],
                  sampsWithGenoCut, minAlleleFreqNuclMito)
chiStats <- as.data.frame(t(chiStats))
colnames(chiStats) <- c('X2', 'pvalue', 'complete')
chiStats <- cbind(NuclVar = rownames(chiStats),       
                  MitoVar = rownames(mitoVcfGeno)[mitoVarIndex],    
                  chiStats)
chiStats <- as.data.table(chiStats)
saveRDS(chiStats, saveTo)
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste0('Time taken: ', time.taken))
