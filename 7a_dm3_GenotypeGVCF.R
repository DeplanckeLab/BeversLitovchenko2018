#!/usr/bin/env Rscript
# FILE DESCRIPTION: 7a_dm3_GenotypeGVCF -----------------------------------------
#
# DESCRIPTION : 
#
# USAGE: 
#
# OPTIONS:  none
# REQUIREMENTS:  data.table, lubridate, stringr
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:
# REVISION:

#setwd('~/Desktop/BitBucket/MitoSeq_RPJB_ML_Aug2017/')
setwd('~/Desktop/MitoSeq_RPJB_ML_Aug2017/')
source('17_functions.R')

# FUNCTIONS -------------------------------------------------------------------
#' deNovoVariants
#' Get do novo variants
#' @param testGR GRanges object with variants of the test sample
#' @param testGeno genotype table with variants of the test sample
#' @param refGR GRanges object with variants of the reference samples
#' @return 
deNovoVariants <- function(testGR, testGeno, refGR) {
  # detection of de novo variants - GRanges
  testGRdeNovo <- testGR[testGR %outside% refGR]
  if (length(testGRdeNovo) != 0) {
    # detection of de novo variants - genotype table
    testGenoDeNovo <- testGeno[names(testGRdeNovo), ]
    
    # heteroplasmy of de-novo variants
    testHeteroplDeNovo <- ifelse(testGenoDeNovo != '0/0' & 
                                 testGenoDeNovo != '1/1', T, F)
    
    # structural type of de novo variants
    testStrTypeDeNovo <- getVariantStructType(testGRdeNovo)
    
    # read support of de-novo variants
    testReadSuppDeNovo <- deplAD[names(testGRdeNovo)]
  } else {
    testGRdeNovo <- GRanges()
    testHeteroplDeNovo <- c()
    testStrTypeDeNovo <- c()
    testReadSuppDeNovo <- data.table(DGRP = character(), Variant = character(), 
                                     ReadCount = integer(), 
                                     ReadCountPerc = numeric())
  }
  result <- list(testGRdeNovo, testHeteroplDeNovo, testStrTypeDeNovo,
                 testReadSuppDeNovo)
  names(result) <- c('deNovoGRanges', 'deNovoHeteroplasmy',
                     'deNovoStructType', 'deNovoReadSupport')
  result
} 

#' compareToOneDGRP
#' Compare (genotype) one deplancke line to one dgrp line
#' @param oneDGRP one DGRP line (reference)
#' @param oneDepl one deplancke line
#' @return 
compareToOneDGRP <- function(oneDGRP, oneDepl) {
  oneDGRP <- gsub('\\.', NA, oneDGRP)
  both <- cbind(as.character(oneDepl), as.character(oneDGRP))
  rownames(both) <- names(oneDepl)
  both <- both[complete.cases(both), ]
  if (is.vector(both)) {
    both <- data.frame(t(both), stringsAsFactors = F)
    result <- list(sum(both[1] == both[2]),
                   data.frame(t(both[both[1] != both[2]])))
  } else {
    both <- data.frame(both, stringsAsFactors = F)
    result <- list(sum(both[, 1] == both[, 2]) / nrow(both),
                   both[both[, 1] != both[, 2], ])
  }
  names(result) <- c('PercentMatch', 'Notmatching')
  result
}

# INPUTS ----------------------------------------------------------------------
args <- commandArgs(trailingOnly = T)
# path to ALL VCFs 
allVCFdir <- 'results/genotyping/vcfForGenotyping'
#allVCFdir <- args[1]
#sampleName <- args[2]

# path to all deplancke vcfs
allDepl <- list.files(allVCFdir, "_GT_in_Depl.vcf.gz", T, T)
names(allDepl) <- sapply(allDepl, 
                         function(x) gsub("_GT_in_Depl.vcf.gz", "", 
                                          gsub(paste0(allVCFdir, '/'), '', x)))
#allDepl <- readRDS('allDepl.Rds')
# path to all DGRP vcfs
allDgrp <- list.files(allVCFdir, "_GT_in_DGRP2.vcf.gz", T, T)
names(allDgrp) <- sapply(allDgrp, 
                         function(x) gsub("_GT_in_DGRP2.vcf.gz", "", 
                                          gsub(paste0(allVCFdir, '/'), '', x)))
# allDgrp <- readRDS('allDgrp.Rds')
# autosomes names
autosomes <- c("chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet",
               "chr3R", "chr3RHet", "chr4", "chrX", "chrXHet", "chrYHet")
# structural type order
strTypeOrder <- c('SNP', 'DEL', 'INS', 'MNP', 'SUBS')

# OUTPUTS ---------------------------------------------------------------------
genotypingOvrView <- c()
deNovoOvrView <- c()
heteroplOvrView <- c()
allelesOvrView <- c()

# PERFORM GENOTYPING ----------------------------------------------------------
# it's modified slightly to run on vital-it
#for (sampleName in names(allDgrp)) {
  message(paste('Working with', sampleName, Sys.time()))
  
  # read-in vcf of dgrps genotyped on the same variants as deplancke line
  dgrpVcf <- readVcf(allDgrp[sampleName])
  dgrpGeno <- geno(dgrpVcf)$GT
  dgrpGR <- rowRanges(dgrpVcf)
  
  # READING IN
  # read-in one deplancke dgrp line genotyping
  deplVcf <- readVcf(allDepl[sampleName])
  deplGR <- rowRanges(deplVcf) # GRanges
  deplGR <- deplGR[deplGR$FILTER == 'PASS'] # filter out bad quality
  deplGR <- deplGR[seqnames(deplGR) %in% autosomes] # restrict to autosomes
  deplGeno <- geno(deplVcf)$GT # genotype table
  deplGeno <- as.data.frame(deplGeno[names(deplGR), ])
  colnames(deplGeno) <- sampleName
  # get read support of alternative allele
  deplAD <- as.data.table(getReadSupportOfAlt(deplGeno, geno(deplVcf)$AD))
  deplAD <- deplAD[deplAD$ReadCount != 0]
  setkey(deplAD, Variant) # read support
  deplGR <- deplGR[deplAD$Variant]
  deplGeno <- as.data.frame(deplGeno[names(deplGR), ])
  colnames(deplGeno) <- sampleName
  rownames(deplGeno) <- names(deplGR)
  message(paste('\t Finished reading', Sys.time()))
  
  # DE NOVO VARIANTS
  deNovoOneLine <- deNovoVariants(deplGR, deplGeno, dgrpGR)
  toAdd <- c(sampleName, length(deplGR), length(deNovoOneLine$deNovoGRanges), 
             sum(deNovoOneLine$deNovoHeteroplasmy), 
             table(deNovoOneLine$deNovoStructType)[strTypeOrder],
             summary(deNovoOneLine$deNovoReadSupport$ReadCount))
  deNovoOvrView <- rbind(deNovoOvrView, toAdd)
  message(paste('\t Finished de-novo variants', Sys.time()))
  
  # HETEROPLASMIC NOT DE NOVO
  deplGRhetero <- deplGR[deplGR %within% dgrpGR]
  deplGenoHetero <- deplGeno[names(deplGRhetero), ]
  names(deplGenoHetero) <- names(deplGRhetero)
  deplGenoHetero <- deplGenoHetero[deplGenoHetero != '0/0' &
                                   deplGenoHetero != '1/1']
  deplGRhetero <- deplGRhetero[names(deplGenoHetero)]
  if (length(deplGRhetero) != 0) {
    # structural type of heteroplasmic variants
    deplStrTypeHetero <- getVariantStructType(deplGRhetero)
    # read support of heteroplasmic variants
    deplReadSuppHetero <- deplAD[names(deplGRhetero)]
  } else {
    deplStrTypeHetero <- c()
    deplReadSuppHetero <- data.frame(DGRP = character(), Variant = character(),
                                     ReadCount = numeric(),
                                     ReadCountPerc = numeric())
  }
  
  toAdd <- c(sampleName, length(deplGR), length(deplGRhetero),
             table(deplStrTypeHetero)[strTypeOrder],
             summary(deplReadSuppHetero$ReadCount))
  heteroplOvrView <- rbind(heteroplOvrView, toAdd)
  message(paste('\t Finished heteroplasmy variants', Sys.time()))
  
  # RESTRICT DEPLANCKE TO THE VALID FOR GENOTYPING SITES
  deplGRovrl <- deplGR[deplGR %within% dgrpGR]
  deplGenoOvrl <- deplGeno[names(deplGRovrl), ]
  names(deplGenoOvrl) <- names(deplGRovrl)
  deplGenoOvrl <- deplGenoOvrl[deplGenoOvrl == '0/0' | deplGenoOvrl == '1/1']
  deplGRovrl <- deplGRovrl[names(deplGenoOvrl)]
  
  # RESTRICT REF DGRPs TO THE VALID FOR GENOTYPING SITES
  dgrpGRovrl <- dgrpGR[dgrpGR %within% deplGRovrl]
  dgrpGenoOvrl <- dgrpGeno[names(dgrpGRovrl), ]
  
  # CHECK ON ALLELES MATCH
  # check alt allele, sum is necesary in case there are MNPs, so if 1 allele 
  # overlaps, it's already enough
  alleles <- sapply(1:length(deplGRovrl), 
                    function(x) sum(deplGRovrl$ALT[[x]] == dgrpGRovrl$ALT[[x]]))
  if (length(deplGRovrl[alleles == 0]) != 0) {
    deplStrTypeAlleles <- getVariantStructType(deplGRovrl[alleles == 0]) 
    deplReadSuppAlleles <- deplAD[names(deplGRovrl[alleles == 0])]
  } else {
    deplStrTypeAlleles <- c()
    deplReadSuppAlleles <- data.frame(DGRP = character(), Variant = character(),
                                     ReadCount = numeric(),
                                     ReadCountPerc = numeric())
  }
  toAdd <- c(sampleName, length(deplGR), sum(alleles == 0),
             table(deplStrTypeAlleles)[strTypeOrder],
             summary(deplReadSuppAlleles$ReadCount))
  allelesOvrView <- rbind(allelesOvrView, toAdd)
  message(paste('\t Finished allele match', Sys.time()))
  
  # RESTRICT DEPLANCKE AND DGRPs TO THE CLEAN VARIANTS
  deplGRovrlAlleles <- deplGRovrl[alleles >= 1, ]
  deplGenoOvrlAlleles <- deplGenoOvrl[names(deplGRovrlAlleles)]
  deplGenoOvrlAlleles <- as.character(deplGenoOvrlAlleles)
  names(deplGenoOvrlAlleles) <- names(deplGRovrlAlleles)
  dgrpGRovrlAlleles <- dgrpGRovrl[alleles >= 1, ]
  dgrpGenoOvrlAlleles <- dgrpGenoOvrl[names(dgrpGRovrlAlleles), ]
  
  # GENOTYPE
  GTs <- apply(dgrpGenoOvrlAlleles, 2, compareToOneDGRP, deplGenoOvrlAlleles)
  GTsVals <- sapply(GTs, function(x) x$PercentMatch)
  GTsVals <- GTsVals[order(-GTsVals)]
  firstMatch <- names(GTsVals)[1]
  secondMatch <- names(GTsVals)[2]
  thirdMatch <- names(GTsVals)[3]
  # check, where do variants come from if they don't overlap with 1st match
  nonOvrlVars <- rownames(GTs[[firstMatch]]$Notmatching)
  deplGRnonOvrlVars <- deplGRovrlAlleles[nonOvrlVars]
  deplGenoNonOvrlVars <- deplGenoOvrlAlleles[names(deplGRnonOvrlVars)]
  dgrpGRnonOvrlVars <- dgrpGRovrlAlleles[dgrpGRovrlAlleles %within% deplGRnonOvrlVars]
  dgrpGenoNonOvrlVars <- dgrpGenoOvrlAlleles[names(dgrpGRnonOvrlVars), ]
  if (is.vector(dgrpGenoNonOvrlVars)) {
    GTsNotMatch <- lapply(dgrpGenoNonOvrlVars, compareToOneDGRP,
                          deplGenoNonOvrlVars)
  } else {
    GTsNotMatch <- apply(dgrpGenoNonOvrlVars, 2, compareToOneDGRP,
                         deplGenoNonOvrlVars)
  }
  GTsNotMatchVals <- sapply(GTsNotMatch, function(x) x$PercentMatch)
  GTsNotMatchVals <- GTsNotMatchVals[order(-GTsNotMatchVals)]
  firstMatchNotMatch <- names(GTsNotMatchVals)[1]
  
  toAdd <- c(sampleName, length(deplGR), length(deNovoOneLine$deNovoGRanges),
             length(deplGRhetero), sum(alleles == 0), 
             length(deplGR) - length(deNovoOneLine$deNovoGRanges) - 
                              length(deplGRhetero) - sum(alleles == 0),
             summary(deplAD[names(deplGRovrlAlleles)]$ReadCount), 
             firstMatch, secondMatch, thirdMatch,  GTsVals[1], GTsVals[2], 
             GTsVals[3], length(nonOvrlVars),
             GTsNotMatchVals[1], firstMatchNotMatch)
  genotypingOvrView <- rbind(genotypingOvrView, toAdd)
  message(paste('\t Finished fully', Sys.time()))
#}
colnames(deNovoOvrView)[1:9] <- c('Name', 'NumberOfGTsites', 'NumberOfDeNovo',
                                  'NumberOfHeteropl', strTypeOrder)
colnames(heteroplOvrView)[1:8] <- c('Name', 'NumberOfGTsites', 
                                    'NumberOfHeteropl', strTypeOrder)
colnames(genotypingOvrView)[1:6] <- c('Name', 'NumberOfGTsites', 
                                      'NumberOfDeNovo', 'NumberOfHeteropl',
                                      'NumberOfAlleleNotMatch', 'CleanSites')
colnames(genotypingOvrView)[13:21] <- c('firstMatch', 'secondMatch',
                                        'thirdMatch', 'firstMatchPerc', 
                                        'secondMatchPerc', 'thirdMatchPerc',
                                        'NumbNotOvrlVars', 
                                        'PercOfNotOvrlVarsOvrlWithOtherLine',
                                        'OtherLine')
saveRDS(deNovoOvrView, paste0('deNovoOvrView_', sampleName, '.Rds'))
saveRDS(heteroplOvrView, paste0('heteroplOvrView_', sampleName, '.Rds'))
saveRDS(genotypingOvrView, paste0('genotypingOvrView_', sampleName, '.Rds'))
message('Done!')

# BASH CODE TO RUN ON VITAL-IT ------------------------------------------------
#!/bin/bash

#BSUB -L /bin/bash
#BSUB -e GENO.%I.err
#BSUB -o GENO.%I.out
#BSUB -J GENO[1-241]
#BSUB -M 20000000
#BSUB -R rusage[mem=20000]
#BSUB -n 1
#BSUB -u maria.litovchenko@epfl.ch

export PATH=/software/bin:$PATH;
module use /software/module/;
module add R/3.3.2;

vcfDIR=/scratch/el/monthly/mlitovch/MitoSeq_RPJB_ML_Aug2017/vcfForGenotyping

samp=($(ls $vcfDIR | grep _Depl_DGRP2.vcf.gz | sed 's/_Depl_DGRP2.vcf.gz//g'))
zeroArr=( zero )
samples=("${zeroArr[@]}" "${samp[@]}")

sample=${samples[${LSB_JOBINDEX}]};

echo $sample
# run R script

Rscript --vanilla 30_genotyping_v2.R $vcfDIR $sample

exit 0;

# MERGE ALL RDS FROM VITAL-IT -------------------------------------------------
genotypingRds <- 'results/genotyping/genotypingRDS/'

deNovoOvrView <- lapply(list.files(genotypingRds, 'deNovoOvrView', T, T), 
                        readRDS)
deNovoOvrView <- as.data.table(do.call(rbind, deNovoOvrView))
setkey(deNovoOvrView, Name)
deNovoOvrView <- deNovoOvrView[names(allDepl)]
saveRDS(deNovoOvrView, 'deNovoOvrView_genotyping_v2.Rds')
write.table(deNovoOvrView, 'deNovoOvrView.csv', quote = F, sep = '\t', 
            col.names = T, row.names = F)

heteroplOvrView <- lapply(list.files(genotypingRds, 'heteroplOvrView',
                                     T, T), readRDS)
heteroplOvrView <- as.data.table(do.call(rbind, heteroplOvrView))
setkey(heteroplOvrView, Name)
heteroplOvrView <- heteroplOvrView[names(allDepl)]
saveRDS(heteroplOvrView, 'heteroplOvrView_genotyping_v2.Rds')
write.table(heteroplOvrView, 'heteroplOvrView.csv', quote = F, sep = '\t', 
            col.names = T, row.names = F)

genotypingOvrView <- lapply(list.files(genotypingRds, 'genotypingOvrView',
                                       T, T), readRDS)
genotypingOvrView <- as.data.table(do.call(rbind, genotypingOvrView))
setkey(genotypingOvrView, Name)
genotypingOvrView <- genotypingOvrView[names(allDepl)]
saveRDS(genotypingOvrView, 'ovrView_genotyping_v2.Rds')
write.table(genotypingOvrView, 'genotypingOvrView.csv', quote = F, sep = '\t', 
            col.names = T, row.names = F)
