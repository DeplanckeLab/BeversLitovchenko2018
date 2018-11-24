# FILE: 2_functions -----------------------------------------------------------
#
# USAGE: 
#
#  DESCRIPTION: contains functions for reading-in VCF file, annotate and charac
#  terize variants
#
#  OPTIONS:  none
#  REQUIREMENTS:  none
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  1
#  CREATED:  26.07.2017
#  REVISION: 18.01.2017

library(biomaRt)
library(corrplot)
library(ggplot2)
library(ggsignif)
library(data.table)
library(effsize)
library(factoextra)
library(klaR)
library(NbClust)
library(seqinr)
library(sqldf)
library(plotrix)
library(topGO)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(VariantAnnotation)
library(venn)
library(Vennerable)

txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene

# Reading and processing VCF --------------------------------------------------
#' genoToRefHetAlt
#' Converts genotype table from 0/0, 1/1, 2/2, etc to REF, ALT, HET
#' @param VOIgeno genotype table (as geno(vcf)$GT) for the variants of interest
#' @return genotype table with REF, ALT, HET
genoToRefHetAlt <- function(VOIgeno) {
  VOIgeno <- apply(VOIgeno, 1, function(x) gsub('/', '', x))
  
  # check presence of more than 9 alleles. If value in geno table longer than 2
  # characters, than there was definetelly > 9 alleles
  numbAlleles <- apply(VOIgeno, 1, function(x) sum(nchar(x) > 2))
  if (length(which(numbAlleles > 0)) > 0) {
    stop(paste('ERROR: number of alleles for variant(s)', rownames(VOIgeno),
               'is more than 9!'))
  }
  
  VOIgeno <- gsub('00', 'REF', VOIgeno)
  VOIgeno <- gsub('0[1-9]|[1-9]0', 'HET', VOIgeno)
  VOIgeno <- gsub('[1-9]{2}', 'ALT', VOIgeno)
  VOIgeno
}

#' getAllelicDepth
#' Returns data frame with Variant, Sample, TotalReads for variants, 
#' PercAlt for the percentage of alternative count
#' @param variantCode, i.e. chrM:10088_C/T
#' @param vcfGeno genotype table
#' @param vcfAD allele depth table
#' @return data frame
getAllelicDepth <- function(variantCode, vcfGeno, vcfAD) {
  vcfAdVOI <- vcfAD[variantCode, ]
  if (max(unlist(lapply(vcfAdVOI, length))) == 2) {
    samplesWithVar <- colnames(vcfGeno)[vcfGeno[variantCode, ] != '0/0' &
                                        vcfGeno[variantCode, ] != '.']
    percOfAlt <- sapply(samplesWithVar, 
                      function(x) 100 * vcfAdVOI[x][[1]][2] / 
                                  (vcfAdVOI[x][[1]][1] + vcfAdVOI[x][[1]][2]))
    totalReads <- sapply(samplesWithVar, 
                         function(x) vcfAdVOI[x][[1]][1] + 
                                     vcfAdVOI[x][[1]][2])
    result <- data.frame(Variant = variantCode,
                         Sample = samplesWithVar, 
                         TotalReads = totalReads,
                         PercAlt = percOfAlt)
  } else {
    stop('ERROR: MNP are not supported')
  }
  result
}

#' mergeInterGenVars
#' Merges geno information for the intergenic variants
#' @param VOIgeno genotype table (as geno(vcf)$GT) for the variants of interest
#' @param newVarCode string-code for the new variant
#' @return vector with REF if all were REF and ALT otherwise
mergeInterGenVars <- function(VOIgeno, newVarCode = NA) {
  result <- rep('ALT', ncol(t(genoToRefHetAlt(VOIgeno))))
  result[apply(VOIgeno, 2, 
               function(x) all(x == 'REF'))] <- 'REF'
  result <- as.data.frame(t(result))
  colnames(result) <- colnames(VOIgeno)
  if (!is.na(newVarCode)) {
    rownames(result) <- newVarCode
  }
  result
}

#' readInAlleleDepth
#' Reads in allelic depth info from vcf
#' @param vcfpath path to vcf
#' @param genomeVer version of the genome
#' @return data frame of lists
readInAlleleDepth <- function(vcfpath, genomeVer) {
  vcf <- readVcf(vcfpath, genomeVer)
  vcfAD <- geno(vcf)$AD
  vcfAD
}

#' readInMitoVCF
#' Function to specificly read-in VCF file for Bevers study
#' @param mitoVcfpath path to VCF file
#' @param genomeVer genome version, i.e. dm6
#' @param codingEnds where coding region ends
#' @param removeRefLines whatever reference lines (w- and ore) be removed
#' @param interGen whatever region 5950 - 5975 should be merged (merge), 
#'        removed (remove) or untouched (NA)
#' @return list, containing GRanges table and Genotype table
readInMitoVCF <- function(vcfpath, genomeVer, codingEnds = 14917, 
                          removeRefLines = T, interGen = "merge") {
  vcf <- readVcf(vcfpath, genomeVer)
  vcfGeno <- geno(vcf)$GT
  
  # restrict to the coding part
  vcfGR <- rowRanges(vcf)
  if(!is.na(codingEnds)) {
    vcfGR <- vcfGR[start(vcfGR) < codingEnds]
  }
  vcfGeno <- vcfGeno[names(vcfGR), ]
  
  # remove bad lines and reference lines if requested
  if (removeRefLines) {
    badLines <-"ore|w1|w2|Berk|DGRP-338|DGRP-356|DGRP338|DGRP356"
    vcfGeno <- vcfGeno[, !grepl(badLines, colnames(vcfGeno))]
  } else { # w1118 is bad, always remove it! Berk2 also has a bad quality
    badLines <- "w1118|Berk2|ore3|DGRP-338|DGRP-356|DGRP338|DGRP356"
    vcfGeno <- vcfGeno[, !grepl(badLines, colnames(vcfGeno))]
    message('Removed DGRP-338, DGRP-356, w1118, ore3, Berk2')
  }
  
  # remove variants where all of the DGRPs are not genotyped
  nVarsBefore <- nrow(vcfGeno)
  varsBefore <- rownames(vcfGeno)
  vcfGeno <- vcfGeno[!apply(vcfGeno, 1, 
                            function(x) all(x == '.')), ]
  message(paste('Removed', nVarsBefore - nrow(vcfGeno), 
                'out of', nVarsBefore, 'variants because of all ./. :\n',
                paste(setdiff(varsBefore, rownames(vcfGeno)),
                      collapse = '\n ')))
  
  # and there's no variance because those vars were specific to ore/w-
  nVarsBefore <- nrow(vcfGeno)
  varsBefore <- rownames(vcfGeno)
  vcfGeno <- vcfGeno[!apply(vcfGeno, 1, 
                            function(x) sum(x == '0/0') + 
                                        sum(x == '.') == length(x)), ]
  message(paste('Removed', nVarsBefore - nrow(vcfGeno), 'out of', 
                nVarsBefore, 
                'variants because they were specific to w-/ore: \n',
                paste(setdiff(varsBefore, rownames(vcfGeno)),
                      collapse = '\n ')))
  vcfGR <- vcfGR[rownames(vcfGeno)]
  
  # remove variants where all of the samples have alternative genotype
  nVarsBefore <- nrow(vcfGeno)
  varsBefore <- rownames(vcfGeno)
  vcfGeno <- vcfGeno[!apply(vcfGeno, 1, 
                            function(x) sum(x == '1/1') + 
                              sum(x == '.') == length(x)), ]
  message(paste('Removed', nVarsBefore - nrow(vcfGeno), 'out of', 
                nVarsBefore, 'variants because all were ALT: \n',
                paste(setdiff(varsBefore, rownames(vcfGeno)),
                      collapse = '\n ')))
  vcfGR <- vcfGR[rownames(vcfGeno)]
  
  # Change 838 - 839, because every time 838 appears, 839 is also there
  if (("chrM:838_AT/A" %in% rownames(vcfGeno)) |
      ("chrM:839_T/A" %in% rownames(vcfGeno))) {
    message(paste('Automated check, that 838 and 839 co-occur passed:',
                  variantCoOccurence(vcfGeno[c("chrM:838_AT/A",
                                               "chrM:839_T/A"), ])))
    rownames(vcfGeno)[rownames(vcfGeno) == "chrM:838_AT/A"] <- "chrM:838_AT/TA"
    vcfGeno <- vcfGeno[rownames(vcfGeno) != "chrM:839_T/A", ]
    names(vcfGR)[names(vcfGR) == "chrM:838_AT/A"] <- "chrM:838_AT/TA"
    vcfGR <- vcfGR[names(vcfGR) != "chrM:839_T/A"]
  }
  
  # Change 5950 - 5975, this is intergenic region, it will be considered
  # separately. Now we will reduce all 4 variants found in it to one
  if (interGen == 'remove') {
    message('Removed intergenic region variants 5950 - 5975')
    vcfGR <- vcfGR[start(vcfGR) < 5960 | start(vcfGR) > 5975]
    vcfGeno <- vcfGeno[names(vcfGR), ]
  } 
  
  if (interGen == 'merge') {
    intergenVarsGeno <- vcfGeno[names(vcfGR[start(vcfGR) >= 5960 &
                                              start(vcfGR) < 5975]), ]
    intergenVarsMergeGeno <- mergeInterGenVars(intergenVarsGeno,
                                               'chrM:5960_ATATATTTATATATATATATATAT/TTA')
    vcfGeno <- vcfGeno[names(vcfGR[start(vcfGR) < 5960 |
                                     start(vcfGR) > 5975]), ]
    vcfGeno <- rbind(vcfGeno, intergenVarsMergeGeno)
    
    intergenVarsMergeGT <- vcfGR['chrM:5967_TTATATA/TTA']
    start(intergenVarsMergeGT) <- 5960
    end(intergenVarsMergeGT) <- 5975
    names(intergenVarsMergeGT) <- paste0('chrM:5960_ATATATTTATATATATATATATAT/TTA')
    vcfGR <- vcfGR[start(vcfGR) < 5960 | start(vcfGR) > 5975]
    vcfGR <- c(vcfGR, intergenVarsMergeGT)
    
    vcfGR <- vcfGR[order(start(vcfGR))]
    vcfGeno <- vcfGeno[names(vcfGR), ]
    message('Merged intergenic region variants 5950 - 5975 into one')
  }
  
  # because we removed w-, ore and [maybe] some DGRPs some variants might
  # change their status from MNPs to something else because trird alleles
  # were present just in ore/w-/removed DGRPs. We need to change alles in 
  # Granges too
  maybeMNPs <- vcfGR[getVariantStructType(vcfGR) == 'MNP']
  # 838 and 5960 are MNPs for sure
  maybeMNPs <- maybeMNPs[!names(maybeMNPs) %in% 
                         c('chrM:838_AT/TA', 
                           'chrM:5960_ATATATTTATATATATATATATAT/TTA')]
  if (length(maybeMNPs) != 0) {
    # number of "extra" (not 0/0, '.', '1/1', '0/1', '1/0' alleles in the set)
    if (is.vector(vcfGeno[names(maybeMNPs), ])) {
      numbOfAlleles <- length(unique(vcfGeno[names(maybeMNPs), ][grepl('[2-9]',
                                                vcfGeno[names(maybeMNPs), ])]))
    } else {
      numbOfAlleles <- apply(vcfGeno[names(maybeMNPs), ], 1,
                             function(x) length(unique(x[grepl('[2-9]', x)])))
    }

    # MNPs which became not-MNPs
    if (is.vector(vcfGeno[names(maybeMNPs), ])) {
      if(numbOfAlleles == 0) {notMNPs <- names(maybeMNPs)} 
      else {notMNPs <- NULL}
    } else {
      notMNPs <- names(numbOfAlleles[numbOfAlleles == 0])
    }
    if (length(notMNPs) != 0) {
      message(paste(length(notMNPs), 'variants: ', 
                    paste(notMNPs, collapse = ', '), 'changed their status from', 
                    'MNP to normal because of the removed lines'))
      vcfGR[notMNPs]$ALT <- DNAStringSetList(lapply(vcfGR[notMNPs]$ALT, 
                                                    function(x) x[1]))
    }
  }
  
  result <- list(vcfGR, vcfGeno)
  result
}

#' readInNuclVCF
#' Function to read-in VCF file for nuclear genomes
#' @param mitoVcfpath path to VCF file
#' @param genomeVer genome version, i.e. dm6
#' @param callRate minimul call rate to preseve variant in the set
readInNuclVCF <- function(vcfpath, genomeVer, callRate = NULL) {
  vcf <- readVcf(vcfpath, genomeVer)
  vcfGeno <- geno(vcf)$GT
  if (all(grepl('^DGRP-.*_DGRP-.*', colnames(vcfGeno)))) {
    colnames(vcfGeno) <- gsub('_.*', '', colnames(vcfGeno))
  }
  vcfGR <- rowRanges(vcf)
  
  # remove variants where all of the DGRPs are not genotyped
  nVarsBefore <- nrow(vcfGeno)
  vcfGeno <- vcfGeno[!apply(vcfGeno, 1, 
                            function(x) sum(x == '0/0') + 
                              sum(x == '.') == length(x)), ]
  message(paste('Removed', nVarsBefore - nrow(vcfGeno), 'out of', 
                nVarsBefore, 
                'variants because they were only with 0/0 and .'))
  
  if (!is.null(callRate)) {
    nVarsBefore <- nrow(vcfGeno)
    vcfGeno <- vcfGeno[!apply(vcfGeno, 1, 
                              function(x) sum(x == '.') <= (1 - callRate) * 
                                                           length(x)), ]
    message(paste('Removed', nVarsBefore - nrow(vcfGeno), 'out of', 
                  nVarsBefore, 'variants because call rate was < ', 
                  callRate))
  }
  
  vcfGR <- vcfGR[rownames(vcfGeno)]
  
  # modify IDs, so they correspond to the positions. They might not correspond 
  # to the position id vcf was moved from dm3 to dm6
  newNames <- paste(as.character(seqnames(vcfGR)), start(vcfGR), 
                    sapply(names(vcfGR), 
                           function(x) strsplit(x, '_')[[1]][3]),
                    sep = '_')
  names(vcfGR) <- newNames
  rownames(vcfGeno) <- newNames
   
  result <- list(vcfGR, vcfGeno)
  result
}
 
#' variantCoOccurence
#' Checks, if variants co-occur together. It assumes that there's <= 9 alleles
#' @param VOIgeno genotype table (as geno(vcf)$GT) for the variants of interest
#' @return TRUE, if variants co-occur together
variantCoOccurence <- function(VOIgeno) {
  coOcc <- all(apply(genoToRefHetAlt(VOIgeno), 1,
                     function(x) length(unique(x))) == 1)
  coOcc
}

# Variant annotation ----------------------------------------------------------
#' getAllFlyGenes
#' Returns table of all fly genes  ensemble identifiers and common names
#' @return data frame
#' @author Maria Litovchenko
#' @example  
#' library('biomaRt')
#' library('org.Dm.eg.db')
#' getAllFlyGenes()
getAllFlyGenes <- function() {
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset = "dmelanogaster_gene_ensembl", 
                    host = "jul2015.archive.ensembl.org")
  chroms = c('2L', '2R', '3L', '3R', '4', 'X', 'Y')
  egs <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
               filters = 'chromosome_name', values = chroms,
               ensembl)
  egs
}

#' getDeletedAAs
#' Returns number of deleted AAs from the protein as a result of variant
#' @param VOIanno variant of interest annotation line from annotation table
#' @param vcfGR GRanges for all the variants
#' @param chrMseq chrM sequence
#' @param mitoGenesTab table with the location of mitochondrial genes
getDeletedAAs <- function(VOIanno, vcfGR, chrMseq, mitoGenesTab) {
  message('WARNING: MNPs are not supported!')
  # start/end of the variant
  VOIgr <- vcfGR[VOIanno['Variant']]
  VOIstart <- start(VOIgr)
  VOIend <- VOIstart + nchar(as.character(VOIgr$REF))
  # get gene, where variant of interest falls
  VOIgene <- VOIanno['Gene']
  geneStart <- mitoGenesTab[VOIgene]$start
  geneEnd <- mitoGenesTab[VOIgene]$end
  # sequence of that gene
  geneSeq <- chrMseq[geneStart : geneEnd]
  # sequence with the variant
  geneSeqWithVar <- c(geneSeq[1 : (VOIstart - geneStart)],
                      tolower(unlist(strsplit(as.character(VOIgr$ALT[[1]]),
                                              '*'))),
                      geneSeq[(VOIend - geneStart + 1) : length(geneSeq)])
  
  geneStrand <- mitoGenes[VOIgene]$strand
  geneStrand <- ifelse(geneStrand == "-", 'R', 'F')

  originalProtein <- seqinr::translate(geneSeq, frame = 0, sens = geneStrand,
                                       numcode = 5)
  mutatedProtein <- seqinr::translate(geneSeqWithVar, frame = 0, 
                                      sens = geneStrand, numcode = 5)
  
  mutatedProtein <- mutatedProtein[1: min(which(mutatedProtein == '*'))]
  result <- length(originalProtein) - length(mutatedProtein)
  message(paste('Removes :', result))
  result
}

#' getInfoAboutGenes
#' Returns df with short info about genes: gene short name, full name and 
#' optional GO
#' @param martObj biomart object
#' @param geneNames ensembl IDs
#' @param withGO whatever or not add GO
#' @param withCoord whatever or not add Coordinates
#' @return data frame
getInfoAboutGenes <- function(martObj, geneNames, withGO = F, withCoord = F) {
  attrs <- c('ensembl_gene_id', 'external_gene_name')
  if (withGO) {
    attrs <- c(attrs, 'name_1006')
  } 
  if (withCoord) {
    attrs <- c('chromosome_name', 'start_position', 'end_position', 'strand', 
               attrs)
  }
  egs <- getBM(attributes = attrs, filters = 'external_gene_name',
               values = geneNames, martObj)
  if (withCoord) {
    result <- egs[, c('chromosome_name', 'start_position', 'end_position',
                      'ensembl_gene_id')]
    result <- cbind(egs, score = 0, strand = ifelse(egs$strand == 1, '+', '-'))
    if (withGO) {
      result <- cbind(result, name_1006 = egs$name_1006)
    }
  } else {
    result <- egs
  }
  result
}

#' getVariantStructType
#' Returns structural type of every variant
#' @param GRanges object with variants
#' @return structural type of each variants
getVariantStructType <- function(varsGR) {
  # length of the reference allele
  refLen <- varsGR$REF@ranges@width
  # number of the alternative alleles (determine MNPs)
  altNumbOfAll <- lapply(varsGR$ALT, length)
  # length of the alternative allele
  altLen <- lapply(varsGR$ALT, function(x) x@ranges@width)
  # if alternative allele is empty (for deletions)
  altEmpty <- sapply(1:length(varsGR), 
                     function(x) ifelse(altNumbOfAll[x] == 1, 
                                        unlist(varsGR[x]$ALT) == '',
                                        F))
  
  notOneAllele <- lapply(altNumbOfAll, function(x) x > 1)
  refLenMore1 <- lapply(refLen, function(x) x > 1)
  altLenMore1 <- lapply(altLen, function(x) x > 1)
  # because there could be more than one alternative allele (in case of MNP)
  # I replace such cases with FALSE to avoid problems
  altLenMore1[which(notOneAllele == T)] <- F
  altLenMore1 <- unlist(altLenMore1)
  altLen[which(notOneAllele == T)] <- 0
  altLen <- unlist(altLen)
  
  result <- rep('UNKNOWN', length(varsGR))
  # MNPs first
  result[unlist(notOneAllele)] <- 'MNP'
  result[result == 'UNKNOWN' & refLen == 1 & altLen == 1] <- 'SNP'
  result[result == 'UNKNOWN' & refLen != 1 & altLen == 1] <- 'DEL'
  result[result == 'UNKNOWN' & refLen == 1 & altEmpty == T] <- 'DEL'
  result[result == 'UNKNOWN' & refLen == 1 & altLen != 1] <- 'INS'
  result[result == 'UNKNOWN' & refLen != 1 & altLen != 1  & 
           refLen == altLen] <- 'SUBS'
  names(result) <- names(varsGR)
  
  # protection against BUGS!
  if (sum(result == 'UNKNOWN') != 0) {
    message('ERROR: detected unknown type of variant!')
    print(names(result)[which(result == 'UNKNOWN')])
    stop()
  }
  result
}

#' readChrMannot
#' Reads Manually curated GFT file for dm6, chrM
#' @param dm6ChrMmanualAnnoPath path to chrM annotation file
#' @return data.table with chr, start, end, strand, geneID, geneName, biotype
readChrMannot <- function(dm6ChrMmanualAnnoPath) {
  mitoGenes <- fread(dm6ChrMmanualAnnoPath)
  mitoGenes <- mitoGenes[V3 == 'gene']
  mitoGenes[, V9 := gsub("gene_id", "", V9)]
  mitoGenes[, V9 := gsub("gene_name", "", V9)]
  mitoGenes[, V9 := gsub("gene_biotype", "", V9)]
  mitoGenes[, V9 := gsub("\\s+", "", V9)]
  mitoGenes <- mitoGenes[, c(1, 4, 5, 7, 9), with = F]
  names(mitoGenes) <- c('chr', 'start', 'end', 'strand', 'V9')
  mitoGenes <- cbind(mitoGenes, matrix(unlist(strsplit(gsub("\"", "",
                                                            mitoGenes$V9), ';')),
                                       ncol = 4, byrow = T)[, -3])
  mitoGenes <- mitoGenes[, -5, with = F]
  names(mitoGenes) <- c('chr', 'start', 'end', 'strand', 'geneID', 'geneName', 
                        'biotype')
  setkey(mitoGenes, geneName)
  mitoGenes
}

#' readVariantAnnot
#' Reads-in variant annotation and converts it to the table
#' @param vcfpath path to annotated VCF
#' @param vcfGR GRanges object containing filtered variants for which 
#' annotation will be extracted
#' @return data table with Variant, StructType, FuncType, EffectSize, Gene,
#' AAchange
readVariantAnnot <- function(vcfpath, vcfGR) {
  # read functional annotation of variants
  vcfAnn <- readInfo(vcfpath, 'ANN')
  
  # check, if all the variants in GRanges match one in annotation
  varIsAnnot <- sapply(names(vcfGR), 
                       function(x) sum(grepl(x, names(vcfAnn))) != 0)
  # if we didn't find annotation for some variants
  if (sum(varIsAnnot) != length(vcfGR)) {
    message(paste('For variants', 
                  paste(names(varIsAnnot)[which(varIsAnnot == F)], 
                        collapse = ', '),
                  'annotation was not found, NAs will be added'))
  }
  vcfGRisAnnot <- vcfGR[varIsAnnot]
  # leave only annotation for variants which are in Granges
  varAnnWithGR <- vcfAnn[sapply(names(vcfGRisAnnot), 
                                function(x) grep(x, names(vcfAnn), value = T))]
  
  # reconstruct list to the table
  varAnnTab <- sapply(varAnnWithGR, function(x) strsplit(x, '\\|')[[1]])
  varAnnTab <- data.table(Variant = names(varAnnTab),
                          StructType = as.vector(getVariantStructType(vcfGRisAnnot)),
                          FuncType = sapply(varAnnTab, function(x) x[2]),
                          EffectSize = sapply(varAnnTab, function(x) x[3]),
                          Gene = sapply(varAnnTab, function(x) x[4]),
                          AAchange = sapply(varAnnTab, function(x) x[11]))
  
  # in readInMitoVcf we could delete certain lines, so structural type of the 
  # variants might also change from MNP to something else. However, in 
  # varAnnTab$Variant there will be comma in the name, indicating MNP, we need
  # to remove it
  removeComma <- varAnnTab[StructType != "MNP" & grepl(',', Variant)]
  if (nrow(removeComma) != 0) {
    ind <- which(varAnnTab$StructType != "MNP" & grepl(',', varAnnTab$Variant))
    varAnnTab$Variant[ind] <- gsub(',.*', '', removeComma$Variant)
  }
  
  if (!all(varIsAnnot == T)) {
    vcfGRnoAnno <- vcfGR[!varIsAnnot]
    toAdd <- data.table(Variant = names(vcfGRnoAnno), StructType = NA, 
                        FuncType = NA, EffectSize = NA, Gene = NA, 
                        AAchange = NA)
    varAnnTab <- rbind(varAnnTab, toAdd)
  }
  varAnnTab
}

#' GOenrichment
#' Performs GO enrichment of the selected gene set with TopGO
#' @param selectedGenes FlybaseIDs for the genes of interest
#' @param allGenesList background list of genes, selectedGenes should be
#' part of this list
#' @param ont ontology, BP, MF, CC
GOenricment <- function (selectedGenes, allGenesList, ont = 'BP',
                         topNodes = 10) {
  allGenesList_bg <- rep(1, length(allGenesList))
  names(allGenesList_bg) <- allGenesList
  allGenesList_bg[selectedGenes] <- 0.01
  
  tg.1 <- new("topGOdata", description = 'GO analysis',
              ontology =  ont, allGenes = allGenesList_bg,
              geneSel = function (x) {return (x < 0.05)},
              annot = annFUN.org ,
              nodeSize = 5 , # minimum number of genes in a GO categorie
              ID = "ENSEMBL", mapping = "org.Dm.eg.db")
  GO.res <- runTest(tg.1, algorithm = "elim", statistic = "fisher")
  
  result <- GenTable(tg.1, Fisher = GO.res, orderBy = "Fisher", 
                     ranksOf = "Fisher", topNodes = topNodes)
  result
}

#' readVariantAnnotNucl
#' Reads-in variant annotation and converts it to the table, nuclear variants
#' specific. There's a problem: annotation is given by ID of the variant, which
#' is in dm3 coordiantes, but I need to match it with actual variants, which 
#' are in dm6. Function takes care of it
#' @param vcfpath path to annotated VCF
#' @return data table with Variant, StructType, FuncType, EffectSize, Gene,
#' inGeneLoc, AAchange
readVariantAnnotNucl <- function(vcfpath) {
  # read functional annotation of variants
  vcfAnn <- readInfo(vcfpath, 'ANN')
  # reconstruct list to the table
  varAnnList <- lapply(vcfAnn, function(x) strsplit(x, '\\|'))
  varAnnTab <- data.table()
  for (i in 1:length(varAnnList)) {
    oneVarAnno <- varAnnList[[i]]
    toAdd <- data.table(Variant = rep(names(varAnnList)[i], 
                                      length(oneVarAnno)),
                        StructType = rep(strsplit(names(varAnnList)[i],
                                                  '_')[[1]][3], 
                                         length(oneVarAnno)),
                        EffectSize = sapply(oneVarAnno, function(x) x[3]),
                        Gene = sapply(oneVarAnno, function(x) x[4]),
                        inGeneLoc = sapply(oneVarAnno, function(x) x[2]),
                        AAchange = sapply(oneVarAnno, function(x) x[11]))
    varAnnTab <- rbind(varAnnTab, toAdd)
  }
  # there's a problem: annotation is given by ID of the variant, which is in
  # dm3 coordiantes, but I need to match it with actual GRDs, which are in dm6
  setnames(varAnnTab, 'Variant', 'dm3_ID')
  setkey(varAnnTab, 'dm3_ID')
  # get relationship between ID in dm3 and position in dm6
  dm3Todm6 <- rowRanges(readVcf(vcfpath))
  dm3Todm6 <- data.table(dm3_ID = names(dm3Todm6), 
                         NuclVar = paste0(seqnames(dm3Todm6), '_', 
                                          start(dm3Todm6)))
  setkey(dm3Todm6, 'dm3_ID')
  varAnnTab <- merge(dm3Todm6, varAnnTab)
  varAnnTab
}

#' simplifyAnnot
#' Simplifies annnotation table
#' @param varAnnTab annotation table from readVariantAnnot
#' @return simplified version
simplifyAnnot <- function(varAnnTab) {
  varAnnTab[, circosType := integer()]
  # let's consider variants by the biological gene type they affect
  # MNPs
  varAnnTab[StructType == 'MNP', 
            c('FuncType', 'EffectSize') := list('MNP', 'Unknown')]
  # tRNA
  varAnnTab[StructType == 'SNP' & grepl('tRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize', 
              'circosType') := list('SNP in tRNA', 'LOW', 3)]
  varAnnTab[(StructType == 'INS' | StructType == 'DEL') &  
            grepl('tRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize') := list('INDEL in tRNA', 'MODERATE')]
  # lrRNA
  varAnnTab[StructType == 'SNP' & grepl('lrRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize', 
              'circosType') := list('SNP in lrRNA', 'LOW', 4)]
  varAnnTab[(StructType == 'INS' | StructType == 'DEL') & 
             grepl('lrRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize') := list('INDEL in lrRNA', 'MODERATE')]
  # srRNA
  varAnnTab[StructType == 'SNP' &  grepl('srRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize', 
              'circosType') := list('SNP in srRNA', 'LOW', 5)]
  varAnnTab[(StructType == 'INS' | StructType == 'DEL') & 
            grepl('srRNA', Gene) & !grepl('-mt:', Gene),
            c('FuncType', 'EffectSize') := list('INDEL in srRNA', 'MODERATE')]
  # intergenic
  varAnnTab[StructType == 'SNP' & grepl('-mt:', Gene),
            c('FuncType', 'circosType') := list('intergenic SNP', 6)]
  varAnnTab[(StructType == 'INS' | StructType == 'DEL') &  grepl('-mt:', Gene),
            FuncType := 'intergenic INDEL']
  # in protein coding part
  varAnnTab[FuncType == 'synonymous_variant', 
            c('FuncType', 'circosType') := list('Synonymous', 1)]
  varAnnTab[FuncType == 'missense_variant', 
            c('FuncType', 'circosType') := list('Missense', 2)]
  varAnnTab[grepl('frameshift_variant', FuncType), FuncType := 'Frameshift']
  varAnnTab[FuncType == 'inframe_deletion', FuncType := 'Inframe DEL']
  varAnnTab[FuncType == 'inframe_insertion', FuncType := 'Inframe INS']
  varAnnTab[StructType == 'SNP' & grepl('stop_retained_variant', FuncType),
            c('FuncType', 'circosType') := list('Synonymous', 1)]
  varAnnTab
}

# Variant MAC/MAF/etc ---------------------------------------------------------
#' calcMacMaf
#' Calculates minor allele count (MAC) and minor allele frequnecy (MAF) for one
#' variant
#' @param oneVarFromGenoTab one line of GT from vcf
#' @return named vector with mac and maf 
calcMacMaf <- function(oneVarFromGenoTab) {
  # alternative allele count
  altCount <- length(oneVarFromGenoTab[oneVarFromGenoTab != '.' &
                                       oneVarFromGenoTab != '0/0'])
  # count also reference allele, to circumvent troubles with NAs
  refCount <- length(oneVarFromGenoTab[oneVarFromGenoTab == '0/0'])
  # minor allele count
  mac <- min(altCount, refCount)
  result <- c(mac, mac/(altCount + refCount))
  names(result) <- c('MAC', 'MAF')
  result
}

#' findVarClusters
#' Finds clusters (LD) of variants
#' @param vcfGenoTab genotype table
#' @return vector with names of the clusters and empty space for the 
#' clusterless variants
findVarClusters <- function(vcfGenoTab) {
  distMatr <- dist(vcfGenoTab)
  distMatr <- as.matrix(distMatr)
  # create list of clusters
  clusters <- list()
  for (j in 1:ncol(distMatr)) {
    if (!colnames(distMatr)[j] %in% unlist(clusters)) {
      newCluster <- rownames(distMatr)[which(distMatr[, j] == 0)]
      clusters[[length(clusters) + 1]] <- newCluster
    }
  }
  clusters <- clusters[sapply(clusters, function(x) length(x) > 1)]
  names(clusters) <- paste0('CL:', 1:length(clusters))
  # convert them to vector
  result <- rep("", nrow(vcfGenoTab))
  names(result) <- rownames(vcfGenoTab)
  for (varID in rownames(vcfGenoTab)) {
    inCluster <- sapply(clusters, function(x) ifelse(varID %in% x, T, F))
    result[varID] <- ifelse(sum(inCluster) != 0, 
                            paste0(names(clusters)[inCluster], ','), '')
  }
  result
}

# GRDs ------------------------------------------------------------------------
#' autoChiSq
#' Performs Fisher eexact test for detection of mito-nuclear incompatibility
#' @param nuclGeno VECTOR with genotype of nuclear variant
#' @param mitoGeno VECTOR with genotype of mitochondrial variant
#' @param cutoffOnLines minimal number of mito-nuclear combinations to compute
#' statistics
#' @param cutOffOnAllFreq minimal allele frequency to compute statistics
#' @return 
autoChiSq <- function(nuclGeno, mitoGeno, cutoffOnLines, cutOffOnAllFreq) {
  mitoNuclDf <- cbind(nuclGeno, mitoGeno)
  mitoNuclDf <- mitoNuclDf[apply(mitoNuclDf, 1, 
                                 function(x) all(x != "0" & x != "2")), ]
  if (!is.null(mitoNuclDf) & 
      (is.data.frame(mitoNuclDf) | is.matrix(mitoNuclDf))) {
    if (nrow(mitoNuclDf) >= cutoffOnLines) {
      mitoNuclDf <- as.data.frame(mitoNuclDf)
      if ((min(table(mitoNuclDf[, 1])) > cutOffOnAllFreq * nrow(mitoNuclDf)) &
          (min(table(mitoNuclDf[, 2])) > cutOffOnAllFreq * nrow(mitoNuclDf))) {
        # contigency table
        contTab <- table(mitoNuclDf)
        if (ncol(contTab) == 2 & nrow(contTab) == 2) {
          chiSq <- chisq.test(contTab, correct = F)
          result <- c(chiSq$statistic, chiSq$`p.value`, min(contTab) == 0)
        } else {
          result <- rep(NA, 3)
        }
      } else {
        result <- rep(NA, 3)
      }
    } else {
      result <- rep(NA, 3)
    }
  } else {
    result <- rep(NA, 3)
  }
  result
}

#' contigencyTable
#' Returns contigency table for nuclear - mitochondrial variant pair
#' @param nuclGeno nuclear geno string
#' @param mitoGeno mitochondrial geno string
#' @return contingency table
contigencyTable <- function(nuclGeno, mitoGeno, cutoffOnLines = 150) {
  mitoNuclDf <- cbind(nuclGeno, mitoGeno)
  mitoNuclDf <- mitoNuclDf[apply(mitoNuclDf, 1, 
                                 function(x) all(x != "0" & x != "2")), ]
  # contigency table
  contTab <- table(as.data.frame(mitoNuclDf))
  if (nrow(mitoNuclDf) < cutoffOnLines) {
    message('Number of samples do not pass cut off on samples')
  }
  contTab
}

#' GRDtabToGranges
#' Coverts data table with SIGNIFICANT GRDs to Granges
#' @param GRDdataTab data table with SIGNIFICANT GRDs
#' @return GRanges object
GRDtabToGranges <- function(GRDdataTab) {
  # convert it to GRanges object 
  GRDdataTab[, NuclVar := as.character(NuclVar)]
  GRDdataTab[, chr := sapply(GRDdataTab$NuclVar, 
                             function(x) strsplit(x, '_')[[1]][1])]
  GRDdataTab[, start := as.integer(sapply(GRDdataTab$NuclVar, 
                               function(x) strsplit(x, '_')[[1]][2]))]
  GRDdataTab[, end := as.integer(sapply(GRDdataTab$NuclVar, 
                                 function(x) strsplit(x, '_')[[1]][2]))]
  GRD_gr <- makeGRangesFromDataFrame(GRDdataTab, keep.extra.columns = T)
  GRD_gr <- sort(sortSeqlevels(GRD_gr))
  GRD_gr
}

plotGRDbyVar <- function(GRDs_VartoPlot, labelCol = NULL, ...) {
  if (grepl('Mito', colnames(GRDs_VartoPlot)[1])) {
    result <-  ggplot(GRDs_VartoPlot, aes(x = MitoVar, y = N)) +
               geom_bar(colour = pallete(5)[5], stat = "identity", 
                        fill = pallete(5)[5]) + coord_flip() +
               xlab('mtDNA variant') + 
               ylab('# of incompatible nuclear variants') +
               mashaGgplot2Theme + ...
    if (!is.null(labelCol)) {
      result <- result + theme(axis.text.y = element_text(colour = labelCol))
    }
  } else {
    result <- ggplot(GRDs_VartoPlot, aes(x = 1:nrow(GRDs_VartoPlot), y = N)) +
              geom_bar(colour = pallete(5)[5], stat = "identity", 
              fill = pallete(5)[5]) + coord_flip() + xlab('Nuclear variant') +
              ylab('# of incompatible mitochondial variants')  + 
              mashaGgplot2Theme + 
              theme(axis.text.y = element_blank(),
                    axis.ticks = element_blank()) + ...
  }
  result
}

#' selectVarsWithSignNeighbors
#' @param GRDdataTab data table with SIGNIFICANT GRDs
#' @return Granges object with nuclear variants which have significant GRD 
#' and have at least one variant from each side closer than distance which is
#' also significant GRD
selectVarsWithSignNeighbors <- function(GRDdataTab, LDdist) {
  GRD_gr <- GRDtabToGranges(GRDdataTab)
  GRD_gr <- reduce(GRD_gr)
  # proceed into chromosomes, because otherwise there will be suspicious calls
  GRD_grSignNeigh <- GRanges()
  print(paste('Started calculating neighbors at', Sys.time()))
  for (chrom in seqlevels(GRD_gr)) {
    # get all starts and calculate difference to the neighbors
    # GRDdist_p1 - distance to the first preceeding variant
    # GRDdist_p2 - distance to the second preceeding variant
    chromGRD_gr <- GRD_gr[seqnames(GRD_gr) == chrom]
    if (length(chromGRD_gr) >= 3) {
      GRDpos <- start(chromGRD_gr)
      GRDposRev <- rev(GRDpos)
      GRD_gr_coordsDist <- data.frame(GRDdist_p1 = c(NA, diff(GRDpos, 1)),
                                      GRDdist_p2 = c(NA, NA, diff(GRDpos, 2)),
                                      GRDdist_p3 = c(NA, NA, NA,
                                                     diff(GRDpos, 3)),
                                      GRDdist_f1 = -c(rev(diff(GRDposRev, 1)),
                                                      NA),
                                      GRDdist_f2 = -c(rev(diff(GRDposRev, 2)), 
                                                      NA, NA),
                                      GRDdist_f3 = -c(rev(diff(GRDposRev, 3)),
                                                      NA, NA, NA))
      GRD_gr_toTake <- apply(GRD_gr_coordsDist, 1, 
                             function(x) sum(x <= LDdist) >= 2)
      GRD_gr_toTake[is.na(GRD_gr_toTake)] <- F
      GRD_grSignNeigh <- c(GRD_grSignNeigh, chromGRD_gr[GRD_gr_toTake])
    }
  }
  print(paste('Finished calculating neighbors at', Sys.time()))
  GRD_grSignNeigh <- paste0(seqnames(GRD_grSignNeigh), '_', 
                            start(GRD_grSignNeigh))
  GRD_grSignNeigh
}

writeGRDsCircos <- function(GRD_Granges, outputPath) {
  if (class(GRD_Granges) != 'data.table') {
    GRDs_dt <- as.data.table(mcols(GRD_Granges))
  } else {
    GRDs_dt <- GRD_Granges
  }
  # position of incompatibilities on mitochondrial side
  mitoPos <- sapply(GRDs_dt$MitoVar, 
                    function(x) strsplit(gsub('chrM:', '', x), '_')[[1]][1])
  # position of incompatibilities on nuclear side
  nuclChr <- sapply(GRDs_dt$NuclVar, 
                    function(x) gsub('chr', '', strsplit(x, '_')[[1]][1]))
  nuclChr <- paste0('dm', tolower(nuclChr))
  nuclChr <- gsub('dmx', 'dmX', nuclChr) 
  nuclPos <- sapply(GRDs_dt$NuclVar, function(x) strsplit(x, '_')[[1]][2])
  result <- data.frame(mitoChr = 'chrM', mitoStart = mitoPos, 
                       mitoStop = mitoPos, nuclChr, nuclStart = nuclPos,
                       nuclStop = nuclPos)
  write.table(result, outputPath, quote = F, sep = ' ', row.names = F,
              col.names = F)
}

# Plotting functions ----------------------------------------------------------
#' plotVarsPerDGRP
#' Plots number of study-specific variants per DGRP line
#' @param studySpecVcfGeno geno table with study-specific variants
#' @param studyName name of the study
#' @return plot
plotVarsPerDGRP <- function(studySpecVcfGeno, studyName) {
  studySpecDGRPvarCount <- sort(apply(studySpecVcfGeno, 2, 
                                      function(x) sum(x != '.' & x != '0/0')))
  studySpecDGRPvarCount <- as.data.frame(studySpecDGRPvarCount)
  colnames(studySpecDGRPvarCount) <- 'count'
  dgrpNames <- gsub('DGRP-', '', rownames(studySpecDGRPvarCount))
  dgrpOrder <- order(studySpecDGRPvarCount$count)
  studySpecDGRPvarCount$DGRP <- factor(dgrpNames, levels = dgrpNames[dgrpOrder])
  ggplot(studySpecDGRPvarCount, aes(x = DGRP, y = count)) + 
    geom_bar(stat = "identity", col = 'black') + mashaGgplot2Theme +
    ylab('Variant count') + ggtitle(paste('Number of', studyName, 
                                          'specific variants per DGRP line')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
}

#' plotMACPerStudySpecVar
#' Plots MAF of study-specific variants
#' @param studySpecVcfGeno geno table with study-specific variants
#' @param studyName name of the study
#' @return plot
plotMACPerStudySpecVar <- function(studySpecVcfGeno, studyName) {
  if (grepl('/', studySpecVcfGeno[1, 1])) { # in case of vcf-formated
    studySpecMAC <- sort(apply(studySpecVcfGeno, 1, 
                               function(x) sum(x != '.' & x != '0/0')))
    studySpecMAC <- as.data.frame(studySpecMAC)
  } else { # in case of richardson
    studySpecMAC <- apply(studySpecVcfGeno, 1,  function(x) sort(table(x))[1])
    studySpecMAC <- as.data.frame(studySpecMAC)
  }
  colnames(studySpecMAC) <- 'count'
  dgrpNames <- gsub('DGRP-', '', rownames(studySpecMAC))
  dgrpOrder <- order(studySpecMAC$count)
  studySpecMAC$DGRP <- factor(dgrpNames, levels = dgrpNames[dgrpOrder])
  ggplot(studySpecMAC, aes(x = DGRP, y = count)) + 
    geom_bar(stat = "identity", col = 'black') + mashaGgplot2Theme +
    ylab('# of DGRPs') + xlab('Variant') +
    ggtitle(paste('MAC for', studyName, 'specific variants')) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
}

#' getReadSupportOfAlt
#' Returns read support of alternative allele for study-specific variants
#' @param studySpecVcfGeno geno table for study specific variants
#' @param studyVcf vcf object for that study
#' @return data.frame with columns DGRP, Variant, ReadCount, ReadCountPerc
getReadSupportOfAlt <- function(studySpecVcfGeno, studyAD) {
  result <- data.frame(DGRP = character(), Variant = character(), 
                       ReadCount = integer(),  ReadCountPerc = numeric())
  for (i in 1:nrow(studySpecVcfGeno)) {
    for (j in 1:ncol(studySpecVcfGeno)) {
      if (studySpecVcfGeno[i, j] != '.') {
        readCount <- studyAD[rownames(studySpecVcfGeno)[i],
                             colnames(studySpecVcfGeno)[j]]
        for (k in 1:length(readCount)) {
          toAdd <- data.frame(DGRP = colnames(studySpecVcfGeno)[j],
                              Variant = rownames(studySpecVcfGeno)[i],
                              ReadCount = readCount[[k]][2],
                              ReadCountPerc = 100 * readCount[[k]][2] / 
                                sum(readCount[[k]]))
          result <- rbind(result, toAdd)
        }
      }
    }
  }
  result
}

#' plotPhenoGeno
#' Plots GWAS classical box plots
#' @param genoOI genotype of interest, character, should be in the same order
#' as phenoOI
#' @param genoName genotype name
#' @param phenoOI phenotype of interest should be in the same order as genoOI
#' @param phenoName name of the phenotype
#' @param pval p value
#' @return ggplot2
plotPhenoGeno <- function(genoOI, genoName, phenoOI, phenoName, pval) {
  dtToPlot <- data.frame(Genotype = genoOI, Phenotype = phenoOI)
  dtToPlot <- dtToPlot[complete.cases(dtToPlot), ]
  colnames(dtToPlot) <- c('Genotype', 'Phenotype')
  colorsToPlot <- colorRampPalette(c("#BB4444", "#4477AA"))
  colorsToPlot <- colorsToPlot(length(unique(dtToPlot$Genotype)))
  print(ggplot(dtToPlot, aes(x = Genotype, y = Phenotype, fill = Genotype)) +
          geom_boxplot() + geom_jitter(position = position_jitter(0.1)) + 
          mashaGgplot2Theme + ylab(phenoName) +
          ggtitle(paste('GWAS association between', phenoName, '\nand',
                        genoName, 'p =', round(pval, 5))) +
          scale_x_discrete(labels = sort(unique(dtToPlot$Genotype))) +
          scale_fill_manual(values = colorsToPlot) + 
          guides(fill = F))
}

# GWAS ------------------------------------------------------------------------
#' adjustForWolbInsBRB
#' Adjusts expression for the presence of wolbachia, insertions and BRB library
#' effect. Only for NOT factorial phenotypes
#' @param geneExpr data table with column ID and another column with expression
#'                 of ONE gene OR phenotype
#' @param covars data table with columns ID and one column per covariate,
#'               covariats as factors if they are factors
#' @return expression matrix with corrected for wolbachia, insertions and brb
adjustForWolbAndIns <- function(geneExpr, covars) {
  # convert both gene expression and covariats to data table for merging
  geneExprDT <- as.data.table(geneExpr)
  setnames(geneExprDT, colnames(geneExprDT), c('ID','geneExpr'))
  setkey(geneExprDT, ID)
  setkey(covars, ID)
  df <- merge(geneExprDT, covars)
  # we need to remove NAs in order to perform adjustment, we will introduce 
  # them later
  dfNoNA <- df[complete.cases(df), ]
  baseFormula <- as.formula("geneExpr ~ dummy")
  addedCovars <- 0
  
  if (length(unique(dfNoNA$In_2L_t)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + In_2L_t)
    addedCovars <- addedCovars + 1
  }
  if (length(unique(dfNoNA$In_2R_NS)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + In_2R_NS)
    addedCovars <- addedCovars + 1
  }
  if (length(unique(dfNoNA$In_3R_K)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + In_3R_K)
    addedCovars <- addedCovars + 1
  }
  if (length(unique(dfNoNA$In_3R_P)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + In_3R_P)
    addedCovars <- addedCovars + 1
  }
  if (length(unique(dfNoNA$In_3R_Mo)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + In_3R_Mo)
    addedCovars <- addedCovars + 1
  }
  if (length(unique(dfNoNA$Wolb)) >= 2) {
    baseFormula <- update(baseFormula,    ~ . + Wolb)
    addedCovars <- addedCovars + 1
  }
  
  result <- geneExprDT
  if (addedCovars != 0) {
    baseFormula <- update(baseFormula,    ~ . - dummy)
    lmGE <- lm(baseFormula, dfNoNA)
    result <- resid(lmGE) + coef(lmGE)[1]
    result <- data.frame(ID = dfNoNA$ID,
                         adjGeneExpr = result)
  }
  setnames(result, colnames(result), c('ID', 'adjGeneExpr'))
  
  # re-introduce NAs
  dfNA <- df[!complete.cases(df), c('ID', 'geneExpr')]
  setnames(dfNA, colnames(dfNA), c('ID', 'adjGeneExpr'))
  result <- rbind(result, dfNA)
  # return order to the original one
  setkey(result, ID)
  result <- result[geneExpr$ID, ]
  setnames(result, colnames(result), colnames(geneExpr))
  result
}

#' allelesMAFunderPheno
#' Checks, if one variant passes minimum MAC cutoff under phenotype
#' @param phenoGenoDT data frame with phenotype and genotype
#' @param minMAC cut off on minor allele count
#' @return true or false
allelesMAFunderPheno <- function(phenoGenoDT, minMAC) {
  phenoGenoDT <- phenoGenoDT[complete.cases(phenoGenoDT), ]
  min(table(phenoGenoDT[, 2])) >= minMAC
}

#' getPvalMatr
#' @param gwasRes results of gwas, data table
#' @param pvalName name of the column with significance (P or EMP1)
#' @param varName name of the column with variant (SNP by default)
#' @return matrix with P values, with phenotypes in columns and
#'         variants in rows
getPvalMatr <- function(gwasRes, pvalName, varName = 'SNP') {
  pvalMatr <- lapply(unique(gwasRes$PhenoID), 
                     function(x) gwasRes[PhenoID == x, pvalName, with = F])
  pvalMatr <- do.call(cbind, pvalMatr)
  pvalMatr <- as.data.frame(pvalMatr)
  colnames(pvalMatr) <- unique(gwasRes$PhenoID)
  rownames(pvalMatr) <- unlist(unique(gwasRes[, varName, with = F]))
  pvalMatr
}

#' KStest
#' @param phenoName name of the phenotype
#' @param nuclDT data table with results for nuclear variants
#' @param mitoDT data table with results for mitochondrial variants
#' @param ground on what to perform the test (name of the column)
#'               can be STAT, P, EMP1
#' @param useAbs whatever or not absolute values of test statistics should be 
#'               used
#' @return vector with results of test
KStest <- function(phenoName, nuclDT, mitoDT, ground = 'STAT', useAbs = T,
                   ...) {
  # nuclear distribution
  nuclDist <- unlist(nuclDT[PhenoID == phenoName][, ground, with = F])
  if (useAbs) {nuclDist <- abs(nuclDist)}
  # mitochondrial distribution
  mitoDist <- unlist(mitoDT[PhenoID == phenoName][, ground, with = F])
  if (useAbs) {mitoDist <- abs(mitoDist)}
  
  # do the test
  ksOnePhenoSTAT <- ks.test(nuclDist, mitoDist, ...)
  res <- c(phenoName, ksOnePhenoSTAT$statistic, ksOnePhenoSTAT$p.value,
           ground)
  res
}

#' PhenofromDF
#' This function creates PLINK-suitable PHENO file from data frame. 
#' Data frame should have future individual id in the first column and all the
#' rest columns should be phenotypes. Matching will be performed inside the 
#' function
#' @param phenotypes data frame with phenotypes, first column ID, 
#'                   others - phenos
#' @param adjusted indicates, whatever adjustment for covariats should be made
#' @param covars data table with columns ID and one column per covariate
#' @return data frame with Fam_ID, Ind_Id and phenotypes
PhenofromDF <- function(phenotypes, adjusted = F, covars = data.frame()) {
  result <- data.frame(Fam_ID = phenotypes$ID, Ind_ID = phenotypes$ID)
  phenoNames <- colnames(phenotypes)
  phenoNames <- phenoNames[phenoNames != 'ID']
  if (adjusted == T & nrow(covars) != 0) {
    for (phenoName in phenoNames) {
      onePheno <- data.table(ID = phenotypes$ID,
                             Pheno = phenotypes[, phenoName])
      onePheno[, Pheno := suppressWarnings(as.numeric(as.character(Pheno)))]
      result <- cbind(result, adjustForWolbAndIns(onePheno, 
                                                  covars)[, 2, with = F])
      print(paste(phenoName, dim(result)))
    }
    colnames(result) <- c("FID", "IID", phenoNames)
  } else {
    message('Adjustment for covariats was not performed')
    result <- cbind(result, phenotypes[, phenoNames])
  }
  result
}

#' phenoGenoPassMAC 
#' Checks, if phenotype - genotype combinations pass MAC cutoff
#' @param genoDT data table of genotypes
#' @param phenoDT data table of phenotypes
#' @param macCut cut off on minor allele count
#' @param genoGeno whatever or not a table for genotype-genotype combinations
#'        is being calculated
#' @return matrix with TRUE/FALSE values
phenoGenoPassMAC <- function(genoDT, phenoDT, macCut, genoGeno = F) {
  # arrange DGRPs in the same order for both genotypes and phenotypes
  genoDT <- genoDT[order(DGRP), ]
  phenosDT <- phenoDT[order(DGRP), ]
  # all combinations of genotypes and phenotypes
  allComb <- expand.grid(colnames(phenosDT)[-1], 
                         colnames(genoDT)[-ncol(genoDT)])
  if (genoGeno) {
    allComb <- expand.grid(colnames(phenosDT)[-ncol(phenosDT)], 
                           colnames(genoDT)[-ncol(genoDT)])
  }
  passMAC <- apply(allComb, 1, 
                   function(x) allelesMAFunderPheno(data.frame(phenosDT[, x[1], 
                                                                        with = F],
                                                               genoDT[, x[2],
                                                                      with = F]),
                                                    macCut))
  passMAC <- t(matrix(passMAC, nrow = ncol(phenosDT) - 1, 
                      ncol = ncol(genoDT) - 1, byrow = F))
  colnames(passMAC) <- colnames(phenoDT)[-1]
  rownames(passMAC) <- colnames(genoDT)[-ncol(genoDT)]
  passMAC
}

#' readPlinkOutput
#' Reads output of PLINK GWAS 
#' @param pathToFile path to PLINK output file
#' @param phenotypeDB phenotype data table with columns: Fam_ID, Ind_ID, 
#' phenotype1, phenotype2, etc
#' @return data.table with columns CHR, SNP, EMP1,EMP2, PhenoID in case study 
#' was done with permutations and columns CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P,
#' PhenoID in case it's just linear
readPlinkOutput <- function(pathToFile, phenotypeDB) {
  gwasResult <- fread(pathToFile, header = T)
  # phenotypeIndex <- strsplit(pathToFile, '\\.')[[1]]
  # phenotypeIndex <- phenotypeIndex[grepl('^P\\d+$', phenotypeIndex)]
  # phenotypeIndex <- as.integer(gsub('P', '', phenotypeIndex)) + 2
  # gwasResult[, PlinkPhenoID := as.integer(gsub('P', '', phenotypeIndex))]
  # gwasResult[, PhenoID := colnames(phenotypeDB)[phenotypeIndex]]
  
  phenotypeIndex <- strsplit(pathToFile, '\\.')[[1]][2]
  gwasResult[, PlinkPhenoID := phenotypeIndex]
  gwasResult[, PhenoID := phenotypeIndex]
  
  gwasResult
}

#' readPlinkOutFolder
#' Reads and summarizes in one table output from plink for many phenotypes
#' @param dirPath path to the directory with all PLINK results
#' @param pattern pattern of the files names
#' @param phenotypeDB data base of phenotypes
#' @param clustPhenos cluster of phenotypes
#' @param macPass data table which tells if MAC cutoff is passed for the 
#'                each combination of phenotype and genotype
#' @param correctIDs if IDs need to be corrected (only nuclear)
readPlinkOutFolder <- function(dirPath, pattern, phenotypeDB, clustPhenos, 
                               macPass = NULL, correctIDs = F) {
  # read - in all results files
  allGwasPaths <- list.files(dirPath, pattern = pattern, full.names = T)
  gwas <- lapply(allGwasPaths, readPlinkOutput, phenotypeDB)
  gwas <- do.call("rbind", gwas)
  
  if (correctIDs) {
    # during the reading in of vcf file used for the check of phenotype-
    # genotype MAF (see readInNuclVCF) I replaced wrong bp in ID to the
    # right one (because of dm3 -> dm6 transition), i do here the same
    newNames <- do.call(rbind, sapply(gwas$SNP, function(x) strsplit(x, '_')))
    newNames[, 1] <- sapply(newNames[, 1], 
                            function(x) ifelse(grepl('^4', x), x,
                                               paste0('chr', x)))
    newNames <- paste(newNames[, 1], gwas$BP, newNames[, 3], sep = '_')
    gwas$SNP <- newNames
  }
  
  # this line of code does nothing if GWAS was performed on all phenotypes and
  # it selects for independent phenotypes if it wasn't
  gwas <- gwas[PhenoID %in% colnames(phenotypeDB), ]
  if (!is.null(macPass)) {
    if ("P" %in% colnames(gwas)) {
      gwas$P <- apply(gwas, 1, function(x) ifelse(macPass[x['SNP'],
                                                          x['PhenoID']], 
                                                  as.numeric(x['P']), "noDATA"))
    }
    if ("EMP1" %in% colnames(gwas)) {
      gwas$EMP1 <- apply(gwas, 1, function(x) ifelse(macPass[x['SNP'],
                                                             x['PhenoID']], 
                                                     as.numeric(x['EMP1']), 
                                                     "noDATA"))
    }
  }
  # add phenotype cluster
  gwas[, PhenoClust := clustPhenos[PhenoID]$cluster]
  gwas
}

#' selectRandomNuclVars
#' Selects random variants 
#' @param varsPosGR granges object with positions of ALL variants to select 
#'                  from
#' @param numbOfVars number of random variants to select
#' @param minDist minimal distance between selected variants
#' @return GRanges object with selected random variants 
selectRandomNuclVars <- function(varsPosGR, numbOfVars, minDist) {
  minVarDist <- 0
  while (minVarDist < minDist) {
    # get random variants
    randomVars <- sample(1:length(varsPosGR), numbOfVars)
    randomVars <- varsPosGR[randomVars]
    # calculate the distance between them
    varDist <- sapply(1:numbOfVars, 
                      function(x) distanceToNearest(randomVars[x],
                                                    randomVars[-x]))
    varDist <- sapply(varDist, function(x) x@elementMetadata@listData$distance)
    minVarDist <- min(unlist(varDist))
  }
  randomVars
}

#' showGWASphenoStats
#' Prints out statistics about phenotypes in GWAS
#' @param gwasRes results of GWAS
#' @param pvalName name of the field which contains p value
#' @param phenoInfoDT data table with phenotypes
#' @return VOOOOOID!
showGWASphenoStats <- function(gwasRes, pvalName, phenoInfoDT) {
  # add phenotypes full names
  gwasRes$PhenoFullName <- phenoInfoDT[gwasRes$PhenoID]$FullPhenotype
  gwasResSign <- gwasRes[unlist(gwasRes[, pvalName, with = F]) < 0.05]
  
  # select significant phenotypes
  signPhenos <- unique(gwasResSign$PhenoID)
  signPhenosFullName <- unique(gwasResSign$PhenoFullName)
  message(paste('Number of independent phenotypes with association < 0.05:',
                length(unique(gwasResSign$PhenoID))))
  message(paste('Their names:', paste(signPhenosFullName, collapse = ', ')))
  
  # number of hits per phenotype:
  hitsPerPheno <- table(gwasResSign$PhenoID)
  hitsPerPhenoTab <- table(table(gwasResSign$PhenoID))
  hitsPerPhenoTab <- data.frame(hitsPerPhenoTab)
  message('Number of hits per phenotype:')
  print(hitsPerPhenoTab)
  
  # phenotype with the most number of hits
  freqPheno <- names(hitsPerPheno)[which(hitsPerPheno == max(hitsPerPheno))]
  freqPheno <- phenoInfoDT[freqPheno]$FullPhenotype
  message(paste0('Phenotype with maximum number of hits (', max(hitsPerPheno), ')',
                 ' is ', freqPheno))
  
  # phenotype with the second most number of hits
  secondFreq <- as.integer(hitsPerPhenoTab[nrow(hitsPerPhenoTab) - 1, 1])
  secondFreqPheno <- names(hitsPerPheno[which(hitsPerPheno == secondFreq)])
  secondFreqPheno <- phenoInfoDT[secondFreqPheno]$FullPhenotype
  message(paste0('Phenotype with second max number of hits (', secondFreq,
                 ') is ', paste(secondFreqPheno, collapse = ', ')))
  
  # correlated phenotypes
  corrSignPhenosClust <- unique(gwasResSign$PhenoClust)
  corrSignPhenos <- phenoClusters[cluster %in% corrSignPhenosClust]$phenotype
  message(paste0('Correlated phenotypes ', 
                 paste(phenoInfoDT[corrSignPhenos]$FullPhenotype,
                       collapse = ', ')))
  
  result <- c(paste(signPhenosFullName, collapse = ', '),
              paste(phenoInfoDT[corrSignPhenos]$FullPhenotype,
                    collapse = ', '))
  result
}

#' simes.test
#' Performs simes test on all p-values
#' @param x vector of p-values
#' @return p-value
simes.test <- function (x) {
  x <- na.omit(x)
  r <- rank(x)
  T <- min(length(x) * x/r)
  T <- ifelse(T < 1, T, 1)
  T
}

#' pAdjMultPhenBenjamini2013
#' Performs adjustment of p-values and selection of the significant 
#' associations for multiple phenotypes GWAS/tests according to
#' Selective inference on multiple families of hypotheses by
#' Yoav Benjamini & Marina Bogomolov 2013
#' @param pvalMatr matrix of RAW p values, where SNPs are in rows and 
#'                 phenotypes are in columns. SNPs IDs are rownames,
#'                 phenotype IDs are column names
#' @param adjMeth adjustment method to use: "holm", "hochberg", "hommel", 
#'                "bonferroni", "BH", "BY", "fdr", "none". Default : "BH"
#' @param pvalCut cut off on adjusted p values to apply. Default : 0.05
#' @return data frame with 4 columns Vars, Pheno, pvalue, pvalueAdj
pAdjMultPhenBenjamini2013 <- function(pvalMatr, adjMeth = 'BH', 
                                      pvalCut = 0.05) {
  # P - number of phenotypes, S - number of snps
  # Step 2: calculate the intersection hypothesis for each SNP by combining all 
  # the P p-values calculated for that SNP by using Simes's test
  intersectHyp <- apply(pvalMatr, 1, simes.test)
  # Step 3: test all the S intersection hypotheses by using the BH procedure
  # at level pvalCut by using the p-values calculated at step 2. 
  # Let R be the number of rejected intersection hypotheses.
  intersectHypCorr <- p.adjust(intersectHyp, method = adjMeth)
  intersectHypCorrPass <- which(intersectHypCorr <= pvalCut)
  R <- length(intersectHypCorrPass)
  message(paste0('R = ', R))
  if (R != 0) {
    # Step 4: select the R SNPs for which the intersection hypothesis was 
    # rejected at stage 3. 
    pvalMatr_R <- pvalMatr[intersectHypCorrPass, ]
    # For each selected SNP, apply the BH procedure at level R×0.05/S on all the 
    # P p-values calculated at stage 1.
    pvalCutoff <- R * pvalCut / nrow(pvalMatr) # level of p-value
    message(paste('P value cutoff calculated at stage 3 is', pvalCutoff))
    pvalMatrAdj <- apply(pvalMatr_R, 1, p.adjust, adjMeth)
    # now SNPs are in columns, and phenotypes are in rows
    res <- apply(pvalMatrAdj, 2, function(x) which(x <= pvalCutoff))
    res <- res[sapply(res, function(x) length(x) > 0)]
    if (!is.null(res)) {
      if (is.list(res)) {
        res <- lapply(names(res), 
                      function(x) data.frame(Vars = x,
                                             Pheno = names(res[[x]]),
                                             pvalue = unlist(pvalMatr_R[x, names(res[[x]])]),
                                             pvalueAdj = pvalMatrAdj[names(res[[x]]), x]))
        res <- do.call(rbind, res)
        rownames(res) <- NULL
      }
      if (is.vector(res)) {
        res <- lapply(names(res), 
                      function(x) data.frame(Vars = names(res[x]),
                                             Pheno = colnames(pvalMatr)[res[x]],
                                             pvalue = pvalMatr_R[names(res[x]),
                                                                 colnames(pvalMatr)[res[x]]],
                                             pvalueAdj = pvalMatrAdj[colnames(pvalMatr)[res[x]],
                                                                     names(res[x])]))
        res <- do.call(rbind, res)
        rownames(res) <- NULL
      }
    } 
  } else {
    res <- NULL
    message('No intersection hyposthesis pass selected p-value cutoff')
  }
  res
}

#' writePlinkRandomVarsScript
#' Writes PLINK script to run PLINK on random variants
#' @param varList list of GRanges containing variants
#' @param outDir output directory, there to put all the files and all the 
#'               script
#' @param baseVCF vcf from which to extract variants, not gz
#' @param phenoFile path to phenotype file
#' @param scriptName name of the output script
writePlinkRandomVarsScript <- function(varList, outDir, baseVCF, phenoFile,
                                       scriptName) {
  # init the script
  write("#!/bin/bash", scriptName, append = F)
  for (i in 1:length(varList)) {
    oneRandVars <- varList[[i]]
    # write down a positions file for the extraction from base vcf
    oneRandVarsDF <- as.data.frame(oneRandVars)[, 1:3]
    write.table(oneRandVarsDF, paste0('random_variants_', i, '.position'),
                col.names = F, row.names = F, quote = F, sep = '\t')
    # command to extract random positions from base vcf
    extrRandVars <- paste0('vcftools --vcf ', baseVCF, ' --positions ', 
                           'random_variants_', i, '.position --recode',
                           ' --recode-INFO-all --out ', outDir, 
                           'random_variants_', i)
    # translate to PLINK format
    trPlink1 <- paste0("plink --vcf ", outDir, "random_variants_", i, 
                       ".recode.vcf --allow-extra-chr --out ", outDir, 
                       "random_variants_", i)
    trPlink2 <- paste0("vcftools --vcf ", outDir, "random_variants_", i, 
                       ".recode.vcf --plink --out ", outDir, 
                       "random_variants_", i)
    # run PLINK
    makeOutDir <- paste0('outputdir=results/PLINK/PLINK_random/random_variants_PLINK_res_', i)
    makeOutDir1 <- "mkdir $outputdir"
    runPLINK <- paste0("plink --file ", outDir, "random_variants_", i, 
                       " --linear --pheno ", phenoFile, " --all-pheno ",
                       "--allow-no-sex --missing-phenotype 10121992 ", 
                       " --allow-extra-chr --out ", 
                       "$outputdir'/'random_variants_PLINK_res_", i)
    
    # write all commands to the script
    write(extrRandVars, scriptName, append = T)
    write(trPlink1, scriptName, append = T)
    write(trPlink2, scriptName, append = T)
    write(makeOutDir, scriptName, append = T)
    write(makeOutDir1, scriptName, append = T)
    write(runPLINK, scriptName, append = T)
  }
  message('Finished!')
}

# Plotting themes and pallets -------------------------------------------------
# ggplot2 theme
mashaGgplot2Theme <- list(
  theme_classic(base_size = 18) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "grey", size = 0.5, 
                                          linetype = 2))
)

# DO NOT CHANGE NAMES OR COLORS OF THIS PALLETE
# THEY ARE VERY IMPORTANT FOR CIRCOS PLOT
# THEY ARE HARD CODED IN /home/litovche/bin/circos-0.69-3/etc/colors.conf
rainbowPallete <- c('#C02F1E', '#D84E1C', '#F16C21', # reds
                    '#EF8B2D', '#ECAA39', '#ECC83E', # orange/yellow
                    '#CFF09E', '#A8DBA8', '#79BD9A', '#3B8686',  # green
                    '#1395B9', '#117899', '#0F5C78', #blues
                    '#d3d3d3', '#bdbdbd', '#a8a8a8', '#939393')  # greys
names(rainbowPallete) <- c('Frameshift', 'INDEL in tRNA', 'INDEL in lrRNA',
                           'INDEL in srRNA', 'Missense', 'intergenic INDEL',
                           'Synonymous', 'SNP in tRNA', 'SNP in lrRNA', 'SNP in srRNA',
                           "Inframe INS", "Inframe DEL", "intergenic SNP", 
                           "MNP")
pallete <- colorRampPalette(c('#C02F1E', '#EF8B2D', '#A3B86D', '#1395B9', 
                              '#0E3C52'))
# color pallete for hierarkical cluster of phenotypes
phenoColors <- colorRampPalette(c("#BB4444", "#EE9988", "#77AADD", "#4477AA"))