# FILE: 31_GRD_phenotypes -----------------------------------------------------
# USAGE: 
#
# DESCRIPTION: performs accociation of GRDs with the phenotypes, tests 
#              mito + nucl + mito:nucl model vs mito + nucl
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  11.10.2018
# REVISION: 11.10.2018
source('2_functions.R')

# Functions -------------------------------------------------------------------
#' GRDwithOnePheno
#' Checks on influence of GRDs on 1(!) phenotype via full and reduced models
#' @param mitoGeno mitochoondrial genotype
#' @param nuclGeno nuclear genotype
#' @param onePheno values for 1 phenotype
#' @param cutoffOnLines cut off on the number of lines with valid genotypes and
#'                      phenotypes
#' @param minNumbOfIndPerGeno minimal number of individuals per mitoNucl
#'                            genotype
#' @return NULL, if data doesn't pass the cutoff and list of fits (full model,
#' without interaction, anova between them) if it does
GRDwithOnePheno <- function(mitoGeno, nuclGeno, onePheno, cutoffOnLines = 150,
                            minNumbOfIndPerGeno = 5) {
  # get all the data into 1 dataframe
  mitoNuclPhenDF <- cbind(mitoGeno, nuclGeno, onePheno)
  # remove missing data and heteroplasms
  mitoNuclPhenDF <- mitoNuclPhenDF[apply(mitoNuclPhenDF, 1, 
                                         function(x) all(x[1:2] != "0" & x[1:2] != "2")), ]
  mitoNuclPhenDF <- mitoNuclPhenDF[complete.cases(mitoNuclPhenDF), ]
  
  result <- list(full = NA, notInter = NA, modelComp = data.frame("Pr(>F)" = NA))
  colnames(result$modelComp) <- "Pr(>F)"
  
  if (!is.null(mitoNuclPhenDF) & 
      (is.data.frame(mitoNuclPhenDF) | is.matrix(mitoNuclPhenDF))) {
    # if we have enough of cases
    if (nrow(mitoNuclPhenDF) >= cutoffOnLines) {
      mitoNuclPhenDF <- as.data.frame(mitoNuclPhenDF)
      # need to check, that the minimal number of lines per mito-nucl combination
      # with the phenotype is > 5
      contTab <- table(as.data.frame(mitoNuclPhenDF[, 1:2])) 
      contTab <- sort(contTab)
      if ((min(contTab) >= minNumbOfIndPerGeno) | 
          (min(contTab) == 0 & contTab[2] > minNumbOfIndPerGeno)) {
        mitoNuclPhenDF <- as.data.frame(mitoNuclPhenDF)
        mitoNuclPhenDF[, 3] <- as.numeric(as.character(mitoNuclPhenDF[, 3]))
        colnames(mitoNuclPhenDF) <- c('M', 'N', 'P')
        fitFull <- lm(P ~ M*N, data = mitoNuclPhenDF)
        fitNoIntr <- lm(P ~ M + N, data = mitoNuclPhenDF)
        modelComp <- anova(fitFull, fitNoIntr)
        result <- list(full = fitFull, notInter = fitNoIntr, 
                       modelComp = modelComp)
      }
    }
  }
  
  result
}

#' GRDswithAllPheno
#' Applies GRDwithOnePheno to all phenotypes: checks on influence of GRDs on
#' phenotypes via full and reduced models
#' @param phenoMatr matrix of the phenotypes, phenotypes in rows
#' @param mitoGOI genotype of mito variant
#' @param nuclGOI genotype of nuclear variant
#' @return list of results from GRDwithOnePheno
GRDswithAllPheno <- function(phenoMatr, mitoGOI, nuclGOI, ...) {
  apply(phenoMatr, 1,  function(y) GRDwithOnePheno(mitoGOI, nuclGOI, y, ...))
}

#' summarisePvals
#' Summarizes pvalues from results of GRDswithAllPheno
#' @param grdsPhenoRes results of GRDswithAllPheno
#' @return data frame with M (mitovar ID), N (nuclvar ID), P (phenotype ID), 
#' pvalue of anova comparison between models
summarisePvals <- function(grdsPhenoRes) {
  result <- lapply(1:length(grdsPhenoRes), 
                   function(i) do.call(rbind, 
                                       lapply(1:length(grdsPhenoRes[[i]]),
                                              function(j) data.frame(M = strsplit(names(grdsPhenoRes)[i], 
                                                                                  ' ')[[1]][1],
                                                                     N = strsplit(names(grdsPhenoRes)[i], 
                                                                                  ' ')[[1]][2],
                                                                     P = names(grdsPhenoRes[[i]])[j],
                                                                     pvalAOV = grdsPhenoRes[[i]][[j]]$modelComp$`Pr(>F)`[2]))))
  result <- do.call(rbind, result)
  result
}

#' plotMitoNuclPheno
#' Plots gwas-like box plot for the interaction between mito and nucl and 
#' phenotype
#' @param mitoGeno mitochoondrial genotype
#' @param nuclGeno nuclear genotype
#' @param onePheno values for 1 phenotype
#' @return ggplot2 
plotMitoNuclPheno <- function(mitoGOI, nuclGOI, onePheno) {
  dtToPlot <- data.frame(M = mitoGOI, nuclGOI, P = onePheno)
  dtToPlot <- dtToPlot[complete.cases(dtToPlot), ]
  dtToPlot <- dtToPlot[!abu$M %in% c('0', '2') & 
                       !dtToPlot$N %in% c('0', '2'), ]
  ggplot(dtToPlot, aes(x = factor(N), y = P, colour = M)) + 
    geom_boxplot() + geom_point(position = position_jitterdodge()) + 
    mashaGgplot2Theme
}

#' grdToGenoVec
#' Converts mitochondrial and nuclear genotype of a GRD into vector of merged 
#' genotype
#' @param mitoGeno mitochoondrial genotype
#' @param nuclGeno nuclear genotype
#' @return named vector of merged genotype, i.e. 33, 11, 13, 31. Heteroplams
#' and missing data are excluded
grdToGenoVec <- function(mitoGeno, nuclGeno) {
  genoTab <- data.frame(M = mitoGeno, N = nuclGeno, stringsAsFactors = F)
  # remove missing data and heteroplasms
  genoTab <- genoTab[apply(genoTab, 1, function(x) all(x != "0" & x != "2")), ]
  genoTab <- genoTab[complete.cases(genoTab), ]
  genoVec <- paste0(genoTab$M, genoTab$N)
  names(genoVec) <- rownames(genoTab)
  genoVec
}

#' hyperGeomOneGRD
#' Performs hypergeometric test of enrichment of mito-nucl genotype in tops
#' or bottoms of phenotypic distribution
#' @param onePheno values for 1 phenotype
#' @param grdGenoVect named vector of merged genotype, i.e. 33, 11, 13, 31,
#'                    result of grdToGenoVec
#' @param extrNumb number of DGRP lines to take from the top/bottom of of the
#'                 phenotypic distribution
#' @param minNumbCases minimal number of cases (GRD+Pheno) in order for test to
#'                     be performed
#' @return data.table with columns Extrem, Geno and hgPval
hyperGeomOneGRD <- function(onePheno, grdGenoVect, extrNumb = 10, 
                            minNumbCases = 100) {
  # select not NA phenotype values
  onePheno <- onePheno[!is.na(onePheno)]
  # select the once which have GRDs as well
  onePheno <- onePheno[intersect(names(grdGenoVect), names(onePheno))]
  # restrict GRDs to the same individuals
  grdGenoVect <- grdGenoVect[names(onePheno)]
  
  # init results
  result <- data.table(Extrem = rep(c('top', 'bottom'), 4),
                       Geno = rep(c('11', '13', '31', '33'), each = 2),
                       hgPval = rep(NA, 8))
  
  # check, that there are at least 100 cases
  if (length(grdGenoVect) >= minNumbCases) {
    # sort phenotype values
    onePheno <- sort(onePheno)
    
    # get dgrps at bottom and top 
    #bottom10 <- names(onePheno[1:round(perc * length(onePheno))])
    #top10 <- names(onePheno[round((1 - perc) * length(onePheno)):
    #                        length(onePheno)])
    bottom10 <- names(onePheno[1:extrNumb])
    top10 <- names(onePheno[(length(onePheno) - extrNumb):length(onePheno)])
    
    bottomHG <- sapply(unique(grdGenoVect), hyperGeomOneMNgeno, grdGenoVect,
                       bottom10, onePheno)
    topHG <- sapply(unique(grdGenoVect),  hyperGeomOneMNgeno, grdGenoVect,
                    top10, onePheno)
    
    # all possible genotypes in results - I do that to get uniformity of results
    result <- data.table(Extrem = rep(c('top', 'bottom'), 4),
                         Geno = rep(c('11', '13', '31', '33'), each = 2),
                         hgPval = c(topHG['11'], bottomHG['11'],
                                    topHG['13'], bottomHG['13'],
                                    topHG['31'], bottomHG['31'],
                                    topHG['33'], bottomHG['33']))
  }
  result                     
}

#' hyperGeomOneMNgeno
#' Calculates hypergeometric test for the enrichment of one mito-nuclear geno (i.e. 33) 
#' @param mnGeno
#' @param vectOfGeno result of grdToGenoVec: named vector of merged genotype,
#' i.e. 33, 11, 13, 31. Heteroplams and missing data are excluded
#' @param sampledIndvFromPheno
#' @param aviablIndvWithPheno
#' @return 
hyperGeomOneMNgeno <- function(mnGeno, vectOfGeno, sampledIndvFromPheno,
                               aviablIndvWithPheno) {
  # 1st : number of dgrp lines with the mnGeno which are found in 
  # sampledIndvFromPheno
  # 2nd: total number of DGRP lines posessing MN geno
  # 3rd: total number of DGRP lines NOT posessing MN geno
  # 4th: number of dgrp lines derived from phenotype distribution : it's 
  # basically sampledIndvFromPheno
  phyper(sum(names(vectOfGeno[vectOfGeno == mnGeno]) %in% 
             sampledIndvFromPheno) - 1,
         sum(vectOfGeno == mnGeno), 
         length(aviablIndvWithPheno) - sum(vectOfGeno == mnGeno), 
         length(sampledIndvFromPheno), lower.tail = F)
}

# Inputs ----------------------------------------------------------------------
# read in GRDs
signGRDs <- readRDS('results/GRD_DGRP2_Indep12/GRDs_DGRP2_Indep12_FDR01.Rds')
signGRDsAnno <- signGRDs[[2]]
signGRDs <- signGRDs$GRDs[passNeigh == 1]

# read-in mito variants
beversDm6Path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
mitoVcfdm6 <- readInMitoVCF(beversDm6Path, 'dm6', 14917, 
                            removeRefLines = T, interGen = "merge")
mitoVcfGeno <- mitoVcfdm6[[2]]
# leave only variants tested for GRDs
mitoClust <- readRDS('results/mitoClust.Rds')
mitoVcfGeno <- mitoVcfGeno[mitoClust[tested == 1]$MitoVar, ]
mitoVcfGeno <- mitoVcfGeno[unique(signGRDs$MitoVar), ]
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/0", 1, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/1", 3, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/0", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/1", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("\\.", 0, x))
mitoVcfGeno <- as.matrix(mitoVcfGeno)

# read-in nuclear variants in GRDs
nuclVCFGeno <- readInNuclVCF('results/GRD_DGRP2_Indep12/GRDs_DGRP2_Indep12_FDR01_NuclVars.annot.vcf',
                             'dm6')
nuclVCFGeno <- nuclVCFGeno[[2]]
rownames(nuclVCFGeno) <- gsub('_SNP|_DEL|_INS|_MNP', '', rownames(nuclVCFGeno))
nuclVCFGeno <- nuclVCFGeno[unique(signGRDs$NuclVar), ]
nuclVCFGeno <- apply(nuclVCFGeno, 2, function(x) gsub("0/0", 1, x))
nuclVCFGeno <- apply(nuclVCFGeno, 2, function(x) gsub("1/1", 3, x))
nuclVCFGeno <- apply(nuclVCFGeno, 2, function(x) gsub("1/0", 2, x))
nuclVCFGeno <- apply(nuclVCFGeno, 2, function(x) gsub("0/1", 2, x))
nuclVCFGeno <- apply(nuclVCFGeno, 2, function(x) gsub("\\.", 0, x))
nuclVCFGeno <- as.matrix(nuclVCFGeno)
addZero <- grepl('DGRP-..$', colnames(nuclVCFGeno))
colnames(nuclVCFGeno)[addZero] <- gsub('DGRP-', 'DGRP-0',
                                       colnames(nuclVCFGeno)[addZero])

# read-in phenotypes
phenosPath <- 'results/phenotype.collection.v2.2311.indep.cont.adj.phe'
phenos <- read.table(phenosPath, header = T, row.names = 1)
phenos[phenos == 10121992] <- NA
phenos <- phenos[, -1]
phenos <- t(phenos)

# read-in fitness phenotypes
load('results/obj_FitnessPhenotypeData.out')
fitnessPhenos <- colnames(FitnessPhenotypeData)
FitnessPhenosReduced <- c("Phe00005","Phe00010","Phe00040","Phe00041","Phe00042","Phe00043",
                          "Phe00149","Phe00150","Phe00151","Phe00152","Phe00153","Phe00154","Phe00155","Phe00156",
                          "Phe00157","Phe00158","Phe00159","Phe00160","Phe00199","Phe00201","Phe00202","Phe00203",
                          "Phe00204","Phe00205","Phe00206","Phe00207","Phe00208","Phe00209","Phe00210","Phe00211",
                          "Phe00212","Phe00213","Phe00214","Phe00215","Phe00219","Phe00220","Phe00223","Phe00226",
                          "Phe00227","Phe00228","Phe00229","Phe00230","Phe00231","Phe00232","Phe00233","Phe00234",
                          "Phe00235","Phe00236","Phe00237","Phe00238","Phe00239","Phe00240","Phe00287",
                          "Phe00331","Phe00332","Phe00333","Phe00334",
                          "Phe00382","Phe00384","Phe00386",
                          "Phe00693","Phe00695","Phe00696",
                          "Phe00731","Phe00734","Phe00737","Phe00740","Phe00743","Phe00746",
                          "Phe00801","Phe00802","Phe00803","Phe00804","Phe00806","Phe00807","Phe00808")
# read-in extremeness phenotypes
load('results/obj_FitnessExtremes.out')
rownames(FitnessExtremes) <- gsub('_', '-', rownames(FitnessExtremes))
addZero <- grepl('DGRP-..$', rownames(FitnessExtremes))
rownames(FitnessExtremes)[addZero] <- gsub('DGRP-', 'DGRP-0', 
                                           rownames(FitnessExtremes)[addZero])
addNA <- setdiff(colnames(phenos), rownames(FitnessExtremes))
for (i in 1:length(addNA)) {
  FitnessExtremes <- rbind(FitnessExtremes, rep(NA, ncol(FitnessExtremes)))
  rownames(FitnessExtremes)[nrow(FitnessExtremes)] <- addNA[i]
}
FitnessExtremes <- as.data.frame(FitnessExtremes)
# add them to main phenos
phenos <- rbind(phenos,
                FitnessExtremes[colnames(phenos), ]$MediocreFraction,
                FitnessExtremes[colnames(phenos), ]$ExtremeFraction)
rownames(phenos)[(nrow(phenos) - 1) : nrow(phenos)] <- c('MediocreFraction',
                                                         'ExtremeFraction')

# information about phenotypes
load('results/phenotypes/obj_phenotype.information.v2.0902.out')
phenotype.information.v2 <- as.data.table(phenotype.information.v2)
phenotype.information.v2[, PhenoID := gsub('PID', 'Phe', PhenoID)]
setkey(phenotype.information.v2, PhenoID)
# cluster of phenotypes
phenoClusters <- readRDS('results/phenoClusters.Rds')
setkey(phenoClusters, phenotype)
phenoClustPack <- phenoClusters[,.(paste0(phenotype,collapse = ",")),
                                by = cluster]

# make sure that column ordering is the same
nuclVCFGeno <- nuclVCFGeno[, colnames(mitoVcfGeno)]
phenos <- phenos[, colnames(mitoVcfGeno)]

# Number of GRDs having >5 DGRPs per mito-nucl genotypes ----------------------
dgrpPerGeno <- apply(signGRDs, 1, 
                     function(x) {
                       genoTab <- data.frame(mitoVcfGeno[x['MitoVar'], ],
                                             nuclVCFGeno[x['NuclVar'], ],
                                             stringsAsFactors = F)
                       # remove missing data and heteroplasms
                       genoTab <- genoTab[apply(genoTab, 1, 
                                                function(x) all(x[1:2] != "0" & 
                                                                x[1:2] != "2")), ]
                       genoTab <- genoTab[complete.cases(genoTab), ]
                       table(genoTab)})
# check distribution
summary(apply(dgrpPerGeno, 2, min))
# avaible for testing interactions
avaiblInter <- signGRDs[which(apply(dgrpPerGeno, 2, min) > 5)]
# check, that Sxl hacve some GRDs avaible for the interaction testing
'Sxl' %in% signGRDsAnno[NuclVar %in% avaiblInter$NuclVar]$Gene

# Check mito-nuclear interaction on significance for each phenotype -----------
# Actually, it doesn't make sence to run comparison of full model, like:
# phenotype ~ mito + nucl + mito:nucl to the reduced model, like:
# phenotype ~ mito + nucl by anova, because there's not enough of data to back 
# up interactions. Only makes sence to run it on avaiblInter

# check, how many lines in total we have for each avaiblInter GRD
min(colSums(dgrpPerGeno[, which(apply(dgrpPerGeno, 2, min) > 5)]))
# I put cutoffOnLines = 100 to get more GRDs to test
# RUN ON VITAL IT
#grdsPhen3models <- lapply(1:nrow(avaiblInter), 
#                          function(x) GRDswithAllPheno(phenos, 
#                                                       mitoVcfGeno[avaiblInter$MitoVar[x], ],
#                                                       nuclVCFGeno[avaiblInter$NuclVar[x], ], 
#                                                       cutoffOnLines = 100,
#                                                       minNumbOfIndPerGeno = 5))
#names(grdsPhen3models) <- paste(avaiblInter$MitoVar, avaiblInter$NuclVar)

# summarize into data table
grdPhenSumm <- c()                                  
for (i in 1:nrow(avaiblInter)) {
  print(i)
  grdsPhen3models <- readRDS(paste0('results/GRDs_vs_Pheno/grdsPhen3models_',
                                    i, '.Rds'))
  grdPhenSumm <- rbind(grdPhenSumm, summarisePvals(grdsPhen3models))
}
grdPhenSumm <- grdPhenSumm[order(grdPhenSumm$pvalAOV), ]
grdPhenSumm <- as.data.table(grdPhenSumm)
setnames(grdPhenSumm, 'P', 'PhenoID')
setnames(grdPhenSumm, 'M', 'MitoVar')
setnames(grdPhenSumm, 'N', 'NuclVar')
grdPhenSummCC <- grdPhenSumm[complete.cases(grdPhenSumm), ]
grdPhenSummCC[, pvalAOVadj := p.adjust(pvalAOV, method = 'fdr')]
grdPhenSumm <- rbind(grdPhenSummCC, 
                     grdPhenSumm[!complete.cases(grdPhenSumm), ], fill = T)
grdPhenSumm[, ID := paste(MitoVar, NuclVar)]
# add gene
setkey(signGRDsAnno, NuclVar)
geneInfo <- signGRDsAnno[as.character(unique(unlist(grdPhenSumm$NuclVar))), ]
geneInfo <- geneInfo[, c('NuclVar', 'Gene', 'inGeneLoc', 'EffectSize')]
setkey(grdPhenSumm, NuclVar)
setkey(geneInfo, NuclVar)
grdPhenSumm <- merge(grdPhenSumm, geneInfo, all = T, allow.cartesian = T)
# add phenotype association
grdPhenSumm[, PhenoID := as.character(PhenoID)]
fullPhen <- phenotype.information.v2[grdPhenSumm$PhenoID]$FullPhenotype
grdPhenSumm[, FullPhenotype := fullPhen]
grdPhenSumm <- grdPhenSumm[order(pvalAOV), ]
saveRDS(grdPhenSumm, 'results/grdPhenSumm.Rds')

# try Benjamini 2013 correction
# remove annotation so there will be no repetative lines
grdPhenSumm <- grdPhenSumm[, 1:6, with = T]
grdPhenSumm <- grdPhenSumm[!duplicated(grdPhenSumm)]
grdPhenSumm <- grdPhenSumm[order(PhenoID, ID)]
grdPhenPvalMatr <- getPvalMatr(grdPhenSumm, pvalName = "pvalAOV", 
                               varName = 'ID')
# remove lines and columns with all NA
grdPhenPvalMatr <- grdPhenPvalMatr[, apply(grdPhenPvalMatr, 2, 
                                           function(x) !all(is.na(x)))]
grdPhenPvalMatr <- grdPhenPvalMatr[apply(grdPhenPvalMatr, 1, 
                                         function(x) !all(is.na(x))), ]
grdPhenSign <- pAdjMultPhenBenjamini2013(grdPhenPvalMatr, 'BH', 0.3)
grdPhenSign <- as.data.table(grdPhenSign)
grdPhenSign[, MitoVar := sapply(as.character(Vars),
                                function(x) strsplit(x, ' ')[[1]][1])]
grdPhenSign[, NuclVar := sapply(as.character(Vars),
                                function(x) strsplit(x, ' ')[[1]][2])]
# add gene anno
geneInfo <- signGRDsAnno[as.character(unique(unlist(grdPhenSign$NuclVar))), ]
geneInfo <- geneInfo[, c('NuclVar', 'Gene', 'inGeneLoc', 'EffectSize')]
setkey(grdPhenSign, NuclVar)
setkey(geneInfo, NuclVar)
grdPhenSign <- merge(grdPhenSign, geneInfo, all = T, allow.cartesian = T)
# add phenotype association
grdPhenSign[, Pheno := as.character(Pheno)]
fullPhen <- phenotype.information.v2[grdPhenSign$Pheno]$FullPhenotype
grdPhenSign[, FullPhenotype := fullPhen]
saveRDS(grdPhenSign, 'results/grdPhenSumm_Benj2013.Rds')

# fitness phenotypes only
grdFitPhenSumm <- grdPhenSumm[PhenoID %in% fitnessPhenos, ]
grdFitPhenSumm <- grdPhenSumm[PhenoID %in% FitnessPhenosReduced, ]
grdFitPhenSummCC <- grdFitPhenSumm[complete.cases(grdFitPhenSumm), ]
grdFitPhenSummCC[, pvalAOVadj := p.adjust(pvalAOV, method = 'fdr')]
grdFitPhenSumm <- rbind(grdFitPhenSummCC, 
                        grdFitPhenSumm[!complete.cases(grdFitPhenSumm), ], fill = T)
# add gene
setkey(signGRDsAnno, NuclVar)
geneInfo <- signGRDsAnno[as.character(unique(unlist(grdFitPhenSumm$NuclVar))), ]
geneInfo <- geneInfo[, c('NuclVar', 'Gene', 'inGeneLoc', 'EffectSize')]
setkey(grdFitPhenSumm, NuclVar)
setkey(geneInfo, NuclVar)
grdFitPhenSumm <- merge(grdFitPhenSumm, geneInfo, all = T, allow.cartesian = T)
# add phenotype association
grdFitPhenSumm[, PhenoID := as.character(PhenoID)]
fullPhen <- phenotype.information.v2[grdFitPhenSumm$PhenoID]$FullPhenotype
grdFitPhenSumm[, FullPhenotype := fullPhen]
grdFitPhenSumm <- grdFitPhenSumm[order(pvalAOV), ]
saveRDS(grdPhenSumm, 'results/grdFitPhenSumm.Rds')
# gives 0 hits
grdPhenSumm <- grdPhenSumm[order(PhenoID, ID)]
grdFitPvalMatr <- getPvalMatr(grdFitPhenSumm, pvalName = "pvalAOV", 
                              varName = 'ID')
grdFitPvalMatr <- grdFitPvalMatr[, apply(grdFitPvalMatr, 2, 
                                         function(x) !all(is.na(x)))]
grdFitPvalMatr <- grdFitPvalMatr[apply(grdFitPvalMatr, 1, 
                                       function(x) !all(is.na(x))), ]
grdFitPhenSign <- pAdjMultPhenBenjamini2013(grdFitPvalMatr, 'BH', 0.3)

plotMitoNuclPheno(mitoVcfGeno['chrM:3892_A/C', ], 
                  nuclVCFGeno['chr3R_10477352', ], phenos['Phe00164', ]) 

# Hypergeometric test for the enrichment in extremes of phenos ----------------
# for every GRD compute genotype per line
GRDsGenoPerLine <- apply(signGRDs, 1,
                         function(x) grdToGenoVec(mitoVcfGeno[x['MitoVar'], ],
                                                  nuclVCFGeno[x['NuclVar'], ]))
names(GRDsGenoPerLine) <- paste(signGRDs$MitoVar, signGRDs$NuclVar,
                                sep = '_vs_')
hgGRDsPhenos <- c()
for (i in 1:length(GRDsGenoPerLine)) {
  print(i)
  hgOneGRDonePhen <- apply(phenos, 1, hyperGeomOneGRD, GRDsGenoPerLine[[i]])
  hgOneGRDonePhen <- lapply(1:length(hgOneGRDonePhen), 
                            function(x) hgOneGRDonePhen[[x]][, Pheno := rownames(phenos)[x]])
  hgOneGRDonePhen <- lapply(1:length(hgOneGRDonePhen), 
                            function(x) hgOneGRDonePhen[[x]][, MitoVar := signGRDs$MitoVar[i]])
  hgOneGRDonePhen <- lapply(1:length(hgOneGRDonePhen), 
                            function(x) hgOneGRDonePhen[[x]][, NuclVar := signGRDs$NuclVar[i]])
  hgOneGRDonePhen <- do.call(rbind, hgOneGRDonePhen)
  hgGRDsPhenos <- rbind(hgGRDsPhenos, hgOneGRDonePhen)
}
saveRDS(hgGRDsPhenos, paste0('results/hgGRDsPhenos.Rds'))

hgGRDsPhenos <- readRDS(paste0('results/hgGRDsPhenos_top20perc.Rds'))
hgGRDsPhenos[, ID := paste(MitoVar, NuclVar, sep = '-')]
setnames(hgGRDsPhenos, 'Pheno', 'PhenoID')

# check p-value distribution
summary(hgGRDsPhenos$hgPval)
summary(p.adjust(hgGRDsPhenos$hgPval, method = 'fdr'))
hgGRDsPhenosPassMTC <- list()
for (mnGeno in unique(hgGRDsPhenos$Geno)) {
  for (extr in unique(hgGRDsPhenos$Extrem)) {
    message(paste(mnGeno, extr))
    hgMnGenoExtrPvalMatr <- getPvalMatr(hgGRDsPhenos[Geno == mnGeno & 
                                                     Extrem == extr], 
                                        pvalName = 'hgPval', varName = 'ID')
    hgMnGenoExtrBenj2013 <- pAdjMultPhenBenjamini2013(hgMnGenoExtrPvalMatr,
                                                      'BH', 0.2)
    if (!is.null(hgMnGenoExtrBenj2013)) {
      hgGRDsPhenosPassMTC[[length(hgGRDsPhenosPassMTC) + 1]] <- hgMnGenoExtrBenj2013
    } else {
      hgGRDsPhenosPassMTC[[length(hgGRDsPhenosPassMTC) + 1]] <- NA
    }
    names(hgGRDsPhenosPassMTC)[length(hgGRDsPhenosPassMTC)] <- paste(mnGeno, 
                                                                     extr,
                                                                     sep = '_')
  }
}

# only extremeness phenotype
summary(hgGRDsPhenos[PhenoID %in% c("MediocreFraction", 
                                  "ExtremeFraction")]$hgPval)
summary(p.adjust(hgGRDsPhenos[PhenoID %in% c("MediocreFraction", 
                                           "ExtremeFraction")]$hgPval,
                 method = 'fdr')) 
summary(hgGRDsPhenos$hgPval)
summary(p.adjust(hgGRDsPhenos$hgPval, method = 'fdr'))
hgGRDsExtremPassMTC <- list()
for (mnGeno in unique(hgGRDsPhenos$Geno)) {
  for (extr in unique(hgGRDsPhenos$Extrem)) {
    message(paste(mnGeno, extr))
    hgMnGenoExtrPvalMatr <- getPvalMatr(hgGRDsPhenos[Geno == mnGeno & 
                                                     Extrem == extr &
                                                     PhenoID %in% c("MediocreFraction", 
                                                                    "ExtremeFraction")], 
                                        pvalName = 'hgPval', varName = 'ID')
    hgMnGenoExtrBenj2013 <- pAdjMultPhenBenjamini2013(hgMnGenoExtrPvalMatr,
                                                      'BH', 0.2)
    if (!is.null(hgMnGenoExtrBenj2013)) {
      hgGRDsExtremPassMTC[[length(hgGRDsExtremPassMTC) + 1]] <- hgMnGenoExtrBenj2013
    } else {
      hgGRDsExtremPassMTC[[length(hgGRDsExtremPassMTC) + 1]] <- NA
    }
    names(hgGRDsExtremPassMTC)[length(hgGRDsExtremPassMTC)] <- paste(mnGeno, 
                                                                     extr,
                                                                     sep = '_')
  }
}

# only fitness phenotypes
summary(hgGRDsPhenos[PhenoID %in% fitnessPhenos]$hgPval)
summary(p.adjust(hgGRDsPhenos[PhenoID %in% fitnessPhenos]$hgPval,
                 method = 'fdr')) 
hgGRDsFitPassMTC <- list()
for (mnGeno in unique(hgGRDsPhenos$Geno)) {
  for (extr in unique(hgGRDsPhenos$Extrem)) {
    message(paste(mnGeno, extr))
    hgMnGenoExtrPvalMatr <- getPvalMatr(hgGRDsPhenos[Geno == mnGeno & 
                                                       Extrem == extr &
                                                       PhenoID %in% fitnessPhenos], 
                                        pvalName = 'hgPval', varName = 'ID')
    hgMnGenoExtrBenj2013 <- pAdjMultPhenBenjamini2013(hgMnGenoExtrPvalMatr,
                                                      'BH', 0.2)
    if (!is.null(hgMnGenoExtrBenj2013)) {
      hgGRDsFitPassMTC[[length(hgGRDsFitPassMTC) + 1]] <- hgMnGenoExtrBenj2013
    } else {
      hgGRDsFitPassMTC[[length(hgGRDsFitPassMTC) + 1]] <- NA
    }
    names(hgGRDsFitPassMTC)[length(hgGRDsFitPassMTC)] <- paste(mnGeno, 
                                                                     extr,
                                                                     sep = '_')
  }
}

# only reduced fitness phenotypes
summary(hgGRDsPhenos[PhenoID %in% FitnessPhenosReduced]$hgPval)
summary(p.adjust(hgGRDsPhenos[PhenoID %in% FitnessPhenosReduced]$hgPval,
                 method = 'fdr')) 
hgGRDsFitRedPassMTC <- list()
for (mnGeno in unique(hgGRDsPhenos$Geno)) {
  for (extr in unique(hgGRDsPhenos$Extrem)) {
    message(paste(mnGeno, extr))
    hgMnGenoExtrPvalMatr <- getPvalMatr(hgGRDsPhenos[Geno == mnGeno & 
                                                       Extrem == extr &
                                                       PhenoID %in% FitnessPhenosReduced], 
                                        pvalName = 'hgPval', varName = 'ID')
    hgMnGenoExtrBenj2013 <- pAdjMultPhenBenjamini2013(hgMnGenoExtrPvalMatr,
                                                      'BH', 0.2)
    if (!is.null(hgMnGenoExtrBenj2013)) {
      hgGRDsFitRedPassMTC[[length(hgGRDsFitRedPassMTC) + 1]] <- hgMnGenoExtrBenj2013
    } else {
      hgGRDsFitRedPassMTC[[length(hgGRDsFitRedPassMTC) + 1]] <- NA
    }
    names(hgGRDsFitRedPassMTC)[length(hgGRDsFitRedPassMTC)] <- paste(mnGeno, 
                                                                     extr,
                                                                     sep = '_')
  }
}

# One little bar plot ---------------------------------------------------------
MitoVarToPlot <- 'chrM:7424_G/A'
NuclVarToPlot <- 'chrX_14792895'
PhenoToPlot <- 'Phe00331'
all(colnames(mitoVcfGeno) == colnames(nuclVCFGeno))
toPlot <- data.frame(M = mitoVcfGeno[MitoVarToPlot, ], 
                     N = nuclVCFGeno[NuclVarToPlot, ])
toPlot <- toPlot[apply(toPlot, 1, 
                function(x) sum(!sapply(x, function(y) y %in% c(0, 2))) == 2), ]
toPlot$M <- gsub(1, 'mREF', toPlot$M)
toPlot$M <- gsub(3, 'mALT', toPlot$M)
toPlot$N <- gsub(1, 'nREF', toPlot$N)
toPlot$N <- gsub(3, 'nALT', toPlot$N)
toPlot$P <- phenos[PhenoToPlot, rownames(toPlot)]
toPlot <- toPlot[complete.cases(toPlot), ]
toPlot$DGRP <- rownames(toPlot)
toPlot <- toPlot[order(toPlot$P), ]
toPlot$DGRP <- factor(toPlot$DGRP, levels = toPlot$DGRP)
ggplot(toPlot, aes(x = DGRP, y = P, fill = interaction(M, N))) + 
  geom_bar(stat = "identity") + mashaGgplot2Theme +
  ggtitle(paste('Pheno00331 - hypergeometric test on enrichment', 
                'of a GRD chrM:7424_G/A - chrX_14792895 in the', 
                'top 20% lines', sep = '\n')) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(fill=guide_legend(title="Genotype"))

# Trash -----------------------------------------------------------------------

abu <- function(x) {
  genoTab <- data.frame(mitoVcfGeno[x['MitoVar'], ],
                        nuclVCFGeno[x['NuclVar'], ],
                        stringsAsFactors = F)
  # remove missing data and heteroplasms
  genoTab <- genoTab[apply(genoTab, 1, 
                           function(x) all(x[1:2] != "0" & 
                                             x[1:2] != "2")), ]
  genoTab <- genoTab[complete.cases(genoTab), ]
  genoVec <- paste0(genoTab[, 1], genoTab[, 2])
  names(genoVec) <- rownames(genoTab)
  # select phenotype, sort, ensure that we have the same indv as in GRDs
  pheno <- sort(phenos[1, intersect(names(genoVec), colnames(phenos))])
  genoVec <- genoVec[names(pheno)]
  phenoBottom10 <- names(pheno[1:round(0.1 * length(pheno))])
  phenoTop10 <- names(pheno[round(0.9 * length(pheno)):length(pheno)])
  
  mnGeno <- '33'
  phyper(sum(names(genoVec[genoVec == mnGeno]) %in% phenoTop10) - 1,
         sum(genoVec == mnGeno), 
         length(pheno) - sum(genoVec == mnGeno), 
         length(phenoTop10), lower.tail = F)
}
mitoGeno <- mitoVcfGeno['chrM:14666_T/TTTA', ]
nuclGeno <- nuclVCFGeno['chr3L_14145740', ]
onePheno <- phenos['Phe00162',]
cutoffOnLines <- 150
minNumbOfIndPerGeno <- 5

set.seed(12)
f1<-gl(n=2,k=30,labels=c("Low","High"))
f2<-as.factor(rep(c("A","B","C"),times=20))
modmat<-model.matrix(~f1*f2,data.frame(f1=f1,f2=f2))
coeff<-c(1,3,-2,-4,1,-1.2)
y<-rnorm(n=60,mean=modmat%*%coeff,sd=0.1)
dat<-data.frame(y=y,f1=f1,f2=f2)

ggplot(dat, aes(x = factor(f1), y = y, colour = f2)) + 
  geom_boxplot() + 
  geom_point(data = dat, aes(y = y)) +
  geom_line(data = dat, aes(y = y, group = f2)) + 
  theme_bw()