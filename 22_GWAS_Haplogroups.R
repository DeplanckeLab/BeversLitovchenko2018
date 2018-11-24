# FILE: 17_GWAS_Haplogroups ---------------------------------------------------
# USAGE: 
#
# DESCRIPTION: performs GWAS on haplogroups
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  22.03.2018
# REVISION: 22.03.2018
# https://brownmath.com/stat/anova1.html
# http://www.blackwellpublishing.com/specialarticles/jcn_8_304.pdf

setwd('~/Desktop/MitoSeq_RPJB_ML_Aug2017/')
setwd('~/Desktop/BitBucket/MitoSeq_RPJB_ML_Aug2017/')
source('2_functions.R')
library(ggpubr)
library(mgcv)
library(multcomp)

# FUNCTIONS -------------------------------------------------------------------
#' prepDataForModel
#' Prepares phenotype - haplogroup pairs for modelling: removes NA, removes 
#' haplogroups which have less than 5 entries per phenotype, translates 
#' haplogroups to factor
#' @param phenoVect vector of phenotypes
#' @param haplGRvect vector of haplogroups
#' @return data.table with columns Phenotype and Haplogroup
prepDataForModel <- function(phenoVect, haplGRvect) {
  result <- cbind(phenoVect, haplGRvect)
  # remove NAs
  result <- result[complete.cases(result), ]
  colnames(result) <- c('Phenotype', 'Haplogroup')
  
  # remove clusters for which number of memebers is < 5
  sampPerHaplGR <- table(result$Haplogroup)
  sampPerHaplGR <- sampPerHaplGR[sampPerHaplGR >= 5]
  result <- result[Haplogroup %in% names(sampPerHaplGR)]
  
  # turn into factors 
  result[, Haplogroup := factor(Haplogroup, levels = sort(unique(Haplogroup)))]
  
  result
}

# INPUTS ----------------------------------------------------------------------
# phenotype collection -----
phenosPath <- 'results/phenotype.collection.v2.2311.indep.cont.adj.phe'
phenos <- fread(phenosPath)
phenos[phenos == 10121992] <- NA
phenosDT <- phenos[, -1, with = F]
names(phenosDT)[1] <- 'DGRP'
setkey(phenosDT, 'DGRP')

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

# Clusters of DGRPs -----
dgrpClust <- fread('results/haplotype.v5.txt', header = T, 
                   stringsAsFactors = F)
dgrpClust <- dgrpClust[!Line %in% c('iso-1_rep1', 'Ore-R_rep1', 'Ore-R-rep2',
                                    'w1118_rep1', 'w1118_rep2')]
dgrpClust[, DGRP := gsub('.*_', '', Line)]
addZero <- grepl('^..$', dgrpClust$DGRP)
dgrpClust$DGRP[addZero] <- paste0('0', dgrpClust$DGRP[addZero])
dgrpClust$DGRP <- paste0('DGRP-', dgrpClust$DGRP)
dgrpClust <- dgrpClust[!DGRP %in% c('DGRP-338', 'DGRP-356')]

# merge phenotypes and haplogroups in one data frame -----
setkey(dgrpClust, DGRP)
phenosAndHaplGR <- merge(dgrpClust, phenosDT)
phenosAndHaplGR[, Clusters_IBS := NULL]
phenosAndHaplGR[, DGRP := NULL]
phenosAndHaplGR[, Line := NULL]
colnames(phenosAndHaplGR)[1] <- 'Haplogroup'

# path to the variants from Garlapow -----
garlapowVars <- 'results/Garlapow_FoodIntake_SignVars.csv'
garlapowVars <- fread(garlapowVars, header = T, sep = '\t')
setkey(garlapowVars, 'Variant ID')

# path to the epistasis for food intake ----
#FIepistMale <- fread('results/PLINK/PLINK_nuclAll_epistasis/PLINK_nuclAll_epistasis.Phe00162.epi.noNAN.qt')
#FIepistMale[, P.ADJ := p.adjust(P, method = 'BH')] 
#FIepistFemale <- fread('results/PLINK/PLINK_nuclAll_epistasis/PLINK_nuclAll_epistasis.Phe00166.epi.noNAN.qt')

# colors ----
dgrpClustColors <- c("#DDCC77", '#B11623', "#88CCEE", "#44AA99", "#008C9E", 
                     "#999933", '#0B486B', "#EFBDC0", "#AA4499", "#8C315D",
                     "#998FB8", "#483078", "#117733") #"#661100")
names(dgrpClustColors) <- sort(unique(dgrpClust$Haplotypes_TCS))

# MULTIPLE COMPARISON PROCEDURE: AOV + TUKEY ----------------------------------
haploGrTukeyGWAS <- c()
haploGrANOVAGWAS <- c()
# takes quite a long time to compute
for (phenoName in colnames(phenosAndHaplGR)[-1]) {
  print(phenoName)
  forModel <- prepDataForModel(phenosAndHaplGR[, phenoName, with = F],
                               phenosAndHaplGR[, 'Haplogroup', with = F])
  # step 1: anova
  anovaRes <- aov(Phenotype ~ Haplogroup, data = forModel)
  haploGrANOVAGWAS[[length(haploGrANOVAGWAS) + 1]] <- summary(anovaRes)
  # step 2: tukey test
  tukeySumm <- summary(glht(anovaRes, linfct = mcp(Haplogroup = "Tukey")))
  haploGrTukeyGWAS[[length(haploGrTukeyGWAS) + 1]] <- tukeySumm
}
names(haploGrANOVAGWAS) <- colnames(phenosAndHaplGR)[-1]
saveRDS(haploGrANOVAGWAS, 'results/haploGrANOVAGWAS.Rds')
names(haploGrTukeyGWAS) <- colnames(phenosAndHaplGR)[-1]
saveRDS(haploGrTukeyGWAS, 'results/haploGrTukeyGWAS.Rds')
haploGrTukeyGWAS <- readRDS('results/haploGrTukeyGWAS.Rds')
haploGrANOVAGWAS <- readRDS('results/haploGrANOVAGWAS.Rds')

# get significant associations
haploGrANOVAPvals <- sapply(haploGrANOVAGWAS, function(x) x[[1]]$`Pr(>F)`[1])
# techically, here we need to restrict to phenotypes passing 0.05 cutoff on 
# anova, but food intake doesn't pass =(
haploGrTukeyPvals <- lapply(haploGrTukeyGWAS,
                            function(x) as.character(x$test$pvalues))
names(haploGrTukeyPvals) <- colnames(phenosAndHaplGR)[-1]
haploGrTukeyPvalsSign <- lapply(haploGrTukeyPvals, function(x) any(x < 0.05))
haploGrTukeyPvalsMin <- as.numeric(sapply(haploGrTukeyPvals, min))
names(haploGrTukeyPvalsMin) <- names(haploGrTukeyPvals)
# get significant phenotypes
haploGrANOVAPvals <- haploGrANOVAPvals[names(haploGrTukeyPvalsMin)]
haploGrSign <- names(which(haploGrANOVAPvals < 0.3 & 
                           haploGrTukeyPvalsMin < 0.05))
haploGrTukeyGWASsign <- haploGrTukeyGWAS[haploGrSign]
# print info about them
message(paste('Significantly associated phenotypes:',
              paste(phenotype.information.v2[haploGrSign]$FullPhenotype,
                    collapse = ', ')))

# make a table for the paper(?)
# calculate the max number of between-groups comparison
numbOfGroups <- sapply(haploGrTukeyGWASsign, 
                       function(x) length(x$test$coefficients))
# get all the possible comparisons between groups names
groupsComps <- haploGrTukeyGWASsign[[which(numbOfGroups == 
                                     max(numbOfGroups))[1]]]
groupsComps <- names(groupsComps$test$tstat)
haploGrTukeyGWASsignPvals <- lapply(haploGrTukeyGWASsign, 
                                    function(x) setNames(x$test$pvalues,
                                                         names(x$test$tstat)))
haploGrSignSumTab <- lapply(haploGrTukeyGWASsignPvals, 
                            function(x) x[groupsComps])
haploGrSignSumTab <- do.call(rbind, haploGrSignSumTab)
colnames(haploGrSignSumTab) <- groupsComps
haploGrSignSumTab <- as.data.frame(haploGrSignSumTab)
haploGrSignSumTab <- cbind(as.character(phenotype.information.v2[rownames(haploGrSignSumTab)]$FullPhenotype),
                           haploGrSignSumTab)
write.table(haploGrSignSumTab, 'results/haploGrSignSumTab.csv', col.names = T,
            row.names = T, sep = '\t')
# restrict to the columns, which have significance:
#haploGrSignSumTab <- haploGrSignSumTab[, apply(haploGrSignSumTab, 2,
#                                               function(x) any(na.omit(as.numeric(x)) < 0.05))]

# plot bar plots 
pdf('plots/haploGrTukeySign.pdf', width = 14, height = 8)
for (phenoName in rownames(haploGrSignSumTab)) {
  print(phenoName)
  forModel <- prepDataForModel(phenosAndHaplGR[, phenoName, with = F],
                               phenosAndHaplGR[, 'Haplogroup', with = F])
  # for the significance bars
  signClust <- which(as.numeric(haploGrSignSumTab[phenoName, ]) < 0.05)
  signClust <- colnames(haploGrSignSumTab)[signClust]
  fromClust <- sapply(signClust, function(x) strsplit(x, ' - ')[[1]][1])
  toClust <- sapply(signClust, function(x) strsplit(x, ' - ')[[1]][2])
  signLevel <- haploGrSignSumTab[phenoName, -1][haploGrSignSumTab[phenoName,
                                                                  -1] < 0.05]
  signLevel <- round(na.omit(signLevel), 2)
  
  print(ggplot(forModel, aes(x = Haplogroup, y = Phenotype, 
                             fill = Haplogroup)) +
          geom_boxplot() + geom_jitter(position = position_jitter(0.1)) + 
          mashaGgplot2Theme + 
          ylab(phenotype.information.v2[phenoName]$FullPhenotype) +
          ggtitle(paste('GWAS association between', 
                        phenotype.information.v2[phenoName]$FullPhenotype, 
                        'and haplogroups - Tukey test')) +
          guides(fill = F) + 
          scale_fill_manual(values = dgrpClustColors[levels(forModel$Haplogroup)]) +
          geom_signif(y_position = seq(1.1, 1.3,
                                       length.out = length(fromClust)) *
                                                    max(forModel$Phenotype),
                      xmin = fromClust, xmax = toClust, 
                      annotation = as.vector(signLevel), tip_length = 0.05))
}
dev.off()
