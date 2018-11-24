# FILE: 11_parse_PhenoDB_forGWAS ----------------------------------------------
# USAGE: 
#
# DESCRIPTION: converts phenoDB from Roel to plink-suitable format
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  08.02.2018
# REVISION: 08.02.2018
setwd('~/Desktop/MitoSeq_RPJB_ML_Aug2017/')
setwd('~/Desktop/BitBucket/MitoSeq_RPJB_ML_Aug2017/')
source('2_functions.R')
library(corrplot)
library(factoextra)
library(NbClust)
library(klaR)

# Functions -------------------------------------------------------------------
#' getUniqueStudies
#' Returns unique studies corresponding to genotypes
#' @param phenosInOneCluster phenotype codes in one cluster
#' @param phenoInfoTab table, containing information about phenotypes, i.e.
#'        phenotype.information.v2
getUniqueStudies <- function(phenosInOneCluster, phenoInfoTab) {
  unique(phenoInfoTab[phenosInOneCluster$phenotype, ]$Study)
}

# Inputs ----------------------------------------------------------------------
# annotated vcf file, mapping on dm6
mitoVcfdm6path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
mitoVcfdm6 <- readInMitoVCF(mitoVcfdm6path, 'dm6', 14917, 
                            removeRefLines = T, interGen = "merge")
mitoVcfGR <- mitoVcfdm6[[1]]
mitoVcfGeno <- mitoVcfdm6[[2]]

sampsWithGenoCut <- 150 # cut off on the number of samples
MAFcutoff <- 0.05 # MAF cut off

# loading phenoDB
load('results/phenotypes/obj_phenotype.collection.v2.2311.out')
load('results/phenotypes/obj_phenotype.information.v2.0902.out')
phenotype.information.v2$PhenoID <- gsub('PID', 'Phe', 
                                         phenotype.information.v2$PhenoID)
# covariates file
wolbInvPath <- 'results/phenotypes/wolb_and_inv_DGRP2.txt'

# Outputs ---------------------------------------------------------------------
outputDir <- 'results/'

# Selecting independent mitochondrial variants --------------------------------
# remove MNPs
mitoVcfGeno <- mitoVcfGeno[getVariantStructType(mitoVcfdm6[[1]]) != 'MNP', ]
# I exclude heteroplasmic variants in individuals, I do not though the whole 
# variant away!
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/0", ".", x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/1", ".", x))
# select variants with number of genotyped samples > cut off
mitoVcfGeno <- mitoVcfGeno[apply(mitoVcfGeno, 1, 
                                 function(x) sum(x != './.' & 
                                                 x != '.') >= sampsWithGenoCut), ]
# select variants with MAF > MAFcutoff
mitoVcfGeno <- mitoVcfGeno[apply(mitoVcfGeno, 1,
                                 function(x) calcMacMaf(x)['MAF'] >= MAFcutoff),]

mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/0", 1, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/1", 3, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("\\.", 0, x))
mitoVcfGeno <- as.matrix(mitoVcfGeno)

# get clusters of mito variants
#mitoClust <- findVarClusters(mitoVcfGeno)
#mitoClust <- data.table(MitoVar = names(mitoClust), Cluster = mitoClust)
# this mito clusters are based on manual selection. 
mitoClust <- data.table(MitoVar = c("chrM:4247_C/T", "chrM:1154_C/T", 
                                    "chrM:12345_C/A",
                                    "chrM:12381_T/C", "chrM:13561_T/C", 
                                    "chrM:1512_C/T", "chrM:4616_T/A", 
                                    "chrM:6047_ATTAAT/A", "chrM:6308_C/A", 
                                    "chrM:8876_A/G", "chrM:9065_A/C", 
                                    "chrM:10226_C/T", "chrM:3583_T/C",
                                    "chrM:5396_C/T", "chrM:6989_G/A", 
                                    "chrM:7871_G/A", "chrM:8982_C/T", 
                                    "chrM:11128_C/T", "chrM:12682_T/TA",
                                    "chrM:12898_G/A", "chrM:2349_C/T",
                                    "chrM:5500_C/T", "chrM:11416_T/C",
                                    "chrM:2223_A/G", "chrM:3892_A/C",
                                    "chrM:4954_T/C", "chrM:12132_C/T",
                                    "chrM:14666_T/TTTA", "chrM:1716_G/A",
                                    "chrM:2071_C/T", "chrM:2661_C/T", 
                                    "chrM:6872_G/A", "chrM:7424_G/A",
                                    "chrM:791_G/A"),
                        Cluster = c("CL:1", "CL:1", "CL:1", "CL:1", "CL:1", "CL:1", 
                                    "CL:1", "CL:1", "CL:1", "CL:1", "CL:1", 
                                    "CL:1", "CL:1", "CL:1", "CL:1", "CL:1", 
                                    "CL:1", "CL:2", "CL:2", "CL:2", "CL:2", 
                                    "CL:2", "CL:3", "CL:3", "CL:4", "CL:4", 
                                    "chrM:12132_C/T", "chrM:14666_T/TTTA", 
                                    "chrM:1716_G/A", "chrM:2071_C/T", 
                                    "chrM:2661_C/T", "chrM:6872_G/A", 
                                    "chrM:7424_G/A", "chrM:791_G/A"),
                        Haplotypes_TCSspec = c(rep(NA, 17), rep('MH3', 5),
                                               rep('MH2', 2), rep('MH1', 2),
                                               NA, 'MH4', 'MH8b', 'MH11', NA,
                                               'MH6', 'MH7', 'MH8'),
                        Symbol = c(rep(18, 17), rep(17, 5), rep(16, 2), 
                                   rep(15, 2), 7:0))
setkey(mitoClust, 'MitoVar')

# Get independent variants: take one per mtDNA var cluster + all clusterless
independVars <- c(mitoClust[Cluster != "", ][,.(MitoVar[1]), by = Cluster]$V1,
                  mitoClust[Cluster == "", ]$MitoVar)
independVars <- sapply(independVars, function(x) strsplit(x, '_')[[1]][1])

# write bed file for vcf-tools
independVars <- data.frame(chr = sapply(independVars, 
                                        function(x) strsplit(x, ':')[[1]][1]),
                           start = sapply(independVars, 
                                          function(x) strsplit(x, ':')[[1]][2]),
                           end = sapply(names(independVars), 
                                        function(x) nchar(strsplit(x, '\\/')[[1]][2])))
independVars$start <- as.integer(as.character(independVars$start)) - 1
independVars$end <- independVars$start + 
                    as.integer(as.character(independVars$end))
write.table(independVars, 'results/independVars.bed', col.names = T, quote = F,
            row.names = F)
mitoClust[, tested := ifelse(MitoVar %in% rownames(independVars), 1, 0)]
saveRDS(mitoClust, 'results/mitoClust.Rds')

# create a separate VCF with independent variants only
system(paste('vcftools --vcf', mitoVcfdm6path, '--bed',
             'results/independVars.bed --remove-indv w1 --remove-indv w2',
             '--remove-indv w1118 --remove-indv Berk1',
             '--remove-indv Berk2 --remove-indv ore --remove-indv ore2',
             '--remove-indv ore2 --remove-indv ore3',
             '--recode --recode-INFO-all --out',
             'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.INDEPNDT'))
file.remove('results/independVars.bed')
system('sed -i "s@0/1@./.@g" results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.INDEPNDT.recode.vcf')
system('sed -i "s@1/0@./.@g" results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.INDEPNDT.recode.vcf')

# Parcing phenotype DB --------------------------------------------------------
phenotype.collection.v2 <- as.data.frame(phenotype.collection.v2)
# removing ral and stock numbers
phenotype.collection.v2 <- phenotype.collection.v2[, -3:-2]
# add zeros to two-digit DGRP names
phenotype.collection.v2$DGRP <- gsub('DGRP_', 'DGRP-',
                                     phenotype.collection.v2$DGRP)
addZero <- grepl('DGRP-..$', phenotype.collection.v2$DGRP)
phenotype.collection.v2$DGRP[addZero] <- gsub('DGRP-', 'DGRP-0', 
                                              phenotype.collection.v2$DGRP[addZero])
colnames(phenotype.collection.v2)[1] <- 'ID'
rownames(phenotype.collection.v2) <- phenotype.collection.v2$ID
# restricting to the lines present in our mitochondrial study
phenotype.collection.v2 <- phenotype.collection.v2[colnames(mitoVcfGeno), ]
# some phenotypes in the table are not real: they are SE or number of reads,etc
realPheno <- phenotype.information.v2[!phenotype.information.v2$PhenotypeClass %in% 
                                      c('Phenotype-like', 
                                        'Phenotype-Stat-Info'), ]$PhenoID
phenotype.collection.v2 <- phenotype.collection.v2[, c('ID', realPheno)]
# remove phenotypes with too much of NAs
less100NAs <- apply(phenotype.collection.v2, 2, function(x) sum(!is.na(x)))
phenotype.collection.v2 <- phenotype.collection.v2[, less100NAs >= 100]
# remove % sign
phenotype.collection.v2 <- apply(phenotype.collection.v2, 2, 
                                 function(x) gsub('%', '', x))
# put factorial phenotypes into separate table, because it's not possible to
# adjust them
continPheno <- phenotype.information.v2[phenotype.information.v2$Class == 
                                        'Continuous', ]$PhenoID
continPheno <- c('ID',
                 intersect(continPheno, colnames(phenotype.collection.v2)))
phenotype.collection.v2.cont <- phenotype.collection.v2[, continPheno]
phenotype.collection.v2.cont <- as.data.frame(phenotype.collection.v2.cont,
                                              stringsAsFactors =  F)
factPheno <- phenotype.information.v2[phenotype.information.v2$Class == 
                                     'Factorial', ]$PhenoID
factPheno <- c('ID', intersect(factPheno, colnames(phenotype.collection.v2)))
phenotype.collection.v2.fact <- phenotype.collection.v2[, factPheno]
phenotype.collection.v2.fact <- as.data.frame(phenotype.collection.v2.fact,
                                              stringsAsFactors =  F)

# Parcing covariats file ------------------------------------------------------
covarDF <- fread(wolbInvPath, header = T)
# remove bloomington 
covarDF <- covarDF[ , -2, with = F]
setnames(covarDF, 'line', 'ID')
covarDF[, Wolb := as.factor(Wolb)]
covarDF[, In_2L_t := as.factor(In_2L_t)]
covarDF[, In_2R_NS := as.factor(In_2R_NS)]
covarDF[, In_3R_P := as.factor(In_3R_P)]
covarDF[, In_3R_K := as.factor(In_3R_K)]
covarDF[, In_3R_Mo := as.factor(In_3R_Mo)]
# restrict to lines which we have
setkey(covarDF, ID)
covarDF <- covarDF[phenotype.collection.v2.cont$ID, ]

# Perform adjustment ----------------------------------------------------------
plinkPheno.cont <- PhenofromDF(phenotype.collection.v2.cont, 
                               adjusted = T, covarDF)
# set up missing phenotype value, rhe same should be indicated in plink command
plinkPheno.cont[is.na(plinkPheno.cont)] <- 10121992

write.table(plinkPheno.cont, paste0(outputDir, 
                                    'phenotype.collection.v2.2311.cont.adj.phe'),
            quote = F, row.names = F, sep = '\t')
write.table(phenotype.collection.v2.fact, 
            paste0(outputDir, 'phenotype.collection.v2.2311.fact.adj.phe'),
            quote = F, row.names = F, sep = '\t')

# Correlation of phenotypes ---------------------------------------------------
allPhenos <- plinkPheno.cont
rownames(allPhenos) <- allPhenos$Ind_ID
allPhenos <- allPhenos[, -2:-1]
allPhenos[allPhenos == 10121992] <- NA # I inserted my BD as NA value for PLINK
phenoCor <- cor(allPhenos, use = "pairwise.complete.obs")
phenoCorPval <- cor.mtest(phenoCor)$p
# Combine with significance
corrplot(phenoCor,  method = "color", col = phenoColors(200), type = "upper",
         order = "hclust", number.cex = .7, p.mat = phenoCorPval, 
         sig.level = 0.01, insig = "blank", diag = F, tl.pos = 'n', 
         title = paste('Correlation between', ncol(allPhenos), 'phenotypes'))

# Hierarchical clustering of the phenotypes -----------------------------------
# Compute hierarchical clustering and cut into 15 clusters
phenoCorCl <- corclust(allPhenos, mincor = 0.2)
# number of clusters if we cut at 0.8 correlation between phenotypes
phenoClusters <- cutree(phenoCorCl$clustering, h = 0.2)
phenoClusters <- data.table(phenotype = names(phenoClusters),
                            cluster = phenoClusters)

# plot
fviz_dend(phenoCorCl$clustering, cex = 0.5, xlab = 'Phenotypes', h = 0.2, 
          k_colors = phenoColors(length(unique(phenoClusters$cluster)) + 1),
          main = paste("Cluster Dendrogram for", ncol(allPhenos), 
                       'phenotypes divided in', 
                       length(unique(phenoClusters$cluster)), 'clusters'),
          ylab = '1 - minimum absolute correlation within cluster',
          horiz = T)
clustSize <- phenoClusters[,.N, by = cluster][order(-N), ]

# Select not correlated phenotypes --------------------------------------------
phenotype.information.v2 <- as.data.table(phenotype.information.v2)
setkey(phenotype.information.v2, PhenoID)

# clusters of size 1 - not actually a cluster
indepPheno <- phenoClusters[cluster %in% clustSize[N == 1, ]$cluster]$phenotype

# cluster of size 2: if phenos come from different studies, preserve both,
# if from one - take any
twoPhenInClust <- phenoClusters[cluster %in% clustSize[N == 2, ]$cluster]
twoPhenInClust[, cluster := as.character(cluster)]
setkey(twoPhenInClust, cluster)
twoPhenInClust <- lapply(unique(twoPhenInClust$cluster), 
                         function(x) ifelse(length(getUniqueStudies(twoPhenInClust[x], 
                                                   phenotype.information.v2)) == 1,
                                            twoPhenInClust[x][1],
                                            twoPhenInClust[x]))
intersect(indepPheno, unlist(twoPhenInClust)) # check
indepPheno <- c(indepPheno, unlist(twoPhenInClust))

# for the clusters of size 3 and more select the phenotype which is correlated
# with the rest the most
moreTwoPhenInClust <- c()
for (clust in unique(clustSize[N > 2, ]$cluster)) {
  corrInClust <- phenoCor[phenoClusters[cluster == clust]$phenotype, 
                          phenoClusters[cluster == clust]$phenotype]
  studiesInClust <- phenotype.information.v2[colnames(corrInClust), ]$Study
  numbStudiesInClust <- table(as.character(studiesInClust))
  for (study in names(numbStudiesInClust)) {
    if (numbStudiesInClust[study] == 1) {
      indepPheno <- c(indepPheno, phenotype.information.v2[phenotype.information.v2[PhenoID %in% rownames(corrInClust) &
                                                                                    Study == study]$PhenoIDStudy == study]$PhenoID)
    } else {
      phenosForStudy <- phenotype.information.v2[PhenoID %in% rownames(corrInClust) &
                                                 Study == study]$PhenoID
      corrInClustOnePheno <- corrInClust[phenosForStudy, phenosForStudy]
      corrInClustOnePheno[corrInClustOnePheno == 1] <- NA
      selectPheno <- which(corrInClustOnePheno == max(corrInClustOnePheno, 
                                                      na.rm = T), arr.ind = T)
      selectPheno <- rownames(selectPheno)[1]
      moreTwoPhenInClust <- c(moreTwoPhenInClust, unlist(selectPheno))
    }
  }
}
indepPheno <- c(indepPheno,moreTwoPhenInClust)

# write files for PLINK with independent phenotypes only - adjusted
write.table(plinkPheno.cont[, c('FID', 'IID', indepPheno)], 
            paste0(outputDir, 
                   'phenotype.collection.v2.2311.indep.cont.adj.phe'),
            quote = F, row.names = F, sep = '\t')
mitoClust[, tested := ifelse(MitoVar %in% rownames(independVars), 1, 0)]
saveRDS(phenoClusters, 'results/phenoClusters.Rds')

# write files for Fast LMM with independent phenotypes only - UNadjusted
# Put your phenotype in the Inputs/ folder and name it 
# *Pheno_Name*_Phenotype_Full.txt. It should follow the formatting 
# line_id \t phenotype value \t sex (m/f), and should not have a header.
phenotype.collection.v2.cont.FastLMM <- phenotype.collection.v2.cont[, c('ID',
                                                                  indepPheno)]
for (i in 2:ncol(phenotype.collection.v2.cont.FastLMM)) {
  # I remove zeros, because like that they will mach the ones in vcfs
  # create two identical with males and females to mimic Michael's input
  # structure
  toPrint_f <- data.frame(line.id = gsub('DGRP-0', 'DGRP-', 
                                         phenotype.collection.v2.cont.FastLMM$ID),
                          phenotype = phenotype.collection.v2.cont.FastLMM[, i], 
                          sex = 'f')
  toPrint_m <- data.frame(line.id = gsub('DGRP-0', 'DGRP-', 
                                    phenotype.collection.v2.cont.FastLMM$ID),
                          phenotype = phenotype.collection.v2.cont.FastLMM[, i], 
                          sex = 'm')
  toPrint <- rbind(toPrint_f, toPrint_m)
  toPrint$line.id <- gsub('DGRP-', 'line_', toPrint$line.id)
  outName <- paste0('results/Fast_LMM/Inputs/', 
                    colnames(phenotype.collection.v2.cont.FastLMM)[i],
                    '_Phenotype_Full.txt')
  write.table(toPrint, outName, col.names = T, row.names = F, sep = '\t',
              quote = F)
  message(paste0("pheno=", colnames(phenotype.collection.v2.cont.FastLMM)[i],
                 " sex=Female reproducibility=FALSE ", 
                 "use_official_gsm=TRUE ",
                 "./Bash_GWAS_Pipeline_Full.sh"))
}