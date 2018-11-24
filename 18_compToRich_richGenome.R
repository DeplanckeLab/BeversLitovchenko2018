# FILE: 4_compToRich_richGenome -----------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: compares Bevers et al data set aligned to chrU, chr3L, etc as in
# the Richardson paper methods and compares it to Richardson variants
#
#  OPTIONS:  none
#  REQUIREMENTS:  none
#  BUGS: --
#  NOTES:  ---
#  AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
#  COMPANY:  EPFL, Lausanne, Switzerland
#  VERSION:  1
#  CREATED:  26.07.2017
#  REVISION: 14.06.2018

# FUNCTIONS -------------------------------------------------------------------
#' readRichardson
#' Reads richardson file
#' @param pathToTab path to the table from richardson paper
#' @return data.frame
readRichardson <- function(pathToTab) {
  richardson <- read.table(pathToTab, header  = T, stringsAsFactors = F)
  richardson <- richardson[-2:-1, ]
  richardson <- rbind(colnames(richardson), richardson)
  richardson <- t(richardson)
  colnames(richardson) <- richardson[1, ]
  colnames(richardson)[1] <- 'pos'
  richardson <- richardson[-1, ]
  richardson <- as.data.frame(richardson, stringsAsFactor = F)
  richardson$pos <- as.integer(gsub('X', '', richardson$pos))
  richardson
}

# INPUTS ----------------------------------------------------------------------
baseDir <- '~/Desktop/MitoSeq_RPJB_ML_Aug2017/'
setwd(baseDir)
source('2_functions.R')

richardsonDir <- 'comparison_with_richardson/'
# coverage information for Bevers data, mapped on richardson genome
deplCovPath <- paste0(richardsonDir, 'align_to_rich_stats.txt')
# table to convert to right GT
gtConvTabPath <- 'results/genotyping/genotyping_DGRP2_dm3.csv'
# Variants called from Bevers data mapped on richardson genome(chrU, wolbachia, 
# part of 3L)
deplVcfForComparPath <- paste0(richardsonDir, 
                               'MitoSeq_AllRuns_richardson_chrU_rigthGT.fltr.vcf')

# coverage information for Richardson data, downloaded from Supplement
richCovPath <- paste0(richardsonDir, 'Richardson_Dataset_S1.txt')
# Variants from Richardson paper, just donloaded from Supplement
richNotInfPath <- paste0(richardsonDir, 
                         'DGRP_DPGP_mtDNA_annotated_noindels.variants.txt')
richInfPath <- paste0(richardsonDir, 
                      'DGRP_DPGP_mtDNA_infected_annotated_noindels.variants.txt')

pallete <- colorRampPalette(c('#C02F1E', '#EF8B2D', '#A3B86D', '#1395B9', 
                              '#0E3C52'))

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

# OUTPUTS ---------------------------------------------------------------------
plotOutDir <- 'plots'

# READ-IN ---------------------------------------------------------------------
# Read-in and preprocess Richardson coverage
richCov <- fread(richCovPath)
richCov <- richCov[project == 'DGRP']
colnames(richCov)[3] <- 'Name'
richCov[, Name := gsub('DGRP', 'DGRP-', gsub(".*-", '', Name))]
# select only informative columns
richCov <- richCov[, c('Name', 'num_reads', 'num_bp', 'mtdna_num_reads_mapped',
                       'mtdna_depth_of_coverage')]
colnames(richCov) <- c('Name', 'NR_total', 'num_bp', 'NR_U', 'chrMmedianCov')
addZero <- grepl('DGRP-..$', richCov$Name)
richCov[addZero]$Name <- gsub('DGRP-', 'DGRP-0', richCov[addZero]$Name)
richCov[, Study := 'Richardson']
richCov <- richCov[order(chrMmedianCov)]
richCov[, Rank := 1:nrow(richCov)]

# Read-in and preprocess Bevers coverage
deplCov <- fread(deplCovPath, sep = ' ', fill = T)
setkey(deplCov, Name)
# assign right GT and select good samples
gtConvTab <- fread(gtConvTabPath)
# exclude lines which will not be used anyway and w- and ore
gtConvTab <- gtConvTab[will_file_be_used == 'Y' & DGRP_1stMatch != 'NA']
gtConvTab <- gtConvTab[, c('DEPL_line', 'DGRP_1stMatch')]
gtConvTab[, DEPL_line := gsub("_dm3", '', DEPL_line)]
gtConvTab[, DGRP_1stMatch := gsub('line_', 'DGRP-', DGRP_1stMatch)]
addZero <- grepl('DGRP-..$', gtConvTab$DGRP_1stMatch)
gtConvTab[addZero]$DGRP_1stMatch <- gsub('DGRP-', 'DGRP-0', 
                                         gtConvTab[addZero]$DGRP_1stMatch)
colnames(gtConvTab)[1] <- 'Name'
setkey(gtConvTab, Name)
# add real gt to coverage table
deplCov <- merge(gtConvTab, deplCov)
deplCov <- deplCov[, -c('Name')]
colnames(deplCov)[1] <- 'Name'
deplCov[, Study := 'Bevers']
deplCov <- deplCov[order(chrMmedianCov)]
deplCov[, Rank := 1:nrow(deplCov)]

# restrict tables of coverage to the common lines only
setkey(richCov, Name)
setkey(deplCov, Name)
richCov <- richCov[intersect(richCov$Name, deplCov$Name)]
deplCov <- deplCov[richCov$Name]

# COMPARE COVERAGE ------------------------------------------------------------
# in general
bothCov <- rbind(deplCov, richCov, fill = T)
png(paste0(plotOutDir, '/chrUcoverage_BeversVsRichardson.png'),
    height = 800, width = 1600)
ggplot(bothCov, aes(x = Rank, y = chrMmedianCov)) +
  geom_point(shape = 20, colour = pallete(5)[1], size = 3) +
  xlab('DGRP line') + ylab('chrU coverage') + 
  ggtitle(paste0('Comparisoon of coverage of chrU between Bevers et al and ',
                 "Richardson et al\n",
                 "Mapping to chrU, chr3L:10Mb–11.2Mb, Wolbachia")) + 
  theme_classic(base_size = 20) + mashaGgplot2Theme + facet_grid(. ~ Study)
dev.off()

# per DGRP 
setkey(richCov, Name)
setkey(deplCov, Name)
mergedCov <- merge(deplCov, richCov)
colnames(mergedCov) <- gsub('.x', 'Bevers', colnames(mergedCov))
colnames(mergedCov) <- gsub('.y', 'Rich', colnames(mergedCov))
png(paste0(plotOutDir, '/chrUcoverageAcrossDGRP_BeversVsRichardson.png'),
    height = 800, width = 800)
ggplot(mergedCov, aes(x = chrMmedianCovBevers, y = chrMmedianCovRich)) +
  geom_point(shape = 20, size = 3, colour = pallete(5)[1]) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  xlab('chrU coverage, Bevers et al') +
  ylab('chrU coverage, Richardson et al') + 
  ggtitle(paste0('Comparison of coverage of chrU between Bevers et al\nand ',
                 "Richardson et al\n",
                 "Mapping to chrU, chr3L:10Mb–11.2Mb, Wolbachia")) +
  mashaGgplot2Theme
dev.off()

# REFORMAT VCF TO RICHARDSON FORMAT -------------------------------------------
deplVCF <- readInMitoVCF(deplVcfForComparPath, 'dm3', codingEnds = 14917, 
                         removeRefLines = T, interGen = "noChange")
deplVCFGR <- deplVCF[[1]]
# select only SNPs and SNP-like MNPs. MNPs because in Richardosn set they are
# present
deplVCFGR <- deplVCFGR[!getVariantStructType(deplVCFGR) %in% c('DEL', 'INS')]
indelMNP <- c('chrU:5960_AATATAT/A', 'chrU:5967_TTATA/TTA', 'chrU:8185_A/AT',
              'chrU:13136_TA/TAA')
deplVCFGR <- deplVCFGR[!names(deplVCFGR) %in% indelMNP]
# genotype info
deplVCFGeno <- deplVCF[[2]][names(deplVCFGR), ]

deplNewFormat <- c()
for (i in 1:nrow(deplVCFGeno)) {
  ref <- tolower(as.character(mcols(deplVCFGR[i])$REF))
  alt <- tolower(as.character(mcols(deplVCFGR[i])$ALT[[1]]))
  deplNewFormat <- rbind(deplNewFormat, ifelse(deplVCFGeno[i, ] == '.' |
                                                 deplVCFGeno[i, ] == '0/0', ref,
                                               ifelse(deplVCFGeno[i, ] == '1/1',
                                                      alt, 'HET')))
}
deplNewFormat <- cbind(start(deplVCFGR), deplNewFormat)

colnames(deplNewFormat)[1] <- 'pos'
colnames(deplNewFormat) <- gsub('DGRP0', 'DGRP', colnames(deplNewFormat))

deplNewFormat <- as.data.frame(deplNewFormat,stringsAsFoctors = F)
deplNewFormat$pos <- as.integer(as.character(deplNewFormat$pos))

# COMPARE VARIANTS BETWEEN STUDIES --------------------------------------------
richardsonNotInf <- readRichardson(richNotInfPath)
richardsonInf <- readRichardson(richInfPath)
setdiff(richardsonInf$pos, richardsonNotInf$pos)
setdiff(richardsonNotInf$pos, richardsonInf$pos)
# not infected has more positions
setdiff(colnames(richardsonNotInf), colnames(richardsonInf))
# and more DGRPlines
setdiff(colnames(richardsonInf), colnames(richardsonNotInf))
# so, basically, we can use only non-infected

# restrict to common lines only
richardsonNotInf <- richardsonNotInf[, intersect(colnames(deplNewFormat),
                                                 colnames(richardsonNotInf))]
deplNewFormat <- deplNewFormat[, colnames(richardsonNotInf)]

# since we restricted to common lines only, it might be that the position is
# not polymorphic anymore
richardsonNotInf <- richardsonNotInf[apply(richardsonNotInf[, -1], 1, 
                                           function(x) length(unique(x))) != 1, ]
deplNewFormat <- deplNewFormat[apply(deplNewFormat[, -1], 1, 
                                     function(x) length(unique(x))) != 1, ]

# sort according to possition, just for convinience
richardsonNotInf <- richardsonNotInf[order(richardsonNotInf$pos), ]
deplNewFormat <- deplNewFormat[order(deplNewFormat$pos), ]

# restrict data sets only to the coding part
richardsonNotInf <- richardsonNotInf[richardsonNotInf$pos < 14917, ]
deplNewFormat <- deplNewFormat[deplNewFormat$pos < 14917, ]

forVenn <- list(Richardson = richardsonNotInf$pos,
                Deplancke = deplNewFormat$pos)
plotInVenn <- Venn(forVenn)

png(paste0(plotOutDir, '/chrU_UniqSites_BeversVsRichardson.png'), 
    height = 800, width = 800)
plot(plotInVenn, doWeights = T)
dev.off()

print('Unique to Richardson set')
richSpecPos <- print(plotInVenn@IntersectionSets$`10`)
print('Unique to Deplancke set')
deplSpecPos <- print(plotInVenn@IntersectionSets$`01`)

# taking into account the allele  - not run
# deplAllAlleles <- c()
# for (i in 1:nrow(deplNewFormat)) {
#  for (j in 2:ncol(deplNewFormat)) {
#    deplAllAlleles <- c(deplAllAlleles, paste(deplNewFormat[i, 1],
#                                              deplNewFormat[i, j],
#                                              colnames(deplNewFormat)[j]))
#  }
#}

#richAllAlleles <- c()
#for (i in 1:nrow(richardsonNotInf)) {
#  for (j in 2:ncol(richardsonNotInf)) {
#    richAllAlleles <- c(richAllAlleles, paste(richardsonNotInf[i, 1],
#                                              richardsonNotInf[i, j],
#                                              colnames(richardsonNotInf)[j]))
#  }
#}

#forVenn <- list(Richardson = richAllAlleles, Deplancke = deplAllAlleles)
#plotInVenn <- Venn(forVenn)

#png(paste0(plotOutDir, '/chrU_SitesAlleles_BeversVsRichardson.png'), 
#    height = 800, width = 800)
#plot(plotInVenn, doWeights = T)
#dev.off()

# PLOT NUMBER OF STUDY SPECIFIC VARIANTS VS OBSERVED FREQUENCY ----------------
richSpecGeno <- richardsonNotInf[richardsonNotInf$pos %in% richSpecPos, ]
richSpecMAC <- apply(richSpecGeno[, -1], 1, function(x) sort(table(x))[1])
richSpecMAC <- as.data.frame(table(richSpecMAC))
richSpecMAC$Study <- 'Richardson'
colnames(richSpecMAC) <- c('Freq', 'VarCount', 'Study')

deplSpecGeno <- deplNewFormat[deplNewFormat$pos %in% deplSpecPos, ]
deplSpecMAC <- apply(deplSpecGeno[, -1], 1, function(x) sort(table(x))[1])
deplSpecMAC <- as.data.frame(table(deplSpecMAC))
deplSpecMAC$Study <- 'Bevers'
colnames(deplSpecMAC) <- c('Freq', 'VarCount', 'Study')

studySpecMAC <- rbind(richSpecMAC, deplSpecMAC)
colnames(studySpecMAC) <- c('Freq', 'VarCount', 'Study')
studySpecMAC$Freq <- as.integer(as.character(studySpecMAC$Freq))

png(paste0(plotOutDir, '/chrU_StudySpec_Freq.png'), 
    height = 800, width = 800)
ggplot(studySpecMAC, aes(x = Freq, y = VarCount, fill = Study)) + 
  geom_bar(stat = "identity", col = 'black', position = 'dodge') +
  geom_text(aes(studySpecMAC$Freq, studySpecMAC$VarCount, 
                label = studySpecMAC$VarCount), vjust = -0.5, 
            position = position_dodge(width = 1)) +
  mashaGgplot2Theme + ylab('Number of variants') + xlab('Frequency observed') +
  ggtitle(paste('Frequency of study-specific variants in DGRPs')) +
  scale_x_continuous(breaks = seq(1, 10, 1)) + ylim(c(0, 60)) +
  scale_fill_manual(values = c("#999999", "#E69F00"))
dev.off()

# OUTPUT TABLES WITH STUDY-SPECIFIC VARIANTS ----------------------------------
# check on max number of alleles
max(apply(richSpecGeno[, -1], 1, function(x) length(table(x))))
richSpecTab <- data.frame(Position = richSpecGeno$pos,
                          A1 = apply(richSpecGeno[, -1], 1, 
                                     function(x) names(sort(-table(x)))[1]),
                          A2 = apply(richSpecGeno[, -1], 1, 
                                     function(x) names(sort(-table(x)))[2]))
richSpecTab$DGRP <- sapply(1:nrow(richSpecTab), 
                           function(x) paste(colnames(richSpecGeno)[richSpecGeno[x, ] == 
                           as.character(richSpecTab[x, 3])], collapse = ', '))
write.table(richSpecTab, 'richardsonSpec_chrU.csv', row.names = F, quote = F,
            sep = '\t')

# check on max number of alleles
max(apply(deplSpecGeno[, -1], 1, function(x) length(table(x))))
deplSpecTab <- data.frame(Position = deplSpecGeno$pos,
                          A1 = apply(deplSpecGeno[, -1], 1, 
                                     function(x) names(sort(-table(x)))[1]),
                          A2 = apply(deplSpecGeno[, -1], 1, 
                                     function(x) names(sort(-table(x)))[2]),
                          A3 = apply(deplSpecGeno[, -1], 1, 
                                     function(x) names(sort(-table(x)))[3]))
deplSpecTab$DGRP_A2 <- sapply(1:nrow(deplSpecTab), 
                              function(x) paste(colnames(deplSpecGeno)[deplSpecGeno[x, ] == 
                              as.character(deplSpecTab[x, 3])], collapse = ', '))
deplSpecTab$DGRP_A3 <- sapply(1:nrow(deplSpecTab), 
                              function(x) paste(na.omit(colnames(deplSpecGeno)[deplSpecGeno[x, ] == 
                              as.character(deplSpecTab[x, 4])]), collapse = ', '))
write.table(deplSpecTab, 'beversSpec_chrU.csv', row.names = F, quote = F,
            sep = '\t')

# Number of study-specific variants per DGRP line -----------------------------
plotVarsPerDGRP(richSpecGeno, "Richardson")

# number of DGRPs with deplancke-specific variants
deplcSpecDGRPs <- apply(deplNewFormat[deplNewFormat$pos %in% 
                                      plotInVenn@IntersectionSets$`01`, -1],
                        1, function(x) table(x)[c('a', 'c', 'g', 't', 'HET')])
deplcSpecDGRPs <- rbind(deplNewFormat[deplNewFormat$pos %in% 
                                      plotInVenn@IntersectionSets$`01`, ]$pos,
                        deplcSpecDGRPs)
rownames(deplcSpecDGRPs) <- c('pos', 'a', 'c', 'g', 't', 'HET')
deplcSpecDGRPs <- as.data.table(t(deplcSpecDGRPs))
table(apply(deplcSpecDGRPs[is.na(HET), ], 1, function(x) min(na.omit(x))))