# FILE: 3_baseAnalysis --------------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: plots overview of chrM coverage
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

# PATHS TO INPUTS -------------------------------------------------------------
baseDir <- '~/Desktop/MitoSeq_RPJB_ML_Aug2017/'
setwd(baseDir)
source('2_functions.R')

# there coding part ends
codingEnds <- 14917

# sequenceing stats, dm6
# file with sequencing statistics, such as coverage per chr, min, med and max
# coverage, create with script 7a_getMap_stats.sh
covPerChrPath <- 'results/coverage_stats/align_to_dm6_with_runs.csv'

# genotyping table
genotypeTabPath <- 'results/genotyping/genotyping_DGRP2_dm3.csv'

# result vcf, dm6
# annotated vcf file, mapping on dm6
mitoVcfdm6path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'

# dm6
# path to the genome, to which data were aligned
dm6Path <- '~/Documents/RefGen/dm6/dm6.Wolb.fa'
# annotation of the genome to which data were aligned
dm6MitoAnnotPath <- 'dm6_annotation/genes_CUSTOM_MANUAL_CHECK.gtf'

# repeats analysis, the table is from running java code MitoRepeats
repeatsPath <- 'results/intergenicRepeats/intergenicRepeats.csv'

# PATHS TO OUTPUTS ------------------------------------------------------------
svgWidth <- 8
pngWidth <- 1000
plotOutDir <- 'plots/'

# input files for circos plot
circosSNPs <- paste0(plotOutDir, 'circos_plot/mitoVarsLocation/SNPs_noRefLines.txt')
circosINDELs <- paste0(plotOutDir, 'circos_plot/mitoVarsLocation/INDELs_noRefLines.txt')
circosHIST <- paste0(plotOutDir, 'circos_plot/mitoVarsLocation/HIST_noRefLines.txt')
  
# READ IN VCF & ANNOTATION ----------------------------------------------------
mitoVcfdm6 <- readInMitoVCF(mitoVcfdm6path, 'dm6', codingEnds, 
                            removeRefLines = T, interGen = "merge")
mitoVcfGR <- mitoVcfdm6[[1]]
mitoVcfGeno <- mitoVcfdm6[[2]]
# annotation
mitoAnn <- readVariantAnnot(mitoVcfdm6path, mitoVcfGR)
mitoAnn <- simplifyAnnot(mitoAnn)
setkey(mitoAnn, Variant)
mitoAnn['chrM:838_AT/TA', c('StructType', 'FuncType', 
                            'EffectSize', 'Gene') := list('MNP', 'MNP', 
                                                          'MODERATE',
                                                          'mt:ND2')]
mitoAnn['chrM:5960_ATATATTTATATATATATATATAT/TTA', 
        c('StructType', 'FuncType', 
          'EffectSize', 'Gene') := list('MNP', 'MNP', 'MODIFIER', 
                                        'mt:ND3-mt:tRNA:A')]
# gene location
mitoGenes <- readChrMannot(dm6MitoAnnotPath)

# PLOT PIE CHART OF GENOTYPING ------------------------------------------------
# read table of genotyping
genotypeTab <- fread(genotypeTabPath, header = T, stringsAsFactors = F)
genotypeStats <- genotypeTab[, .N, Decision]
genotypeStats[is.na(Decision), Decision := 'w-, berk, ore']
setkey(genotypeStats, Decision)
genotypeStats <- genotypeStats[c('I DONT KNOW', 'w-, berk, ore', 
                                 'WEAK MATCH', 'MATCH', 'WEAK SWAP', 'SWAP')]

slices <- genotypeStats$N
lbls <- paste0(genotypeStats$Decision, ',\n', genotypeStats$N, ' samples')

X11 ()
pie(slices, labels = lbls, col = pallete(nrow(genotypeStats)), radius = 1,
    cex = 1.2, main = 'Results of genotyping (aligned on dm3, DGRP2 vcf)',
    cex.main = 2)
dev.print(png, filename = paste0(plotOutDir, '/1_genotypingStats.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/1_genotypingStats.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

message('Number of samples processed')
length(unique(genotypeTab[Decision %in% c('WEAK MATCH', 'MATCH', 
                                          'WEAK SWAP', 'SWAP')]$DGRP_1stMatch))

# PLOT COVERAGE RANGE PER DGRP ------------------------------------------------
# read min/max/median coverage file
covPerChr <- fread(covPerChrPath, header = T, stringsAsFactors = F)
covPerChr[RunType == "PE"]$RunType <- 'paired-end sequenced'
covPerChr[RunType == "SE"]$RunType <- 'single-end sequenced'
if (length(setdiff(covPerChr$Name, colnames(mitoVcfGeno)) != 0)) {
  message(paste("Removed",
                paste(setdiff(covPerChr$Name, colnames(mitoVcfGeno)), 
                      collapse = ', '), 'because they are not present in',
                'genotype table'))
  covPerChr <- covPerChr[Name %in% colnames(mitoVcfGeno)]
}
# if median coverage is < 10 - sequencing failed
failedSeq <- which(covPerChr$chrMmedianCov < 10)
if (length(failedSeq) > 0) {
  message("This samples failed: ")
  print(covPerChr[failedSeq]$DGRP)
  covPerChr <- covPerChr[chrMmedianCov > 10]
}
covPerChr <- covPerChr[order(covPerChr$chrMmedianCov), ]

X11()
ggplot(covPerChr, aes(x = 1:nrow(covPerChr), y = chrMmedianCov / 1000,
                      fill = RunType)) + 
  geom_bar(stat = "identity", color = "white") +
  scale_x_continuous(breaks = seq(1, nrow(covPerChr), by = 1), 
                     labels = NULL) +
  scale_y_continuous(breaks = round(seq(0, max(covPerChr$chrMmedianCov / 1000),
                                        by = 0.5), 1)) +
  xlab('Samples') + ylab('Coverage (x1000)') +
  scale_fill_manual(values = c('grey', 'steelblue')) + 
  mashaGgplot2Theme + theme(panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.75, 0.75), legend.title = element_blank())
#dev.print(png, paste0(plotOutDir, '/Supplemental_Figure_S3.png'), 
#          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/Supplemental_Figure_S3.svg'), 
          width = 12, height = 8)
dev.off()

X11()
ggplot(covPerChr[1:56, ], aes(x = 1:nrow(covPerChr[1:56, ]),
                              y = chrMmedianCov / 1000,
                      fill = RunType)) + 
  geom_bar(stat = "identity", color = "white") +
  scale_x_continuous(breaks = seq(1, nrow(covPerChr), by = 1), 
                     labels = NULL) +
  scale_y_continuous(breaks = round(seq(0, 0.4, by = 0.1), 1)) +
  scale_fill_manual(values = c('grey', 'steelblue')) +
  mashaGgplot2Theme + 
  theme(panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  guides(fill = F)
#dev.print(png, paste0(plotOutDir, '/Supplemental_Figure_S3.png'), 
#          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/Supplemental_Figure_S3_inlet.svg'), 
          width = 12, height = 8)
dev.off()

# PLOT MTDNA / AUTOSOMES COVERAGE ---------------------------------------------
covPerChr[, WolbStatus := ifelse(NR_Wolb / 1144540 >= 1, 'yes', 'no')]
covPerChr[, autoReadsNumb := sum(NR_2L, NR_2R, NR_3L, NR_3R, NR_4, NR_X), 
            by = Name]
covPerChr[, ratio := (80 * NR_mtDNA / 19514) / 
                     (80 * autoReadsNumb / 132532477)]
covPerChr <- covPerChr[order(covPerChr$ratio), ]

X11()
ggplot(covPerChr, aes(x = 1:nrow(covPerChr),  y = ratio)) + 
       geom_bar(colour = pallete(5)[5], stat = "identity", 
                fill = pallete(5)[5]) +
       xlab('DGRP line') + ylab('mtDNA / autosomes coverage ratio') +
       ggtitle('mtDNA:nDNA coverage') + mashaGgplot2Theme
dev.print(png, paste0(plotOutDir, '/3_mtDNAtoAutosomeCoverage.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/3_mtDNAtoAutosomeCoverage.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

# IF YOU WANT, YOU CAN CREATE XLS VERSION OF THE covPerChr BY ADDING COLUMNS
# WITH WOLBACHIA STATUS, AUTOSOMAL NUMBER OF READS AND RATIO 

# PLOT NUMBER OF VARIANTS FOUND IN mtDNA PER DGRP LINE ------------------------
numbVarsPerSample <- apply(mitoVcfGeno, 2, 
                           function(x) length(which(x != '0/0' & x != '.')))
numbVarsPerSamplePlot <- data.frame(table(numbVarsPerSample))

X11()
ggplot(data = numbVarsPerSamplePlot,  aes(x = numbVarsPerSample, y = Freq)) +
  geom_bar(colour = pallete(5)[5], stat = "identity", fill = 'black') +
  xlab('Number of variants') + ylab('Number of lines') +
  mashaGgplot2Theme + theme(panel.grid.minor = element_blank())
dev.print(png, paste0(plotOutDir, '/4_numbOfVarsPerLine.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/4_numbOfVarsPerLine.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

# PLOT MAC --------------------------------------------------------------------
# minor allele count
MAC <- apply(mitoVcfGeno, 1, function(x) calcMacMaf(x)[1])
MAF <- apply(mitoVcfGeno, 1, function(x) calcMacMaf(x)[2])
X11()
hist(MAC, breaks = 100, col = 'black', xlim = c(0, 50), 
     main = 'Minor Allele Count (MAC) for mitochondrial variants')
abline(v = 0.05 * ncol(mitoVcfGeno), col = 'red', lwd = 2)
text(0.5 * max(MAC), 0.85 * max(table(MAC)), 
     paste0('MAC cutoff (5%): ', 0.05 * ncol(mitoVcfGeno), '\n',
            'Variants with MAF > 0.05: ', 
            length(which(MAC > 0.05 * ncol(mitoVcfGeno)))))
dev.print(png, paste0(plotOutDir, '/5_MAChist.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/5_MAChist.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

# PLOT STRUCTURAL TYPE --------------------------------------------------------
# plot according to the structural type
structType <- as.data.frame(table(mitoAnn$StructType))
X11()
ggplot(data = structType, aes(x = Var1, y = Freq)) +
        geom_bar(stat = "identity", 
                 fill = c('red', 'orange', 'grey', 'forestgreen'),
                 col = 'black') +
        geom_text(aes(label = Freq, y = Freq + 10), size = 5) +
        xlab('Structural type') + ylab('Number') + mashaGgplot2Theme +
        ggtitle("Variants by structural type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.print(png, paste0(plotOutDir, '/6_structType.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/6_structType.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

# PRINT FILES FOR CIRCOS ------------------------------------------------------
# commented because i manually modified files with mnps
onlySNPs <- mitoAnn[StructType == 'SNP']
onlySNPs <- data.frame(chr = 'chrM', 
                       start = start(mitoVcfGR[onlySNPs$Variant]),
                       end =  end(mitoVcfGR[onlySNPs$Variant]) + 1,
                       type = onlySNPs$circosType,
                       stringsAsFactors = F)

onlyINDELs <- mitoAnn[StructType == 'INS' | StructType == 'DEL']
refLen <- nchar(as.character(mitoVcfGR[onlyINDELs$Variant]$REF))
altLen <- nchar(as.character(unlist(mitoVcfGR[onlyINDELs$Variant]$ALT)))
onlyINDELs <- data.frame(chr = 'chrM', 
                         start = start(mitoVcfGR[onlyINDELs$Variant]),
                         end = start(mitoVcfGR[onlyINDELs$Variant]) + 
                               ifelse(refLen > altLen, refLen + 1, altLen + 1),
                         type = onlyINDELs$FuncType,
                         stringsAsFactors = F)
onlyINDELs$type <- gsub('^', 'color=', onlyINDELs$type)
onlyINDELs$type <- gsub(' ', '_', onlyINDELs$type)

onlyMNPs <- mitoAnn[StructType == 'MNP']
onlyMNPs <- onlyMNPs[!Variant %in% c('chrM:838_AT/TA', 
                                     'chrM:5960_ATATATTTATATATATATATATAT/TTA')]
for (varID in onlyMNPs$Variant) {
  varRefLen <- nchar(mitoVcfGR[gsub(',.*', '', varID)]$REF)
  varAltLen <- max(sapply(as.character(mitoVcfGR[gsub(',.*', '', 
                                                      varID)]$ALT[[1]]),
                          nchar))
  if (varRefLen == 1 & varAltLen == 1) {
    onlySNPs <- rbind(onlySNPs, 
                      data.frame(chr = 'chrM', 
                                 start = start(mitoVcfGR[gsub(',.*', '', varID)]),
                                 end = start(mitoVcfGR[gsub(',.*', '', varID)]) + 1,
                                 type = 7))
  } else {
    onlyINDELs <- rbind(onlyINDELs, 
                        data.frame(chr = 'chrM', 
                                   start = start(mitoVcfGR[gsub(',.*', '', varID)]),
                                   end = start(mitoVcfGR[gsub(',.*', '', varID)]) + 
                                         max(varRefLen, varAltLen) + 1,
                                   type = 'color=MNP'))
  }
}
onlyINDELs <- rbind(onlyINDELs,
                    data.frame(chr = c('chrM', 'chrM'), start = c(838, 5960),
                               end = c(840, 5977), 
                               type = rep('color=MNP', 2)))

# HISTOGRAM
refLen <- nchar(as.character(mitoVcfGR$REF))
altLen <- sapply(mitoVcfGR$ALT, function(x) max(sapply(unlist(x), nchar)))
varsFreq <- data.frame(chr = 'chrM', 
                       start = start(mitoVcfGR),
                       end = start(mitoVcfGR) + ifelse(refLen > altLen, 
                                                       refLen + 1, altLen + 1),
                       numb = apply(mitoVcfGeno, 1, 
                                    function(x) sum(x != '0/0' & x != '.')))

write.table(onlySNPs, circosSNPs, row.names = F, col.names = F, quote = F, 
            sep ='\t')
write.table(onlyINDELs, circosINDELs, row.names = F, col.names = F, quote = F, 
            sep ='\t')
write.table(varsFreq, circosHIST, row.names = F, col.names = F, quote = F, 
            sep ='\t')
# run circos
#system(paste0('./', plotOutDir, 'mitoVarsLocation/circos.sh'))

# PLOT NUMBER OF VARIANTS PER FUNCTIONAL TYPE ---------------------------------
funcTypeStats <- table(mitoAnn$FuncType)
funcTypeStats <- sort(funcTypeStats, decreasing = T)

X11()
ggplot(data = data.frame(funcTypeStats),
       aes(x = data.frame(funcTypeStats)$Var1,
           y = data.frame(funcTypeStats)$Freq)) +
       geom_bar(stat = "identity",
       fill = rainbowPallete[names(funcTypeStats)]) +
       geom_text(aes(label = Freq, y = Freq + 2), size = 5) +
       xlab('Functional type') + ylab('Number') + mashaGgplot2Theme +
       ggtitle("Variants by functional type") +
       theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.print(png, paste0(plotOutDir, '/7_varsByFuncType.png'), 
          width = pngWidth, height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/7_varsByFuncType.svg'), 
          width = svgWidth, height = svgWidth)
dev.off()

# CHECK, HOW HARMFULL THE HIGH IMPACT VARIANTS ARE ----------------------------
frameShifts <- mitoAnn[FuncType == 'Frameshift']
setkey(frameShifts, Variant)
# reference sequence of mtDNA, dm6
mtDNAseq <- read.fasta(dm6Path, seqtype = "DNA")
mtDNAseq <- mtDNAseq[['chrM']]
frameShifts[, removedAA := apply(frameShifts, 1, getDeletedAAs, mitoVcfGR, 
                                 mtDNAseq, mitoGenes)]

# read variants themselves
mitoVcfAD <- readInAlleleDepth(mitoVcfdm6path, 'dm6')
frameShiftsAD <- do.call("rbind",
                         lapply(frameShifts$Variant, getAllelicDepth,
                                mitoVcfGeno, mitoVcfAD))
frameShiftsAD <- as.data.table(frameShiftsAD)
setkey(frameShiftsAD, Variant)
frameShifts <- merge(frameShifts, frameShiftsAD)

# ANALYSE PATTERNS OF INTERGENIC REPEATS BETWEEN ND3 AND TRNA-A ---------------
repeats <- read.table(repeatsPath, header = T, row.names = 1, sep = '\t',
                      stringsAsFactors = F)
repeats <- repeats[, apply(repeats, 2, max) >= 10]

hist(colSums(repeats), breaks = 100)
range(colSums(repeats))

barplot(apply(repeats, 2, function(x) sum(x > 10)))
repeatsStats <- sort(apply(repeats, 2, function(x) sum(x > 10)))

X11()
par(mar=c(8, 4, 4,2)+0.1)
barplot(repeatsStats, names.arg = c('(AT)3_T_(TA)', '(AT)3_T_(TA)13', 
                                    '(AT)2_T_(TA)9', '(AT)2_T_(TA)10',
                                    '(AT)4_T_(TA)6', '(AT)3_T_(TA)12',
                                    '(AT)2_T_(TA)8', '(AT)4_T_(TA)7', 
                                    '(AT)3_T_(TA)11', '(AT)3_T_(TA)10', 
                                    '(AT)3_T_(TA)9', '(AT)3_T_(TA)5',
                                    '(AT)3_T_(TA)8', '(AT)3_T_(TA)6', 
                                    '(AT)3_T_(TA)7'), 
        las = 2, cex.names = 1, ylab = 'Number of lines', col = '#0E3C52',
        main = 'Number of DGRP lines per intergenic repeat pattern')
dev.print(png, paste0(plotOutDir, '/9_repeatsStats.png'),  width = pngWidth, 
          height = pngWidth)
dev.print(svg, paste0(plotOutDir, '/9_repeatsStats.svg'),  width = svgWidth, 
          height = svgWidth)
dev.off()
