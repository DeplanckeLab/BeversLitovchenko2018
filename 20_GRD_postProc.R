# FILE: 8_GRD_DGRP2_postProc -------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: 
#
# OPTIONS:  none
# REQUIREMENTS:  none
# BUGS: 
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, maria.litovchenko@epfl.ch
# COMPANY:  EPFL, Lausanne, Switzerland
# VERSION:  1
# CREATED:  29.01.2018
# REVISION: 29.01.2018
baseDir <- '~/Desktop/BitBucket/MitoSeq_RPJB_ML_Aug2017/'
setwd(baseDir)
source('2_functions.R')
library(refGenome)

# Inputs ----------------------------------------------------------------------
# folder with all RDS computed for GRDs
rdsDir <- 'results/GRD_DGRP_Rds/'
nuclGRDsVCF <- 'results/GRD_DGRP2/GRDs_DGRP2_NuclVars.annot.vcf'
saveRdsTo <- 'results/GRD_DGRP2/'
plotOutDir <- 'plots/'
circosOutDir <- 'plots/circos_plot/GRD_DGRP2/'

rdsDir <- 'results/GRD_DGRP_Rds_Indep12/'
nuclGRDsVCF <- 'results/GRD_DGRP2_Indep12/GRDs_DGRP2_Indep12_FDR01_NuclVars.annot.vcf'
saveRdsTo <- 'results/GRD_DGRP2_Indep12/'
plotOutDir <- 'plots/'
circosOutDir <- 'plots/circos_plot/GRD_DGRP2_Indep12/'

svgWidth <- 8
pngWidth <- 1000
pAdjustMethod <- c('fdr', 'bonferroni')
pvalSoftest <- 0.1
distLD <- 50

beversDm6Path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
codingEnds <- 14917
# cut-off on MAF
MAFcutoff <- 0.05
# how many genotyped samples are required?
sampsWithGenoCut <- 150
chroms <- c('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chrX', 'chr4')

# gene coordinates, in accordance to the genome to which the alignment 
# was perfromed
currentDir <- getwd()
setwd("~/Documents/RefGen/dm6/")
uscsRefGen <- ensemblGenome() # create USCS object to store the genome
# read GTF file into uscsRefGen object
read.gtf(uscsRefGen, "dm6refGene.srt.gtf")
uscsRefGenCorrds <- as.data.table(getGenePositions(uscsRefGen))
setkey(uscsRefGenCorrds, gene_name)
setwd(currentDir)

# Clusters of DGRPs
dgrpClust <- fread('results/haplotype.v5.txt', header = T, 
                   stringsAsFactors = F)
dgrpClust <- dgrpClust[!Line %in% c('iso-1_rep1', 'Ore-R_rep1', 'Ore-R-rep2',
                                    'w1118_rep1', 'w1118_rep2')]
dgrpClust[, DGRP := gsub('.*_', '', Line)]
addZero <- grepl('^..$', dgrpClust$DGRP)
dgrpClust$DGRP[addZero] <- paste0('0', dgrpClust$DGRP[addZero])
dgrpClust$DGRP <- paste0('DGRP-', dgrpClust$DGRP)
dgrpClust <- dgrpClust[!DGRP %in% c('DGRP-338', 'DGRP-356')]
dgrpClustColors <- c("#DDCC77", '#B11623', "#88CCEE", "#44AA99", "#008C9E", 
                     "#999933", '#0B486B', "#EFBDC0", "#AA4499", "#8C315D",
                     "#998FB8", "#483078", "#117733") #"#661100")
names(dgrpClustColors) <- sort(unique(dgrpClust$Haplotypes_TCS))

# clusters of mitochondrial variants
mitoClust <- readRDS('results/mitoClust.Rds')
mitoClust[, color := dgrpClustColors[Haplotypes_TCSspec]]
mitoClust[Haplotypes_TCSspec == 'MH7']$color <- "#8C315D"
mitoClust[Haplotypes_TCSspec == 'MH8']$color <- "#998FB8"

# Read mitochondrial vcf ------------------------------------------------------
mitoVcfdm6 <- readInMitoVCF(beversDm6Path, 'dm6', codingEnds, 
                            removeRefLines = T, interGen = "merge")
mitoVcfGeno <- mitoVcfdm6[[2]]
# remove MNPs
mitoVcfGeno <- mitoVcfGeno[getVariantStructType(mitoVcfdm6[[1]]) != 'MNP', ]

# there are only 33 variants which pass MAF > 5% cut off. However, many of them
# have the same "configuration" in the sequenced DGRP lines aka clusters. I 
# will scan only 1 variant per cluster, which results into 12 variants
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
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/0", 1, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/1", 3, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("1/0", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("0/1", 2, x))
mitoVcfGeno <- apply(mitoVcfGeno, 2, function(x) gsub("\\.", 0, x))
mitoVcfGeno <- as.matrix(mitoVcfGeno)

# annotation
mitoAnn <- readVariantAnnot(beversDm6Path, 
                            mitoVcfdm6[[1]][rownames(mitoVcfGeno)])
mitoAnn <- simplifyAnnot(mitoAnn)
setkey(mitoAnn, Variant)
setnames(mitoAnn, colnames(mitoAnn), paste0("MitoVar_", colnames(mitoAnn)))
setnames(mitoAnn, 'MitoVar_Variant', 'MitoVar')
setkey(mitoAnn, 'MitoVar')
setkey(mitoClust, 'MitoVar')
mitoAnn <- merge(mitoAnn, mitoClust)

# Read-in RDSs for GRDs -------------------------------------------------------
pvalMatrix <- c()
completeDF <- c()
for (i in 1:nrow(mitoVcfGeno)) {
  print(paste('Started', i, 'at', Sys.time()))
  rdsOneVar <- list.files(rdsDir, pattern = paste0("_", i, ".Rds"), 
                          full.names = T)
  GRDsOneVar <- lapply(rdsOneVar, readRDS)
  GRDsOneVar <- do.call("rbind", GRDsOneVar)
  
  pvalMatrix <- cbind(pvalMatrix, GRDsOneVar$pvalue)
  colnames(pvalMatrix)[i] <- as.character(unique(GRDsOneVar$MitoVar))
  
  completeDF <- rbind(completeDF, GRDsOneVar[, c('NuclVar', 'MitoVar',
                                                 'complete')])
}
rownames(pvalMatrix) <- GRDsOneVar$NuclVar
saveRDS(pvalMatrix, paste0(saveRdsTo, 'pvalMatrix_DGRP2_Indep12.Rds'))
GRDsOneVar <- NULL
completeDF <- as.data.table(completeDF)
saveRDS(completeDF, paste0(saveRdsTo, 'completeDF_DGRP2_Indep12.Rds'))

# Select significant GRDs by Benjamini 2013 procedure -------------------------
pvalMatrix <- readRDS(paste0(saveRdsTo, 'pvalMatrix_DGRP2_Indep12.Rds'))
completeDF <- readRDS(paste0(saveRdsTo, 'completeDF_DGRP2_Indep12.Rds'))

pvalMatrix <- pvalMatrix[, mitoClust[tested == 1]$MitoVar]
signGRDs <- pAdjMultPhenBenjamini2013(pvalMatrix, adjMeth = 'BH', 
                                      pvalCut = pvalSoftest)
colnames(signGRDs) <- c('NuclVar', 'MitoVar', 'pvalue', 'padj')
signGRDs <- as.data.table(signGRDs)
# add information about completness
setkey(signGRDs, 'NuclVar', 'MitoVar')
setkey(completeDF, 'NuclVar', 'MitoVar')
signGRDs <- merge(signGRDs, completeDF, all.x = T)
# add info about whatever or not GRDs pass cutoff on neighbors
signGRDsNeigh <- selectVarsWithSignNeighbors(signGRDs, distLD)
signGRDs[, passNeigh := 0]
setkey(signGRDs, NuclVar)
signGRDs[signGRDsNeigh, ]$passNeigh <- 1
signGRDs[, chr := NULL]
signGRDs[, start := NULL]
signGRDs[, end := NULL]
setkey(signGRDs, 'MitoVar')
signGRDs <- merge(signGRDs, mitoAnn, all.x = T)
saveRDS(signGRDs, paste0(saveRdsTo, 'GRDs_DGRP2_Indep12.Rds'))

# Write files for annotation of nuclear variants in GRDs ----------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12.Rds'))
nuclGRDs <- lapply(as.character(unique(signGRDs$NuclVar)), 
                   function(x) strsplit(x, '_')[[1]])
nuclGRDs <- do.call(rbind, nuclGRDs)
write.table(nuclGRDs, paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_NuclVars.txt'),
            col.names = F, row.names = F, quote = F, sep = '\t')
# run script 35_annotate_GRD_DGRP2.sh

# Add annotation of nuclear variants to GRDs as a separate list ---------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12.Rds'))
signGRDsAnn <- readVariantAnnotNucl(nuclGRDsVCF)
# remove duplicates which are introduces then one gene has several transcripts
signGRDsAnn <- signGRDsAnn[!duplicated(signGRDsAnn), ]
# during vcf-tools extraction some extra variants might arise
signGRDsAnn <- signGRDsAnn[NuclVar %in% signGRDs$NuclVar]

# add information about number of nuclear variants per gene which were tested
signGenes <- unique(signGRDsAnn[inGeneLoc != 'intergenic_region']$Gene)
signGeneCorrds <- uscsRefGenCorrds[gene_name %in% signGenes, ]
pvalMatrix <- readRDS(paste0(saveRdsTo, 'pvalMatrix_DGRP2_Indep12.Rds'))
testedVarsCoords <- rownames(pvalMatrix)
testedVarsCoords <- sapply(testedVarsCoords, strsplit, '_')
testedVarsCoords <- data.table(chr = sapply(testedVarsCoords, 
                                            function(x) x[1]),
                               start = sapply(testedVarsCoords,
                                              function(x) x[2]))
testedVarsCoords[, end := start]
testedVarsCoords <- makeGRangesFromDataFrame(testedVarsCoords)
signGeneCorrdsGR <- makeGRangesFromDataFrame(signGeneCorrds,
                                             keep.extra.columns = T)
ovrl <- as.data.table(findOverlaps(signGeneCorrdsGR, testedVarsCoords))
testedVarsPerGene <- ovrl[,.N, by = queryHits]
testedVarsPerGene[, Gene := signGeneCorrdsGR[queryHits, ]$gene_name]
testedVarsPerGene <- testedVarsPerGene[, -1, with = F]
setnames(testedVarsPerGene, 'N', 'numbOfVarsPerGene')
setkey(testedVarsPerGene, Gene)
setkey(signGRDsAnn, Gene)
signGRDsAnn <- merge(signGRDsAnn, testedVarsPerGene, all.x = T)

# how many nuclear-mito variant chi-sq tests were performed per gene
ovrl[, queryHits := as.character(queryHits)] # otherwise it's not keyable
setkey(ovrl, queryHits)
testsPerGene <- sapply(unique(ovrl$queryHits), 
                       function(x) sum(!is.na(pvalMatrix[ovrl[x, subjectHits],])))
names(testsPerGene) <-  signGeneCorrdsGR$gene_name[as.integer(names(testsPerGene))]
testsPerGene <- data.table(Gene = names(testsPerGene), 
                           testsPerGene = testsPerGene)
setkey(testsPerGene, Gene)     
setkey(signGRDsAnn, Gene)
signGRDsAnn <- merge(signGRDsAnn, testsPerGene, all.x = T)

# add information about gene length
geneLen <- data.table(Gene = uscsRefGenCorrds$gene_name, 
                      geneLen = uscsRefGenCorrds$end - uscsRefGenCorrds$start)
setkey(geneLen, Gene)
setkey(signGRDsAnn, Gene)
signGRDsAnn <- merge(signGRDsAnn, geneLen, all.x = T)

# form final list
signGRDs <- list(signGRDs, signGRDsAnn)
names(signGRDs) <- c('GRDs', 'NuclVarsAnnot')
saveRDS(signGRDs, paste0(saveRdsTo, 'GRDs_DGRP2_Indep12.Rds'))

# Distribution of GRDs over chromosomes vs haplogroups ------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_FDR01.Rds'))
haplGrChr <- data.table(MitoVar = signGRDs$GRDs[passNeigh == 1]$MitoVar,
                        chr = gsub('_.*', '', 
                                   signGRDs$GRDs[passNeigh == 1]$NuclVar))
haplGrChr <- haplGrChr[,.N, by = .(chr, MitoVar)]
mitoClustReduced <- mitoClust[MitoVar %in% signGRDs$GRDs[passNeigh == 1]$MitoVar]
mitoClustReduced[MitoVar == 'chrM:10226_C/T']$Haplotypes_TCSspec <- 'MH5'
mitoClustReduced[MitoVar == 'chrM:2661_C/T']$Haplotypes_TCSspec <- 'MH5+MH6'
mitoClustReduced[MitoVar == 'chrM:12132_C/T']$Haplotypes_TCSspec <- 'chrM:12132_C/T'
setkey(haplGrChr, MitoVar)
setkey(mitoClustReduced, MitoVar)
haplGrChr <- merge(haplGrChr, mitoClustReduced, all = T)
haplGrChr[Haplotypes_TCSspec == 'MH5']$color <- "#0B486B"
haplGrChr[Haplotypes_TCSspec == 'chrM:12132_C/T']$color <- 'grey'
haplGrChr[MitoVar == 'chrM:2661_C/T']$color <- 'black'
clustColours <- data.table(haplGrChr$Haplotypes_TCSspec, haplGrChr$color)
clustColours <- clustColours[!duplicated(clustColours)]
clustColoursV <- clustColours$V2
names(clustColoursV) <- clustColours$V1
haplGrChr[, Haplotypes_TCSspec := factor(Haplotypes_TCSspec,
                                         levels = c(sort(unique(haplGrChr$Haplotypes_TCSspec))[-1],
                                                    'chrM:12132_C/T'))]
ggplot(haplGrChr, aes(x = chr, y = N, fill = Haplotypes_TCSspec)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = clustColoursV) +
  scale_x_discrete(name = "") + ylab("Number of GRDs") +
  guides(fill=guide_legend(title = "Haplotypes")) +
  theme_bw() +
  theme(axis.text.y = element_text(size=15, colour="black"),
        axis.text.x = element_text(colour="black", size=rel(2)),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin=unit(c(.1,.5,0.5,.1),"cm"),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.text=element_text(size = 15),
        legend.title=element_text(size=15)) +
  labs(title = "", x = "", y= "", fill="") +
  scale_y_continuous(limits=c(0,550), breaks=seq(0,550,50), expand=c(0,0))

ggplot(haplGrChr[!Haplotypes_TCSspec %in% c('chrM:12132_C/T', 'MH5+MH6')],
       aes(x = chr, y = N, fill = Haplotypes_TCSspec)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = clustColoursV) +
  scale_x_discrete(name = "") + ylab("Number of GRDs") +
  guides(fill=guide_legend(title = "Haplotypes")) + mashaGgplot2Theme

haplGrChrToPrint <- matrix(haplGrChr$N, nrow = length(unique(haplGrChr$chr)),
                           ncol = length(unique(haplGrChr$Haplotypes_TCSspec)),
                           byrow = T)
colnames(haplGrChrToPrint) <- unique(haplGrChr$Haplotypes_TCSspec) 
rownames(haplGrChrToPrint) <- unique(haplGrChr$chr)
haplGrChrToPrint <- rbind(unique(haplGrChr$MitoVar), haplGrChrToPrint)
write.table(haplGrChrToPrint, '~/Desktop/GRDsPerHaplotypePerChr_passNeigh_1.csv',
            quote = F, sep = '\t', row.names = T, col.names = T)

# Plot info about GRDs and compGRDs: nucl and mito vars -----------------------
plotCounter <- 10
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_FDR01.Rds'))
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_FDR01.Rds'))
signGRDs <- signGRDs$GRDs
# for 0.05 version
signGRDs[, MitoVarAndClust := ifelse(Cluster != '', 
                                     paste0(MitoVar, ' (', Cluster, ')'),
                                     paste0(MitoVar, ' (no cluster)'))] 
for (neighStat in c(0, 1)) {
  for (complete in c(0, 1)) {
    GRDsNeigh <- signGRDs
    if (complete == 1) {
      GRDsNeigh <- signGRDs[complete == 1]
    }
    if (neighStat == 1) {
      GRDsNeigh <- GRDsNeigh[passNeigh == neighStat, ]
    }
    # number of nuclear variants not compatible with a certain mitochondrial 
    # variant
    GRDs_Mito <- GRDsNeigh[,.N, by = MitoVar]
    labelCol <- rainbowPallete[mitoAnn[GRDs_Mito$MitoVar, ]$MitoVar_FuncType]
    GRDs_Mito[, MitoVar := sapply(MitoVar, 
                                  function(x) unique(GRDsNeigh[MitoVar == x, ]$MitoVarAndClust))]
    GRDs_Mito <- GRDs_Mito[order(N), ]
    GRDs_Mito[, MitoVar := factor(MitoVar, levels = MitoVar)]
    
    # number of mito variants not compatible with certain nuclear
    GRDs_Nucl <- GRDsNeigh[,.N, by = NuclVar]
    GRDs_Nucl <- GRDs_Nucl[order(N), ]
    GRDs_Nucl[, NuclVar := factor(NuclVar, levels = NuclVar)]
    
    for (typeToPlot in c('GRDs_Nucl', 'GRDs_Mito')){
      toPlot <- GRDs_Mito
      type <- 'MitoVars'
      if (typeToPlot == 'GRDs_Nucl') {
        toPlot <- GRDs_Nucl
        type <- 'NuclVars'
      }
      plotTitle <- paste('Overview of incompatibilities < 0.05,',
                         'cut off on neighbors =', neighStat)
      fileName <- paste0(plotOutDir, '/',  plotCounter, '_GRDs_DGRP2_Indep12_',
                         type, '_withNeight_', neighStat)
      if (complete == 1) {
        plotTitle <- gsub('of incomp', 'of complete incomp', plotTitle)
        fileName <- gsub('GRDs', 'compGRDs', fileName)
      }
      print(plotGRDbyVar(toPlot, labelCol, ggtitle(plotTitle)))
      #dev.print(png, paste0(fileName, '.png'), width = pngWidth,
      #          height = pngWidth)
      #dev.print(svg, paste0(fileName, '.svg'), width = svgWidth, 
      #          height = svgWidth)
      #dev.off()
      plotCounter <- plotCounter + 1
    }
  }
}

# Figure 4A -------------------------------------------------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_FDR01.Rds'))
signGRDs <- signGRDs$GRDs
GRDsNeigh <- signGRDs[passNeigh == 1]

# number of nuclear variants not compatible with a certain mitochondrial 
# variant
GRDs_Mito <- GRDsNeigh[,.N, by = .(MitoVar, Cluster)]
MitoVarAndClust <- sapply(1:nrow(GRDs_Mito), 
                          function(x) if(grepl('CL', GRDs_Mito$Cluster[x])) {
                            paste0(GRDs_Mito$Cluster[x], '(', 
                                   GRDs_Mito$MitoVar[x], ')')
                          } else {GRDs_Mito$MitoVar[x]})
GRDs_Mito[, MitoVarAndClust := MitoVarAndClust] 
GRDs_Mito <- GRDs_Mito[order(N), ]
setkey(mitoClust, MitoVar)
GRDs_Mito[, colour := mitoClust[GRDs_Mito$MitoVar]$color]
GRDs_Mito[is.na(colour)]$colour <- 'grey'
GRDs_Mito[, Haplotypes_TCSspec := mitoClust[GRDs_Mito$MitoVar]$Haplotypes_TCSspec]
GRDs_Mito[, MitoVar := factor(MitoVar, levels = MitoVar)]

ggplot(GRDs_Mito, aes(x = MitoVar, y = N)) +
  geom_bar(colour = GRDs_Mito$colour, stat = "identity", 
           fill = GRDs_Mito$colour) + coord_flip() +
  geom_text(aes(label = Haplotypes_TCSspec), 
            position = position_dodge(width = 0.9), hjust = -0.25) +
  xlab('mtDNA variant') + ylab('# of incompatible nuclear variants') +
  mashaGgplot2Theme

# Plot info about nuclear variants in GRDs ------------------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_FDR01.Rds'))
plotCounter <- 18
ensembl <- useMart("ensembl", "dmelanogaster_gene_ensembl")

for (neighStat in c(0, 1)) {
  for (complete in c(0, 1)) {
    GRDsNeigh <- signGRDs$GRDs
    if (complete == 1) {
      GRDsNeigh <- GRDsNeigh[complete == 1, ]
    }
    if (neighStat == 1) {
      GRDsNeigh <- GRDsNeigh[passNeigh == neighStat, ]
    }
    GRDsNeighAnn <- signGRDs$NuclVarsAnnot[NuclVar %in% GRDsNeigh$NuclVar]
    
    plotTitle <- paste('Overview of incompatibilities < 0.05, cut off on', 
                       'neighbors =', neighStat)
    if (complete == 1) {
      plotTitle <- gsub('of incom', 'of complete inco', plotTitle)
    }
    
    plotList <- list()
    
    # nuclear variants by structural type
    GRDsNeighAnnStrType <- GRDsNeighAnn[, c('NuclVar', 'StructType')]
    GRDsNeighAnnStrType <- GRDsNeighAnnStrType[!duplicated(GRDsNeighAnnStrType), ]
    GRDsStrType <- GRDsNeighAnnStrType[,.N, by = StructType][order(-N)]
    GRDsStrType[, StructType := factor(StructType, levels = StructType)]
    pl <- ggplot(GRDsStrType, aes(x = StructType,  y = N)) +
          geom_bar(colour = pallete(5)[5], stat = "identity", 
                   fill = pallete(5)[5]) +
          xlab('Structural type') +  ylab('# of variants') +
          ggtitle(plotTitle) + mashaGgplot2Theme
    plotList[[length(plotList) + 1]] <- pl
    
    # nuclear variants by functional type
    GRDsGeneLoc <- GRDsNeighAnn[, c('NuclVar', 'inGeneLoc')]
    GRDsGeneLoc <- GRDsGeneLoc[!duplicated(GRDsGeneLoc), ]
    GRDsGeneLoc <- GRDsGeneLoc[,.N, by = inGeneLoc][order(-N)]
    GRDsGeneLoc[, inGeneLoc := gsub('_variant|_region', '', inGeneLoc)]
    GRDsGeneLoc[, inGeneLoc := factor(inGeneLoc, levels = inGeneLoc)]
    pl <- ggplot(GRDsGeneLoc, aes(x = inGeneLoc, y = N)) +
                 geom_bar(colour = pallete(5)[5], stat = "identity", 
                          fill = pallete(5)[5]) +
            xlab('Functional type') +  ylab('# of variants') +
            ggtitle(plotTitle) + mashaGgplot2Theme + 
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plotList[[length(plotList) + 1]] <- pl
    
    # nuclear variants by gene
    GRDsGene <- GRDsNeighAnn[inGeneLoc != 'intergenic_region']
    GRDsGene <- GRDsGene[, c('NuclVar', 'Gene', 'numbOfVarsPerGene', 
                             'testsPerGene', 'geneLen', 'inGeneLoc')]
    GRDsGene[, inGeneLoc := gsub('_variant.*', '', inGeneLoc)]
    GRDsGene[, inGeneLoc := gsub('UTR.*', 'UTR', inGeneLoc)]
    GRDsGene[, inGeneLoc := gsub("_prime_", "'", inGeneLoc)]
    GRDsGene <- GRDsGene[!duplicated(GRDsGene[, c('NuclVar', 'Gene')]), ]
    byGene <- GRDsGene[,.N, by = Gene]
    byGene <- byGene[order(-N), ]
    byGene[, Gene := factor(Gene, levels = Gene)]
    if(nrow(byGene) > 100) {
      plotTitle <- paste(plotTitle, ', top 100 genes out of', nrow(byGene))
      byGene <- byGene[1:100, ]
    }
    pl <- ggplot(byGene, aes(x = Gene, y = N)) +
          geom_bar(colour = 'black', stat = "identity", fill = 'black') +
          xlab('Gene') +  ylab('# of variants') +
          ggtitle(plotTitle) + mashaGgplot2Theme + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plotList[[length(plotList) + 1]] <- pl
    
    # nuclear variants by gene and gene location
    byGeneAndLoc <- GRDsGene[,.N, by = .(Gene, inGeneLoc)]
    setkey(byGeneAndLoc, Gene)
    # restrict to most interesting genes
    byGeneAndLoc <- byGeneAndLoc[as.character(byGene$Gene), ]
    byGeneAndLoc[, Gene := factor(Gene, levels = unique(Gene))]
    locationsLevels <- byGeneAndLoc[,.N, by = inGeneLoc][order(N), ]$inGeneLoc
    byGeneAndLoc[, inGeneLoc := factor(inGeneLoc, levels = locationsLevels)]
    pl <- ggplot(byGeneAndLoc, aes(x = Gene, y = N, fill = inGeneLoc)) +
          geom_bar(stat = "identity") + xlab('Gene') +  ylab('# of variants') +
          ggtitle(plotTitle) + mashaGgplot2Theme + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plotList[[length(plotList) + 1]] <- pl
    
    # nuclear variants by gene, adjusted by number of comparisons per gene
    byGeneAdj <- GRDsGene[,.N, by = .(Gene, numbOfVarsPerGene, testsPerGene,
                                      geneLen)]
    byGeneAdj <- byGeneAdj[complete.cases(byGeneAdj), ]
    byGeneAdj <- byGeneAdj[, Nadj := N * numbOfVarsPerGene / geneLen]
    byGeneAdj <- byGeneAdj[order(-Nadj), ]
    byGeneAdj[, Gene := factor(Gene, levels = Gene)]
    if(nrow(byGeneAdj) > 100) {
      byGeneAdj <- byGeneAdj[1:100, ]
    }
    pl <- ggplot(byGeneAdj, aes(x = Gene, y = Nadj)) +
          geom_bar(colour = 'black', stat = "identity", fill = 'black') +
          xlab('Gene') + ggtitle(plotTitle) + mashaGgplot2Theme +
          ylab('# of variants adjusted for the 
               # of tested variants and gene length') +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plotList[[length(plotList) + 1]] <- pl
    
    # number of mito variants per gene
    mitoVarsInGRDs <- signGRDs$GRDs[, c('Cluster', 'NuclVar')]
    setnames(mitoVarsInGRDs, 'Cluster', 'MitoVar')
    mitoVarsInGRDs <- mitoVarsInGRDs[!duplicated(mitoVarsInGRDs), ]
    nuclVarsInGRDs <- GRDsNeighAnn[inGeneLoc != 'intergenic_region']
    nuclVarsInGRDs <- nuclVarsInGRDs[, c('NuclVar', 'Gene', 'geneLen')]
    nuclVarsInGRDs <- nuclVarsInGRDs[!duplicated(nuclVarsInGRDs), ]
    nuclVarsInGRDs[, Gene := as.character(Gene)]
    geneAndLen <- nuclVarsInGRDs[, c('Gene', 'geneLen')]
    geneAndLen <- geneAndLen[!duplicated(geneAndLen),]
    setkey(geneAndLen, Gene)
    setkey(mitoVarsInGRDs, NuclVar)
    setkey(nuclVarsInGRDs, Gene)
    mitoVarsPerGene <- sapply(unique(nuclVarsInGRDs$Gene), 
                              function(x) length(unique(mitoVarsInGRDs[nuclVarsInGRDs[x, ]$NuclVar, ]$MitoVar)))
    mitoVarsPerGene <- data.table(Gene = names(mitoVarsPerGene),
                                  NumbOfMitoVars = mitoVarsPerGene,
                                  GeneLen = geneAndLen[names(mitoVarsPerGene), ]$geneLen)
    mitoVarsPerGene <- mitoVarsPerGene[order(-NumbOfMitoVars), ]
    mitoVarsPerGene[, Gene := factor(Gene, unique(Gene))]
    pl <- ggplot(mitoVarsPerGene[1:100], aes(x = Gene, y = NumbOfMitoVars)) +
      geom_bar(colour = 'black', stat = "identity", fill = 'black') +
      xlab('Gene') + ggtitle(plotTitle) + mashaGgplot2Theme +
      ylab('# of incomatible mitochondrial variants') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plotList[[length(plotList) + 1]] <- pl
    
    for (onePlot in plotList) {
      fileName <- paste0(plotOutDir, '/',  plotCounter, '_GRDs_DGRP2_NuclVars',
                         '_withNeight_', neighStat)
      if (complete == 1) {
        fileName <- gsub('GRDs', 'compGRDs', fileName)
      }
      print(onePlot)
      #dev.print(png, paste0(fileName, '.png'), width = pngWidth,
      #          height = pngWidth)
      #dev.print(svg, paste0(fileName, '.svg'), width = svgWidth, 
      #          height = svgWidth)
      #dev.off()
      plotCounter <- plotCounter + 1
    }
    
    GRDgeneNames <- unique(GRDsNeighAnn[inGeneLoc != 'intergenic_region']$Gene)
    GRDgeneNames <- gsub("'", "", GRDgeneNames) # otherwise biomaRt gives error
    GRDsNeighAnnGenesFBids <- getInfoAboutGenes(ensembl, GRDgeneNames)
    for (ont in c('BP', 'MF', 'CC')) {
      GOenr <- GOenricment(GRDsNeighAnnGenesFBids$ensembl_gene_id,
                           getAllFlyGenes()$ensembl_gene_id, ont = ont, 
                           topNodes = 50)
      write.table(GOenr, sep = '\t', row.names = F, quote = F,
                  file = paste0(saveRdsTo, 'GRDs_DGRP2_NuclVars_GOenrich_',
                                ont, '_withNeight_', neighStat, '_complete_',
                                complete, '.csv'))
    }
  }
}

# Figure 4B -------------------------------------------------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_Indep12_FDR01.Rds'))
GRDsNeigh <- signGRDs$GRDs[passNeigh == 1]
GRDsNeighAnn <- signGRDs$NuclVarsAnnot[NuclVar %in% GRDsNeigh$NuclVar]

GRDsGene <- GRDsNeighAnn[inGeneLoc != 'intergenic_region']
GRDsGene <- GRDsGene[, c('NuclVar', 'Gene', 'numbOfVarsPerGene', 
                         'testsPerGene', 'geneLen', 'inGeneLoc')]
GRDsGene[, inGeneLoc := gsub('_variant.*', '', inGeneLoc)]
GRDsGene[, inGeneLoc := gsub('UTR.*', 'UTR', inGeneLoc)]
GRDsGene[, inGeneLoc := gsub("_prime_", "'", inGeneLoc)]
GRDsGene <- GRDsGene[!duplicated(GRDsGene[, c('NuclVar', 'Gene')]), ]
byGene <- GRDsGene[,.N, by = Gene]
byGene <- byGene[order(-N), ]
byGene[, Gene := factor(Gene, levels = Gene)]
if(nrow(byGene) > 100) {
  byGene <- byGene[1:100, ]
}

byGeneAndLoc <- GRDsGene[,.N, by = .(Gene, inGeneLoc)]
setkey(byGeneAndLoc, Gene)
# restrict to most interesting genes
byGeneAndLoc <- byGeneAndLoc[as.character(byGene$Gene), ]
byGeneAndLoc[, Gene := factor(Gene, levels = unique(Gene))]
locationsLevels <- byGeneAndLoc[,.N, by = inGeneLoc][order(N), ]$inGeneLoc
byGeneAndLoc[, inGeneLoc := factor(inGeneLoc, 
                                   levels = c("downstream_gene", "5'UTR", 
                                              "synonymous", "missense", 
                                              "non_coding_exon", "intron",
                                              "splice_region", "3'UTR", 
                                              "upstream_gene"))]
ggplot(byGeneAndLoc, aes(x = Gene, y = N, fill = inGeneLoc)) +
  geom_bar(stat = "identity") + xlab('Gene') + ylab('# of variants') +
  mashaGgplot2Theme + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_manual(values = c("#B9B39E", "#D4C073", "#7797CE", "#306AB2",
                               "#3D475D", "#A18D2A", "#786B3F", "#C4C5C9",
                               "#333B4E"))

# Plot variants and lines of interest to check. Only complete GRDs ------------
# for GRDs which were tested in the lab
plotCounter <- 30
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2.Rds'))
signGRDs <- signGRDs$GRDs[complete == 0]
dgrp2dm6Geno <- readRDS('dgrp2dm6Geno.Rds')
setkey(dgrp2dm6Geno, 'ID')
COI <- 'CL:4,'
LOI <- c('DGRP-153', 'DGRP-287', 'DGRP-306', 'DGRP-189', 'DGRP-318',
         'DGRP-340', 'DGRP-365', 'DGRP-801', 'DGRP-882')
dgrp2dm6Geno <- dgrp2dm6Geno[, c('ID', LOI), with = F] # reduces memory

for (neighStat in c(F, T)) {
  GRDsNeigh <- signGRDs
  if (neighStat) {
    GRDsNeigh <- signGRDs[passNeigh == neighStat, ]
  }
  
  NuclVOI <- unique(GRDsNeigh[MitoVar_cluster == COI, ]$NuclVar)
  MitoVOI <- unique(GRDsNeigh[MitoVar_cluster == COI, ]$MitoVar)
  
  nuclGeno <- as.matrix(dgrp2dm6Geno[NuclVOI, ])
  rownames(nuclGeno) <- nuclGeno[, 1]
  nuclGeno <- nuclGeno[, -1]
  NuclVOIclust <- findVarClusters(nuclGeno)
  NuclVOIgeno <- nuclGeno[names(sort(NuclVOIclust)), ]
  matrToVisual <- rbind(NuclVOIgeno, mitoVcfGeno[as.character(MitoVOI), LOI])
  
  if (nrow(matrToVisual) != 0) {
    matrToVisual <- melt(matrToVisual)
    colnames(matrToVisual) <- c('Variant', 'DGRP', 'Status')
    print(ggplot(matrToVisual, aes(x = DGRP, y = Variant)) +
           geom_raster(aes(fill = as.numeric(Status))) +
            scale_x_discrete(drop = FALSE) +
            scale_y_discrete(drop = FALSE) + mashaGgplot2Theme +
            theme(axis.text.x = element_text(angle = 270, hjust = 0),
                  aspect.ratio = 1, legend.position = "none") + 
            ggtitle(paste('Overview of complete incompatibilities,',
                           'DGRPs of interest\n < 0.05,', 
                           'cut off on neighbors =', neighStat)))
    fileName <- paste0(plotOutDir, '/', plotCounter, '_compGRDs_DGRPtotest_',
                       'withNeight_', neighStat)
    dev.print(png, paste0(fileName, '.png'), width = pngWidth, height = pngWidth)
    dev.print(svg, paste0(fileName, '.svg'), width = svgWidth, height = svgWidth)
    dev.off()
    plotCounter <- plotCounter + 1
  }
}

# Sxl gene with the most GRDs -------------------------------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_FDR01.Rds'))
NuclVOI <- unique(signGRDs$NuclVarsAnnot[Gene  == 'Sxl']$NuclVar)
MitoVOI <- unique(signGRDs$GRDs[NuclVar %in% NuclVOI]$MitoVar)
dgrp2dm6Geno <- readRDS('dgrp2dm6Geno.Rds')
setkey(dgrp2dm6Geno, 'ID')
nuclGeno <- as.matrix(dgrp2dm6Geno[NuclVOI, ])
rownames(nuclGeno) <- nuclGeno[, 1]
nuclGeno <- nuclGeno[, -1]
NuclVOIclust <- findVarClusters(nuclGeno)
#NuclVOIclust <- nuclGeno
NuclVOIgeno <- nuclGeno[names(sort(NuclVOIclust)), ]
matrToVisual <- rbind(NuclVOIgeno, mitoVcfGeno[as.character(MitoVOI), ])
matrToVisual <- t(matrToVisual)
matrToVisual1 <- apply(matrToVisual, 1, as.integer)
rownames(matrToVisual1) <- colnames(matrToVisual)
matrToVisual <- matrToVisual1
matrToVisual[matrToVisual == 0] <- NA

matrToVisual <- melt(matrToVisual)
colnames(matrToVisual) <- c('Variant', 'DGRP', 'Status')
matrToVisual$Status[matrToVisual$Status == 0] <- NA
print(ggplot(matrToVisual, aes(x = DGRP, y = Variant)) +
        geom_raster(aes(fill = as.numeric(Status))) +
        scale_x_discrete(drop = FALSE) +
        scale_y_discrete(drop = FALSE) + mashaGgplot2Theme +
        theme(axis.text.x = element_text(angle = 270, hjust = 0),
              aspect.ratio = 1, legend.position = "none") + 
        ggtitle('Overview of incompatibilities for Sxl, introns'))

my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green
heatmap.2(apply(matrToVisual, 1, as.integer),
          #cellnote = mat_data,  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="col",     # only draw a row dendrogram
          Rowv="NA")            # turn off column clustering

# Write file for circos ------------------------------------------------------
signGRDs <- readRDS(paste0(saveRdsTo, 'GRDs_DGRP2_FDR01.Rds'))
for (neighStat in c(F, T)) {
  GRDsNeigh <- signGRDs$GRDs
  if (neighStat) {
    GRDsNeigh <- signGRDs$GRDs[passNeigh == neighStat, ]
  }
  writeGRDsCircos(GRDsNeigh,
                  paste0(circosOutDir, 'GRDs_withNeight_', 
                         neighStat, '.txt'))
  writeGRDsCircos(GRDsNeigh[complete == 1],
                  paste0(circosOutDir, 'compGRDs_withNeight_', 
                         neighStat, '.txt'))
}

# EXTRA : READ GENOTYPE INFO FOR DGRP2 ----------------------------------------
# put here all nuclear variants with GRDs
dgrp2dm6Geno <- c()
for (chr in chroms) {
  print(chr)
  # read-in dgrp genotype info 
  #dgrp2dm6GenoChr <- fread(paste0('dgrp2_dm6_BeversLinesOnly_MAF005_', chr, 
  #                                '.GT.FORMAT'))
  #addZero <- grepl('DGRP-..$', colnames(dgrp2dm6GenoChr))
  #setnames(dgrp2dm6GenoChr, colnames(dgrp2dm6GenoChr)[addZero],
  #         gsub('DGRP-', 'DGRP-0', colnames(dgrp2dm6GenoChr)[addZero]))
  #dgrp2dm6GenoChr <- cbind(ID = paste0(dgrp2dm6GenoChr$CHROM, "_",
  #                                     dgrp2dm6GenoChr$POS), dgrp2dm6GenoChr)
  #dgrp2dm6GenoChr[, CHROM := NULL]
  #dgrp2dm6GenoChr[, POS := NULL]
  #setkey(dgrp2dm6GenoChr, ID)
  #dgrp2dm6Geno <- rbind(dgrp2dm6Geno, dgrp2dm6GenoChr)
}
#saveRDS(dgrp2dm6Geno, 'dgrp2dm6Geno.Rds')
dgrp2dm6Geno <- readRDS('dgrp2dm6Geno.Rds')
# convert to matrix
dgrp2dm6Geno <- as.matrix(dgrp2dm6Geno)
rownames(dgrp2dm6Geno) <- dgrp2dm6Geno[, 1]
dgrp2dm6Geno <- dgrp2dm6Geno[, -1]