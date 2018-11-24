# FILE: 4_heteroplasmy --------------------------------------------------------
#
# USAGE: 
#
# DESCRIPTION: summarizes results of heteroplasmy calculations
#
#  OPTIONS:  none
#  REQUIREMENTS:  none
#  BUGS: --
#  NOTES:  python script from Mitochondrial heteroplasmy in vertebrates using 
#          ChIP-sequencing data
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
# result vcf, dm6
# annotated vcf file, mapping on dm6
mitoVcfdm6path <- 'results/Bevers_dm6_chrM_GenotypeGVCF.flrt.annot.vcf'
# annotation of the genome to which data were aligned
dm6MitoAnnotPath <- 'dm6_annotation/genes_CUSTOM_MANUAL_CHECK.gtf'

# path to the folder with heteroplasmy results
hetepoplDir <- 'results/heteroplasmy'

# FUNCTIONS -------------------------------------------------------------------
readInHeteroplams <- function(pathToFile, mitoGR) {
  hetInf <- fread(pathToFile, sep = ',', header = T)
  hetInf <- hetInf[Position < 14917]
  if (nrow(hetInf) != 0) {
    hetInf[, `:=`(Species = NULL)]
    # converting to GRanges
    setnames(hetInf, 'Position', 'start')
    hetInf[, end := start]
    hetInf[, chr := 'chrM']
    hetInf <- makeGRangesFromDataFrame(hetInf, keep.extra.columns = T)
    # assign heteroplamy status to the existing variants
    ovrl <- as.data.frame(findOverlaps(hetInf, mitoGR))
    if (length(unique(ovrl[, 1])) != length(hetInf)) {
      warning('Error: heteroplasmy was detected where variant was not')
      print(hetInf[which(!(1:length(hetInf)) %in% ovrl[, 1])])
    } 
    result <- data.table(DGRP = gsub('_BWA_dm6', '',
                                     hetInf[ovrl[, 1]]$Individual),
                         Variant = names(mitoGR)[ovrl[, 2]], 
                         HetRatio = round(100 * hetInf[ovrl[, 1]]$Ratio, 2),
                         Bases = hetInf[ovrl[, 1]]$Bases)
    result <- result[!duplicated(result[, .(DGRP, Variant)]), ]
  } else {
    result <- NULL
  }
}

# READ IN VCF & ANNOTATION ----------------------------------------------------
mitoVcfdm6 <- readInMitoVCF(mitoVcfdm6path, 'dm6', codingEnds, 
                            removeRefLines = F, interGen = "merge")
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

# Read-in heteroplasmy files --------------------------------------------------
allhetFiles <- list.files(hetepoplDir, '.het', full.names = T)
hetInf <- lapply(allhetFiles, readInHeteroplams, mitoVcfGR)
hetInf <- do.call(rbind, hetInf)
hetInf[, Bases := NULL]
# convert it into matrix
hetMatrix <- matrix(NA, nrow = nrow(mitoVcfGeno), ncol = ncol(mitoVcfGeno))
rownames(hetMatrix) <- rownames(mitoVcfGeno)
colnames(hetMatrix) <- colnames(mitoVcfGeno)
for (variant in rownames(hetMatrix)) {
  for (dgrp in colnames(hetMatrix)) {
    hetOneVarOneLine <- hetInf[DGRP == dgrp & Variant == variant, ]
    varDetect <- mitoVcfGeno[variant, dgrp]
    # if variant was not detected, and heteroplasmy was - that's bad
    if (varDetect %in% c('0/0', '.') & nrow(hetOneVarOneLine) != 0) {
      stop(paste('Variant', variant, 'WAS NOT detected in', dgrp, ', but',
                 'heteroplamy was'))
    } else {
      if (nrow(hetOneVarOneLine) != 0) {
        hetMatrix[variant, dgrp] <- hetOneVarOneLine$HetRatio
      }
    }
  }
}

# Plots -----------------------------------------------------------------------
numbOfHetPerLine <- apply(hetMatrix, 2, 
                          function(x) sum(!is.na(x)))
numbOfHetPerLine <- data.frame(DGRP = names(numbOfHetPerLine),
                               NumbOfHet = numbOfHetPerLine,
                               NumbOfVars = apply(mitoVcfGeno, 2,
                                                  function(x) sum(x != '0/0' &
                                                                  x != '.')))
dfToPlot <- as.data.frame(table(numbOfHetPerLine$NumbOfHet))
ggplot(dfToPlot, aes(x = Var1, y = Freq)) +
  geom_bar(colour = 'black', stat = "identity", fill = 'black') +
  xlab('# of heteroplasmic variants') +  ylab('# of DGRP lines') +
  ggtitle('Number of heteroplamic variants per DGRP line') + 
  mashaGgplot2Theme

ggplot(hetInf[grepl('DGRP', DGRP), ], aes(x = HetRatio)) + 
  geom_histogram(binwidth = 2, colour = 'black', fill = 'black') +
  xlab('% of reads per minor allele in heteroplamic site') +  
  ylab('# of heteroplasmic sites') +
  ggtitle('Percentage  of reads per minor allele in heteroplamic sites') +
  mashaGgplot2Theme

# Files for circos ------------------------------------------------------------
GRforCircos <- mitoVcfGR[unique(hetInf[grepl('DGRP', DGRP), ]$Variant)]
GRforCircos <- GRforCircos[order(start(GRforCircos))]
GRforCircos <- data.frame(chr = 'chrM', start = start(GRforCircos),
                          end = end(GRforCircos) + 1, 1)
write.table(GRforCircos, sep = '\t', col.names = F, row.names = F, quote = F, 
            file = 'plots/circos_plot/mitoVarsLocation/HETEROPLASM_withRefLines.txt')
# Comparison with Haag-Liautard 2008 ------------------------------------------
# Download table 3 and put it in bed format
# do cross map:
# CrossMap.py bed dm3ToDm6.over.chain.gz \
#             comparison_with_Haag-Liautard/Haag-Liautard_Table_3.csv \
#             comparison_with_Haag-Liautard/Haag-Liautard_Table_3_dm6.csv
# Add manually columns 

# get all heteroplasmic locations for our data and merge with annotation
mitoAnn[, Heteroplasm := apply(hetMatrix, 1, function(x) sum(!is.na(x)) != 0)]
colnames(mitoAnn) <- paste0('Bevers_', colnames(mitoAnn))
mitoAnn[, Position := as.integer(gsub('_.*', '', 
                                      gsub('chrM:', '', Bevers_Variant)))]

# read-in Haag data
HaagHet <- fread('comparison_with_Haag-Liautard/Haag-Liautard_Table_3_dm6.csv',
                 header = T)
HaagHetGR <- makeGRangesFromDataFrame(HaagHet, keep.extra.columns = T)
HaagHet <- HaagHet[, -2:-1, with = F]
colnames(HaagHet) <- paste0('Haag_', colnames(HaagHet))
colnames(HaagHet)[1] <- 'Position'

setkey(mitoAnn, Position)
setkey(HaagHet, Position)

haagBeversHet <- merge(mitoAnn, HaagHet, all = T)
# complete overlap 
compOverls <- haagBeversHet[!is.na(Bevers_Variant) & !is.na(Haag_Line)]
# I insert + 1 for 5960 - 5972 variant
message(paste("Out of", nrow(mitoAnn), "Bevers-variants", 
              length(unique(compOverls$Bevers_Variant)) + 1, 
              "overlap with Haag variants(", length(unique(HaagHet$start)),
              "variants in total):", 
              paste(unique(compOverls$Bevers_Variant), collapse = ','),
              ", chrM:5960_ATATATTTATATATATATATATAT/TTA"))
write.table(haagBeversHet, '~/Desktop/haagBeversHet.csv', col.names = T,
            row.names = F, quote = F, sep = '\t')

# find closest to the Haag variants our variant
View(distanceToNearest(HaagHetGR, mitoVcfGR))
