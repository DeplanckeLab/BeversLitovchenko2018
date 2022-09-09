# FILE: 3_de_novo_NUMTs.R -----------------------------------------------------
#
# DESCRIPTION: Rscript to identify NUMTs de novo from short read sequencing.
#              The sript can be launched in command line.
# USAGE: run in terminal 
#                       Rscript --vanilla 3_de_novo_NUMTs.R -h 
#        to see description of all arguments.
# OPTIONS:  
# EXAMPLE: 
# REQUIREMENTS: R 4.0.2, BSgenome.Dmelanogaster.UCSC.dm6, data.table, ggplot2,
#               Gviz, msa, Rsamtools
#               
# BUGS: --
# NOTES:  ---
# AUTHOR:  Maria Litovchenko, m.litovchenko@ucl.ac.uk
# COMPANY:  UCL, London, the UK
# VERSION:  2
# CREATED:  13.07.2019
# REVISION: 04.09.2022

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(BSgenome.Dmelanogaster.UCSC.dm6))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(Rsamtools))

# Functions -------------------------------------------------------------------
#' printArgs
#' @description Prints submitted to the script arguments as messages.
#' @author Maria Litovchenko
#' @param argsList named list representing submitted to the script arguments
#' @return void
printArgs <- function(argsList) {
  message("Submitted arguments:")
  for (argName in names(argsList)) {
    oneArg <- argsList[[argName]]
    if (length(oneArg) > 1 & !is.null(names(oneArg))) {
      msg <- paste(apply(data.table(names(oneArg), oneArg), 1, paste, 
                         collapse = ' - '), collapse = ',')
    } else {
      msg <- paste(oneArg, collapse = ', ')
    }
    message(argName, ':', msg)
  }
}

#' cigarToVector
#' Converts cigar string to named vector, where values of the vector - number 
#' of bases and names of the vector - cigar codes
#' @param cigar cigar string
#' @return named vector
cigarToVector <- function(cigar) {
  # get cigar codes
  letterCodes <- strsplit(gsub("\\d", "", cigar), '')[[1]]
  # get number of bases per cigar code
  numbsOfBases <- strsplit(cigar, paste(unique(letterCodes), 
                                        collapse = '|'))[[1]]
  numbsOfBases <- as.integer(numbsOfBases)
  names(numbsOfBases) <- letterCodes
  numbsOfBases
}

#' selectSoftClipped
#' @description Selects softclipped reads suitable for NUMTs discovery
#' @param bamData data read from BAM file, result of scanBam function
#' @param overhangLen minimal length of soft clipped and matched part of a read
selectSoftClipped <- function(bamData, overhangLen) {
  proceed <- F
  
  # select against reads having Ns inside them
  scBam <- !grepl('N', bamData$seq)
  if (sum(scBam) > 0) {
    proceed <- T
    scBam <- lapply(bamData, function(x) x[scBam])
  } else {
    proceed <- F
    return(NULL)
  }
  
  # remove reads which contain something else except soft-clipped and match
  # or there are more than 1 segment soft-clipped / matched
  if (proceed) {
    SCMread <- gsub('\\d', '', scBam$cigar) %in% c('SM', 'MS')
    if (sum(SCMread) > 0) {
      proceed <- T
      scBam <- lapply(scBam, function(x) x[SCMread])
    } else {
      proceed <- F
      return(NULL)
    }
  }
  
  # select reads which have both S and M parts longer than overhang
  if (proceed) {
    smBases <- sapply(scBam$cigar, function(x) strsplit(x, 'S|M'))
    smBases <- lapply(smBases, function(x) as.integer(x))
    smBases <- sapply(smBases, function(x) all(x >= overhangLen))
    if (sum(smBases) > 0) {
      scBam <- lapply(scBam, function(x) x[smBases])
      proceed <- T
      return(scBam)
    } else {
      proceed <- F
      return(NULL)
    }
  }
}

mergeBamObj <- function(bamOne, bamTwo) {
  result <- lapply(names(bamOne), function(x) c(bamOne[[x]], bamTwo[[x]]))
  names(result) <- names(bamOne)
  result
}

#' getReadBreakPoint
#' Calculates breakpoint for soft clipped reads
#' @param samRecord sam record of the read, which have only match and
#'                  soft-clipped
#' @note works only for reads which have one part matching and another soft 
#'       clipped
#' @return vector containg break point (absolute coordinate in bp), string with
#'         soft clipped bases, string with matching bases
getReadBreakPoint <- function(samRecord) {
  # get cigar code - will help to determine read direction as well
  cigarCode <- cigarToVector(samRecord$cigar)
  readType <- gsub('\\d', '', samRecord$cigar)
  
  # the START of the read is the position there MATCH STARTS
  if (names(cigarCode)[1] == 'M') {
    # so, in case of match first it's easy
    breakPoint <- as.integer(samRecord$pos) + cigarCode[1] - 1
    # get matched bases
    mBases <- substr(samRecord$seq, 1, cigarCode[1])
    # get soft clipped bases
    scBases <- substr(samRecord$seq, cigarCode[1] + 1, nchar(samRecord$seq))
  } 
  if (names(cigarCode)[1] == 'S') {
    # but in the other case, read starts from not aligned bases
    breakPoint <- as.integer(samRecord$pos)
    # get soft clipped bases
    scBases <- substr(samRecord$seq, 1, cigarCode[1])
    # get matched bases
    mBases <- substr(samRecord$seq, cigarCode[1] + 1, nchar(samRecord$seq))
  }
  
  result <- c(readType, as.integer(breakPoint), scBases, mBases)
  names(result) <- c("type", "breakPoint", "scBases", "mBases")
  result
}

#' getBreakPoint
#' Calculates break point for every read
#' @param samDT data table which contains sam reads
#' @return data table which contains sam reads with added columns breakPoint, 
#'         scBases, mBases: break point (absolute coordinate in bp), 
#'         string with soft clipped bases, string with matching bases
getBreakPoint <- function(bamObj) {
  numbReads <- length(bamObj$qname)
  allBreaks <- sapply(1:numbReads,
                      function(y) getReadBreakPoint(lapply(bamObj,
                                                           function(x) x[y])))
  allBreaks <- as.data.table(t(allBreaks))
  allBreaks[, breakPoint := as.integer(breakPoint)]
  bamObj[['readType']] <- allBreaks$type
  bamObj[['breakPoint']] <- allBreaks$breakPoint
  bamObj[['scBases']] <- allBreaks$scBases
  bamObj[['mBases']] <- allBreaks$mBases
  
  # remove reads which have 1 as a break-point, most likely they have it 
  # because of circularity of mito genome
  result <- lapply(bamObj, function(x) x[bamObj$breakPoint != 1])
  result
}

#' getReadsPerBreakPoint
#' Get number of reads which have break point at a certain position of a genome
#' @param samDT data table which contains sam reads
#' @param minCov minimal number of reads which should have break point at a 
#'               certain bp so it would be preserved
#' @return data table with columns bp and numbReads, bp = base pair where beak
#'         is detected, numbReads = number of reads having break point at this 
#'         bp
getReadsPerBreakPoint <- function(bamObj, minCov = 20) {
  numbReadsPerBreak <- data.table(table(bamObj$breakPoint))
  colnames(numbReadsPerBreak) <- c('bp', 'numbReads')
  numbReadsPerBreak[, bp := as.integer(bp)]
  numbReadsPerBreak <- numbReadsPerBreak[bp != 1 & numbReads >= minCov]
  numbReadsPerBreak
}

#' mostCommonRead
#' Finds most common read amoung given ones, allows mismatches
#' @param reads array of strings (reads)
#' @param readLenCut to which length cut reads - should be the length of the
#'                   shortest reads
#' @param from from which side of the soft clipped bases make cut: from tail,
#'             if it's SM read and start, if it's MS read
#' @param maxMismatch maximum number of allowed mismatches
#' @return data table with columns mainRead, N, suppReads, there mainRead is a 
#'         common read(pattern), N - number of reads with this string pattern
#'         and also supplementary reads, suppReads all supplementary 
#'         reads(patterns)
mostCommonRead <- function(reads, readLenCut = 20, from = 'tail',
                           maxMismatch = 2) {
  # first, cut to the desired length
  if (from == 'start') {
    cutReads <- sapply(reads, substr, 1, readLenCut)
  }
  if (from == 'tail') {
    cutReads <- reverse(sapply(reverse(reads), substr, 1, readLenCut))
  }
  names(cutReads) <- NULL
  
  # 
  readCount <- as.data.table(table(cutReads))
  colnames(readCount) <- c('mainRead', 'N')
  uniqReads <- readCount[order(N, decreasing = T)]$mainRead
  setkey(readCount, mainRead)
  
  # get closest strings to account for the possibility of mistake in read
  readsDist <- adist(uniqReads)
  colnames(readsDist) <- uniqReads
  rownames(readsDist) <- uniqReads
  if (nrow(readsDist) !=  1) {
    closeReads <- data.table(mainRead = character(), suppReads = character())
    for (i in 1:(nrow(readsDist) - 1)) {
      colsToScore <- c((i + 1):ncol(readsDist))
      closeStr <- i + which(readsDist[i, colsToScore] <= maxMismatch)
      if (!identical(closeStr, integer(0))) {
        closeReads <- rbind(closeReads, 
                            data.table(mainRead = uniqReads[i], 
                                       suppReads = uniqReads[closeStr]))
      }
    }
    setkey(closeReads, mainRead)
    
    readCountMain <- readCount[!mainRead %in% 
                                 setdiff(closeReads$suppReads,
                                         closeReads$mainRead)]
    countFromSupp <- apply(readCountMain, 1, 
                           function(x) sum(readCount[closeReads[x['mainRead']]$suppReads]$N))
    countFromSupp[is.na(countFromSupp)] <- 0
    suppReads <- apply(readCountMain, 1, 
                       function(x) closeReads[x['mainRead']]$suppReads)
    readCountMain[, fromSupp := countFromSupp]
    readCountMain[, suppReads := suppReads]
    readCountMain[, N := N + fromSupp]
  } else {
    readCountMain <- readCount
    readCountMain[, fromSupp := 0] 
    readCountMain[, suppReads := NA]
  }
  
  readCountMain <- readCountMain[order(N, decreasing = T), ]
  readCountMain
}

#' getMostCommonReadPerBreakPoint
#' Wrapper around mostCommonRead read function to apply it to data table
#' @param bamObj bam object
#' @param readsPerBreakDT data table, result of getReadsPerBreakPoint
#' @param minNumbReads minimal number of reads which should have break point at
#'                     a certain bp so it would be preserved
#' @param readLenCut to which length cut reads - should be the length of the
#'                   shortest reads
#' @param from from which side of the soft clipped bases make cut: from tail,
#'             if it's SM read and start, if it's MS read
#' @param maxMismatch maximum number of allowed mismatches
#' @return data table with columns mainRead, N, suppReads, there mainRead is a 
#'         common read(pattern), N - number of reads with this string pattern
#'         and also supplementary reads, suppReads all supplementary 
#'         reads(patterns)
getMostCommonReadPerBreakPoint <- function(bamObj, readsPerBreakDT, 
                                           minNumbReads = 20, readLenCut = 20, 
                                           from = 'tail', maxMismatch = 2) {
  result <- lapply(readsPerBreakDT$bp, 
                   function(x) mostCommonRead(bamObj$scBases[bamObj$breakPoint == x], 
                                              readLenCut, from, maxMismatch))
  # add percentage of reads falling into main reads
  lapply(1:nrow(readsPerBreakDT), 
         function(x) result[[x]][, perc := N / readsPerBreakDT$numbReads[x]])
  lapply(1:nrow(readsPerBreakDT), function(x) result[[x]][, perc := 100*perc])
  # add percentage of the second common read
  lapply(result, function(x) x[, secondMainRead := c(x$mainRead[-1], NA)])
  lapply(result, function(x) x[, secondMainReadN := c(x$N[-1], NA)])
  lapply(result, function(x) x[, secondMainReadPerc := c(x$perc[-1], NA)])
  # add bp
  lapply(1:nrow(readsPerBreakDT), 
         function(x) result[[x]][, bp := readsPerBreakDT$bp[x]])
  result <- do.call(rbind, lapply(result, function(x) x[1]))
  result <- result[N >= 20]
  result
}

barplotNconcensus <- function(breakManyReadsDT, allBreakDT) {
  toPlot <- breakManyReadsDT[, c("bp", "N")]
  toPlot[, type := "Most common read"]
  toAdd <- breakManyReadsDT[, c("bp", "secondMainReadN")]
  setnames(toAdd, "secondMainReadN", "N")
  toAdd[, type := "Second read"]
  toPlot <- rbind(toPlot, toAdd)
  toAdd <- data.table(N = allBreakDT[bp %in% toPlot$bp]$numbReads -
                        toPlot[,.(sum(N, na.rm = T)), by = bp]$V1,
                      type = "Rest")
  toAdd <- toAdd[, N := ifelse(N < 0, 0, N)]
  toAdd <- toAdd[, bp := toPlot$bp]
  toPlot <- rbind(toPlot, toAdd)
  toPlot[, bp := as.factor(bp)]
  toPlot[, type := factor(type, 
                          c('Rest','Second read',
                            'Most common read'))]
  ggplot(data = toPlot, aes(x = bp, y = N, fill = type)) +
    geom_bar(stat = "identity") + theme_classic(base_size = 12) +
    ylab('Number of reads')
}

#' compareMainReads
#' Compares main reads (nuclear parts) of the putative NUMTs which breakpoints
#' are close
#' @param breakPointReads data table, result of getMostCommonReadPerBreakPoint
#' @param type SM and MS
#' @param minDist integer, minimal distance between breakpoints so they would 
#'                be counted as overlapping, usually = overhangBP
#' @param minStrDist minimal lewingstein difference between strings
#' @return warning message
compareMainReads <- function(breakPointReads, minDist, minStrDist) {
  # calculate distance between the breakpoints
  bpDist <- breakPointReads$bp - 
    c(0, breakPointReads$bp[-length(breakPointReads$bp)])
  # get the one close to each other
  closeToEachOther <- which(bpDist < minDist)
  closeToEachOther <- closeToEachOther[closeToEachOther != 1]
  # get string distance between main reads
  smStrDist <- sapply(closeToEachOther, 
                      function(x) stringDist(c(breakPointReads[x]$mainRead,
                                               breakPointReads[x - 1]$mainRead)))
  # get ones which are close
  mainReadsSimilar <- which(smStrDist <= minStrDist)
  mainReadsSimilar <- closeToEachOther[mainReadsSimilar]
  # report
  message('WARNING: following nuclear parts of NUMTs are highly similar')
  sapply(mainReadsSimilar, 
         function(x) message(paste(breakPointReads[x - 1]$bp,
                                   breakPointReads[x - 1]$mainRead, ' - ',
                                   breakPointReads[x]$mainRead,
                                   breakPointReads[x]$bp, '| Dist(bp):',
                                   breakPointReads[x]$bp - 
                                     breakPointReads[x - 1]$bp)))
}

#' writeToSam
#' Writes reads to sam file
#' @param samDT data table with sam reads
#' @param fileName name of the output file
#' @return sam file
writeToSam <- function(bamObj, fileName) {
  write("@HD	VN:1.5	SO:coordinate", paste0(fileName, '.sam'), 1, F)
  write("@SQ	SN:chrM	LN:19524", paste0(fileName, '.sam'), 1, T)
  
  # though 1st cigar I determine if I would need to move start
  firstCigar <- names(cigarToVector(bamObj$cigar[1]))[1]
  if (firstCigar == 'S') {
    moveStartBy <- sapply(bamObj$cigar, function(x) cigarToVector(x)[1])
  } else {
    moveStartBy <- 0
  }
  
  toPrint <- data.table(qname = bamObj$qname, flag = bamObj$flag,
                        chr = 'chrM', pos = bamObj$pos - moveStartBy,
                        mapq = bamObj$mapq, 
                        cigar = gsub('S', 'X', bamObj$cigar),
                        x = '=', matePos = bamObj$mpos,
                        insSize = bamObj$isize, 
                        read = as.character(bamObj$seq),
                        qual = as.character(bamObj$qual))
  toPrint[is.na(toPrint)] <- 0
  toPrint <- toPrint[pos > 0]
  write.table(toPrint, paste0(fileName, '.sam'), col.names = F,
              row.names = F, sep = '\t', quote = F, append = T)
}

#' writeToBam
#' Writes reads to Bam file
#' @param samDT data table with aligned reads
#' @param fileName name of the output file
#' @note samtools are required
#' @return bam file
writeToBam <- function(bamObj, fileName) {
  writeToSam(bamObj, fileName)
  system(paste("samtools view -bS", paste0(fileName, '.sam'), ">", 
               paste0(fileName, '.bam')))
  file.remove(paste0(fileName, '.sam'))
  system(paste("samtools sort", paste0(fileName, '.bam'), ">", 
               paste0(fileName, '.srt.bam')))
  file.remove(paste0(fileName, '.bam'))
  system(paste("samtools index", paste0(fileName, '.srt.bam')))
}

#' mergeOverlapNUMTs
#' Performs merging of the main reads corresponding to putative NUMTs in case
#' there is a significant overlap (>= minOvrl) between the coordinates of the
#' reads and substrings of the reads matches with number of mismahes allowed
#' (maxMismatch)
#' @param inGR input GRanges object
#' @param maxMismatch maximum number of mismatches allowed
#' @param minOvrl minimal overlap between the main reads
#' @return GRanges with merged main NUMTs
#' smBreakPointReads <- readRDS('smBreakPointReads.Rds')
#' smBreakPointReads[, start := bp - overhangBP]
#' smBreakPointReads[, end := bp]
#' smBreakPointReads[, chr := 'chrM']
#' smBreakPointReadsGR <- makeGRangesFromDataFrame(smBreakPointReads,
#'                                                 keep.extra.columns = T)
#' smBreakPointReadsGR <- mergeOverlapNUMTs(smBreakPointReadsGR, 
#'                                          maxMismatch = 2,
#'                                          minOvrl = 10)
#' smBreakPointReadsGR <- as.data.table(mcols(smBreakPointReadsGR))
mergeOverlapNUMTs <- function(inGR, maxMismatch = 2, minOvrl = 10) {
  reducedGR <- inGR
  # this variable will save the number of merges done during one iteration 
  # of the cycle. The loop will stop then nothing is merged
  wasMerged <- -1
  while(wasMerged != 0) {
    wasMerged <- 0
    mergedBP <- c() # will save bps which were merged
    result <- GRanges()
    
    # order according to start to simplify
    reducedGR <- reducedGR[order(start(reducedGR))]
    # find overlaps between regions
    ovrl <- findOverlaps(reducedGR, ignore.strand = T, drop.self = T,
                         drop.redundant = T, minoverlap = minOvrl)
    ovrl <- as.data.frame(ovrl)
    # if no overlaps - nothing to do
    if (nrow(ovrl) == 0) {
      return(reducedGR)
    }
    
    # go though each overlap
    for (i in 1:nrow(ovrl)) {
      # we already know that 2 reads overlap, it means that
      # start of a left read is always to the left of a right one
      leftRead <- reducedGR[ovrl[i, 1]]
      rigthRead <- reducedGR[ovrl[i, 2]]
      
      # starts and ends of where to take substrings from reads from
      substBound <- rep(0, 4)
      names(substBound) <- c('leftStart', 'leftEnd', 'rigthStart', 'rigthEnd')
      substBound['leftStart'] <- start(rigthRead) - start(leftRead) + 1
      substBound['rigthStart'] <- 1
      # however, one of the reads maybe included in the other
      if (end(rigthRead) <= end(leftRead)) {
        ovrlStrLen <- nchar(rigthRead$mainRead)
        substBound['leftEnd'] <- end(rigthRead) - start(leftRead)
        substBound['rigthEnd'] <- nchar(rigthRead$mainRead)
      } else {
        ovrlStrLen <- end(leftRead) - start(rigthRead)
        substBound['leftEnd'] <- nchar(leftRead$mainRead)
        substBound['rigthEnd'] <- ovrlStrLen
      }
      
      # take and compare substrings
      leftOvrl <- substr(leftRead$mainRead, substBound['leftStart'],
                         substBound['leftEnd'])
      leftOvrl <- strsplit(leftOvrl, '')[[1]]
      rigthOvrl <- substr(rigthRead$mainRead, substBound['rigthStart'],
                          substBound['rigthEnd'])
      rigthOvrl <- strsplit(rigthOvrl, '')[[1]]
      if (length(leftOvrl) != length(rigthOvrl)) {
        stop("Lehgth of substrings do not match!")
      }
      
      # count number of mismatched between the substrings
      mismatch <- ovrlStrLen - sum(leftOvrl == rigthOvrl)
      
      if (mismatch <= maxMismatch) {
        # perform merging
        wasMerged <- wasMerged + 1
        mergedBP <- c(mergedBP, leftRead$bp, rigthRead$bp)
        newMainRead <- paste0(leftRead$mainRead, 
                              substr(rigthRead$mainRead, ovrlStrLen + 1, 
                                     nchar(rigthRead$mainRead)))
        toAdd <- data.frame(chr = 'chrM', start = start(leftRead), 
                            end = end(rigthRead), mainRead = newMainRead,
                            N = NA, fromSupp = NA, suppReads = NA, 
                            perc = NA, secondMainRead = NA, 
                            secondMainReadN = NA, secondMainReadPerc = NA,
                            bp = paste(leftRead$bp, rigthRead$bp, sep = '_'),
                            stringsAsFactors = F)
        result <- c(result, 
                    makeGRangesFromDataFrame(toAdd, keep.extra.columns = T))
      } else {
        result <- c(result, reducedGR[ovrl[i, 1]])
        result <- c(result, reducedGR[ovrl[i, 2]])
      }
    }
    # add reads which were not overlapped with anything
    notOverlapping <- setdiff(1:length(reducedGR), unique(unlist(ovrl)))
    result <- c(result, reducedGR[notOverlapping])
    reducedGR <- result
    reducedGR <- reducedGR[!duplicated(reducedGR)]
    reducedGR <- reducedGR[!reducedGR$bp %in% mergedBP]
  }
  result <- sort(reducedGR)
  result
}

#' writeAsFasta
#' Outputs most common reads into fasta format 
#' @param samDT data table with aligned reads
#' @param commonReadPerBreak data table with most common reads per break point,
#'                           result of getMostCommonReadPerBreakPoint
#' @param cutMatchTo number of bp to which matching part of the read should be
#'                   cut
#' @param putSto to which end concatinate soft clipped reads. In case of SM 
#'               it's start, and end in case of MS
#' @param filePath name of the file
#' @return fasta file
writeAsFasta <- function(bamObj, commonReadPerBreak, 
                         cutMatchTo = 20, putSto = 'start', filePath) {
  matchingBases <- sapply(commonReadPerBreak$bp,
                          function(x) bamObj$mBases[bamObj$breakPoint == x][1])
  if (putSto == 'start') {
    matchingBases <- sapply(matchingBases, substr, 1, cutMatchTo)
    commonReadPerBreak[, writeToBam := paste0(mainRead, matchingBases)]
  }
  if (putSto == 'end') {
    matchingBases <- reverse(sapply(reverse(matchingBases), substr, 1,
                                    cutMatchTo))
    commonReadPerBreak[, writeToBam := paste0(matchingBases, mainRead)]
  }
  
  file.remove(filePath)
  for (i in 1:nrow(commonReadPerBreak)) {
    write(paste0('>BP_', commonReadPerBreak$bp[i]), filePath, 1, T)
    write(commonReadPerBreak$writeToBam[i], filePath, 1, T)
  }
}

# Inputs ----------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument('--chrMfolder', nargs = 1, required = T, 
                    type = 'character', 
                    help = paste('Path to folder with BAM files.',
                                 'All BAM files in the folder will be',
                                 'analysed. Mapping have to be performed to',
                                 'chrM only. Otherwise memory overflow'))
parser$add_argument('--overhangBP', nargs = 1, required = T, default = 30,
                    type = 'integer', 
                    help = paste('Minimal length of soft clipped and matched',
                                 'parts of a read so that read would be',
                                 'valid for NUMTs extraction'))
parser$add_argument('--minimalBreakPointCoverage', nargs = 1, required = T, 
                    default = 20, type = 'integer', 
                    help = 'Minimal number of reads covering a breakpoint')
parser$add_argument('--mismatchSoftClipped', nargs = 1, required = T, 
                    default = 2, type = 'integer', 
                    help = paste('Maximum number of mismatches in MATCHED',
                                 '(MT) part of a read'))

args <- parser$parse_args()

message("Start time of run: ", Sys.time())
timeStart <- Sys.time()
printArgs(args)

# Read in mapped to chrM Illumina reads ---------------------------------------
# list all chrM bam files
allBamFiles <- list.files(args$chrMfolder, '.bam$', full.names = T)

# I intentionally do loop here and not lapply, because of the putative memory 
# issues
allBam <- list()
for (bamFile in allBamFiles) {
  message('[', Sys.time(), '] Started reading: ', bamFile)
  # read in bam
  oneBam <- scanBam(bamFile)[[1]]
  # add bam path
  oneBam$bamPath <- bamFile
  # select soft clipped reads which are long enough and have only soft clipped 
  # and matching bases
  oneBam <- selectSoftClipped(oneBam, args$overhangBP)
  if (!is.null(oneBam)) {
    allBam <- mergeBamObj(oneBam, allBam)
  }
  oneBam <- NULL
  message('[', Sys.time(), '] Finished reading: ', bamFile)
}

# Get break points for soft clipped Illumina reads ----------------------------
message('[', Sys.time(), '] Started calculating breakpoints')
# get break point for all reads
allBam <- getBreakPoint(allBam)
message('[', Sys.time(), '] Finished calculating breakpoints')

# separate on SM (beginning of NUMT) and MS (end of NUMT) reads
smReads <- allBam$readType == 'SM'
allBamSM <- lapply(allBam, function(x) x[smReads])
msReads <- allBam$readType == 'MS'
allBamMS <- lapply(allBam, function(x) x[msReads])

# plot histogram of number of reads having a break point per base pair
pdf('SM_MS_hist.pdf', width = 12, height = 8)
par(mfrow = c(2, 1))
hist(allBamSM$breakPoint, breaks = 14000, xaxt = "n",
     main = 'SM reads per basepair', ylab = 'Number of reads', 
     xlab = 'chrM, kBp', col = '#7B9971', border = '#7B9971')
axis(1, at = seq(0, 15000, 1000), labels = seq(0, 15, 1))
abline(args$minimalBreakPointCoverage, 0, lty = "dashed")
hist(allBamMS$breakPoint, breaks = 14000, xaxt = "n",
     main = 'MS reads per basepair', ylab = 'Number of reads', 
     xlab = 'chrM, kBp', col = '#CE6B5D', border = '#CE6B5D')
axis(1, at = seq(0, 15000, 1000), labels = seq(0, 15, 1))
abline(args$minimalBreakPointCoverage, 0, lty = "dashed")
dev.off()
par(mfrow = c(1, 1))
message('[', Sys.time(), '] Output histogram of number of reads having a ',
        'break point per base pair of chrM')

# number of reads having a break point per base pair, which have more than
# minimalBreakPointCoverage reads per break point
message('[', Sys.time(), '] Started selecting breakpoints based on coverage')
readsPerBreakPointSM <- getReadsPerBreakPoint(allBamSM,
                                              args$minimalBreakPointCoverage)
readsPerBreakPointMS <- getReadsPerBreakPoint(allBamMS,
                                              args$minimalBreakPointCoverage)
smBreakPointReads <- getMostCommonReadPerBreakPoint(allBamSM, 
                                                    readsPerBreakPointSM,
                                                    minNumbReads = args$minimalBreakPointCoverage, 
                                                    readLenCut = args$overhangBP, 
                                                    from = 'tail',
                                                    maxMismatch = args$mismatchSoftClipped)
msBreakPointReads <- getMostCommonReadPerBreakPoint(allBamMS, 
                                                    readsPerBreakPointMS, 
                                                    minNumbReads = args$minimalBreakPointCoverage, 
                                                    readLenCut = args$overhangBP, 
                                                    from = 'start',
                                                    maxMismatch = args$mismatchSoftClipped)
message('[', Sys.time(), '] Finished selecting breakpoints based on coverage')

pdf('barplotNconcensus.pdf', width = 12, height = 8)
barplotNconcensus(smBreakPointReads, readsPerBreakPointSM) +
  scale_fill_manual(values = c("#bdccb8", "#7b9971", "#3d4c38"))
barplotNconcensus(msBreakPointReads, readsPerBreakPointMS) +
  scale_fill_manual(values = c("#e1a69d", "#ce6b5d", "#67352e"))
message('[', Sys.time(), '] Plotted barplot of concensus reads')

# write to bam for visualization
smBreakPointReadsInd <- allBamSM$breakPoint %in% smBreakPointReads$bp
writeToBam(lapply(allBamSM, function(x) x[smBreakPointReadsInd]),
           'numtsStart')
msBreakPointReadsInd <- allBamMS$breakPoint %in% msBreakPointReads$bp
writeToBam(lapply(allBamMS, function(x) x[msBreakPointReadsInd]), 
           'numtsEnd')
message('[', Sys.time(), '] Wrote BAMs for visualization')

# Visualize break points ------------------------------------------------------
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(BSgenome.Dmelanogaster.UCSC.dm6)
alTrackStart <- AlignmentsTrack('numtsStart.srt.bam', isPaired = F, 
                                name = '# soft clipped reads', 
                                fill = '#7B9971',
                                fill.coverage = '#7B9971')
alTrackEnd <- AlignmentsTrack('numtsEnd.srt.bam', isPaired = F, 
                              name = 'End', fill = '#CE6B5D',
                              fill.coverage = '#CE6B5D')
pdf('SCreads_on_chrM.pdf', width = 12, height = 4)
plotTracks(c(OverlayTrack(trackList=list(alTrackStart, alTrackEnd)),
             sTrack, gtrack),
           chromosome = "chrM", from = 1, to = 14500, min.height = 1,
           type="coverage")
dev.off()
message('[', Sys.time(), '] Visualized breakpoints in SCreads_on_chrM.pdf')

# IMPORTANT NOTE: one shouldn't merge overlapping NUMTs, because:
# 1) some breakpoints can be DGRP-line specific
# 2) there will be no issue of multiple alignment for 454 reads because they 
#    are long and because I'll require almost full match afterwards

# Check similarity of the nuclear parts of SM - MS ----------------------------
compareMainReads(smBreakPointReads, args$overhangBP, 5)
compareMainReads(msBreakPointReads, args$overhangBP, 5)
message('[', Sys.time(), '] Checked similarity of the nuclear parts of SM-MS')

# Write SM and MS fastas to align 454 reads to them ---------------------------
writeAsFasta(allBamSM, smBreakPointReads, cutMatchTo = 30, putSto = 'start',
             'SM_overBreakPoint.fa')
writeAsFasta(allBamMS, msBreakPointReads, cutMatchTo = 30, putSto = 'end',
             'MS_overBreakPoint.fa')
message('[', Sys.time(), '] Written to fasta')