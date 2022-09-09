#!/bin/bash

# DESCRIPTION -----------------------------------------------------------------
# This script is aimed at trimming, mapping, deduplicating and extracting
# reads mapping to chromosome M.
#
# ARGUMENTS -------------------------------------------------------------------
#    1) Tab separated table with three columns: SRR, DGRP, Instrument.
#    2) Path to reference genome, ALIGN TO CHROMOSOME M ONLY!
#    3) Path to folder with fastq-s
#    4) Path to the folder to store trimmed files (will be created)
#    5) Path to the folder to store mapped files (will be created)
#    6) Path to the folder to store deduplicated files (will be created)
#    7) Path to the folder to store chrM files (will be created)
#    8) Number of processes to run in parallel
# Example
#    REF_DIR=/home/litovche/Documents/RefGen/dm6_mitoOnly/
#    chmod +x 2_trimMapDedupl_SoftClip.sh
#    ./2_trimMapDedupl_SoftClip.sh Process_Illumina_Input.csv \
#                                  $REF_DIR'dm6.chrM.fa' \
#                                  fastq/Illumina/ trimmed/Illumina/ \
#                                  mapped/Illumina/ dedupl/Illumina/ \
#                                  chrMbam/Illumina/ 4
#

# SOFTWARE REQUIREMENTS -------------------------------------------------------
PATH_TO_PICARD=/home/litovche/bin/picard.jar

# FUNCTIONS -------------------------------------------------------------------

#===  FUNCTION  ===============================================================
# NAME:  waitall
# DESCRIPTION: waits till all processes finish to run new batch
#==============================================================================
# This function allows to run several processed in batches in the background
waitall() { # PID...
  ## Wait for children to exit and indicate whether all exited with 0 status.
  local errors=0
  while :; do
    #debug "Processes remaining: $*"
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
           set -- "$@" "$pid"
      elif wait "$pid"; then
        debug "$pid exited with zero exit status."
      else
        debug "$pid exited with non-zero exit status."
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    # TODO: how to interrupt this sleep when a child terminates?
    sleep ${WAITALL_DELAY:-1}
   done
  ((errors == 0))
}
debug() { echo "DEBUG: $*" >&2; }

#===  FUNCTION  ===============================================================
# NAME:  trimFastq
# DESCRIPTION: trims reads with trim-galore
# PARAMETER  1: path to folder to put trimmed files
# PARAMETER  2: path to fastq file - read 1
# PARAMETER  3: path to fastq file - read 2
#==============================================================================
trimFastq() {
   if (( $# == 2 )); then
      echo "Started trimming of "$2" at "$( date +%Y-%m-%d,%H-%M-%S )
      trim_galore -q 30 $2 -o $1 --fastqc
      echo "Finished trimming of "$2" at "$( date +%Y-%m-%d,%H-%M-%S )
   fi

   if (( $# == 3 )); then
      echo "Started trimming of "$2" and "$3" at "$( date +%Y-%m-%d,%H-%M-%S )
      trim_galore -q 30 --paired $2 $3 -o $1 --fastqc
      echo "Finished trimming of "$2" and "$3" at "$( date +%Y-%m-%d,%H-%M-%S )
   fi
}

#===  FUNCTION  ===============================================================
# NAME:  mapWithBWA
# DESCRIPTION: maps reads with BWA
# PARAMETER  1: path to reference genome
# PARAMETER  2: ID of the sequencing, i.e. DGRP_100
# PARAMETER  3: sample info, i.e. A1
# PARAMETER  4: library name, i.e lib1 or run1
# PARAMETER  5: folder to put mapped files after BWA
# PARAMETER  6: path to first read pair file
# PARAMETER  7: path to second read pair file (if pair end)
#==============================================================================
mapWithBWA() {
  echo "Started mapping of "$2" at "$( date +%Y-%m-%d,%H-%M-%S )
  # constract all names
  samName=$5$2"_"$3"_"$4".sam"
  bamName=$5$2"_"$3"_"$4".bam"
  bamSrtName=$5$2"_"$3"_"$4".srt.bam"
  bamSrtFlagstat=$5$2"_"$3"_"$4".srt.flagstat"

  # group info
  groupInfo="@RG\tID:"$2"\tSM:"$3"\tPL:illumina\tLB:"$4"\tPU:unit1"
  if (( $# == 7 )); then
    bwa mem -M -R $groupInfo $1 $6 $7 > $samName
  fi 
  if (( $# == 6 )); then
  	bwa mem -M -R $groupInfo $1 $6 > $samName
  fi   

  samtools view -F 0x4 -bS $samName > $bamName
  rm $samName
  samtools sort $bamName -o $bamSrtName
  rm $bamName
  samtools index $bamSrtName
  samtools flagstat $bamSrtName > $bamSrtFlagstat
 
 echo "Finished mapping of "$2" at "$( date +%Y-%m-%d,%H-%M-%S )
}

#===  FUNCTION  ===============================================================
# NAME:  deduplWithPicard
# DESCRIPTION: remove duplicates with picard
# PARAMETER  1: path to aligned file
# PARAMETER  2: folder to put deduplicated files
#==============================================================================
deduplWithPicard() {
  echo "Started deduplications of "$1" at "$( date +%Y-%m-%d,%H-%M-%S )
  outputFile=( $( echo $1 | sed 's@.*/@@g' | sed 's@.srt.bam@.dedupl.bam@g') )
  outputFile=$2$outputFile
  metrixFile=( $( echo $outputFile | sed 's@.dedupl.bam@.dedupl.metrix@g') )
  flagstatFile=( $( echo $outputFile | sed 's@.dedupl.bam@.dedupl.flagstat@g') )

  java -jar $PATH_TO_PICARD MarkDuplicates I=$1 \
          O=$outputFile \
          CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
          M=$metrixFile \
          REMOVE_DUPLICATES=TRUE
  samtools flagstat $outputFile > $flagstatFile

  echo "Finished deduplications of "$1" at "$( date +%Y-%m-%d,%H-%M-%S )
}

#===  FUNCTION  ===============================================================
# NAME:  getSoftCLippedOfChrM
# DESCRIPTION:
# PARAMETER  1: path to deduplicated file
# PARAMETER  2: folder to put chrM files
# PARAMETER  3: folder to put soft clipped reads
#==============================================================================
getReadsOnChrM() {
   echo "Started extraction of chrM of "$1" at "$( date +%Y-%m-%d,%H-%M-%S )

   echo "-----------------------------"
   echo $1 $2
 
   outputFile=( $( echo $1 | sed 's@.*/@@g' | sed 's@.bam@.chrM.bam@g') )
   outputFile=$2$outputFile
   outputSrtFile=( $( echo $1 | sed 's@.*/@@g' | sed 's@.bam@.chrM.srt.bam@g') )
   outputSrtFile=$2$outputSrtFile
   
   # extract chrM
   samtools view -F 0x4 -bS $1 chrM:1-14916 > $outputFile
   samtools sort  $outputFile -o $outputSrtFile
   samtools index $outputSrtFile
   rm $outputFile

   echo "Finished extraction of chrM of "$1" at "$( date +%Y-%m-%d,%H-%M-%S )
}

#===  FUNCTION  ===============================================================
# NAME:  trimMapDeduplExtrChrM
# DESCRIPTION: wrapper around functions to trim, map, deduplicate and extract
#              chrM reads
# PARAMETER  1: path to reference genome
# PARAMETER  2: SRR number
# PARAMETER  3: DGRP number
# PARAMETER  4: instrument name, i.e. Illumina
# PARAMETER  5: fastq directory
# PARAMETER  6: trimmed directory
# PARAMETER  7: mapped directory
# PARAMETER  8: deduplicated directory
# PARAMETER  9: chrM directory 
#==============================================================================
trimMapDeduplExtrChrM() {
   # get fastq files for srr
   #fastqs=( $(ls -l $5 | grep $2 | sed 's/.*SRR/SRR/g'))

   # trim
   #if [ ${#fastqs[@]} == 1 ]; then
   #   trimFastq $6 $5${fastqs[0]}
   #fi
   #if [ ${#fastqs[@]} == 2 ]; then
   #   trimFastq $6 $5${fastqs[0]} $5${fastqs[1]}
   #fi
   #trimmedFastq=( $(ls -l $6 | grep $currSRR | sed 's/.*SRR/SRR/g' | grep fq.gz))

   # map
   if [ ${#trimmedFastq[@]} == 1 ]; then
     mapWithBWA $1 $2 $3 $4 $7 $6${trimmedFastq[0]}
   fi
   if [ ${#trimmedFastq[@]} == 2 ]; then
      mapWithBWA $1 $2 $3 $4 $7 $6${trimmedFastq[0]} $6${trimmedFastq[1]}
   fi

   # deduplicate
   mappedBam=( $(ls -l $7 | grep $2 | sed 's/.*SRR/SRR/g' | grep .srt.bam$))
   deduplWithPicard $7$mappedBam $8

   # cromosome M and soft clipped
   deduplBam=( $(ls -l $8 | grep $2 | sed 's/.*SRR/SRR/g' | grep .dedupl.bam$))
   getReadsOnChrM $8$deduplBam $9
}

# READ INPUTS -----------------------------------------------------------------
inTable=$1
refGen=$2
fastqDir=$3
trimmedDir=$4
bwaDir=$5
deduplDir=$6
chrMdir=$7
numbOfProc=$8

mkdir $trimmedDir $bwaDir $deduplDir $chrMdir

SRRs=('zero')
DGRPs=('zero')
Intrument=('zero')

echo $inTable

# read in the input table
while read col1 col2 col3
do
  SRRs+=($col1)
  DGRPs+=($col2)
  Intrument+=($col3)
done < $inTable

# PROCESSING ------------------------------------------------------------------
srrCount=${#SRRs[@]}

for (( i=1; i<$srrCount; i+=$numbOfProc )); do
   # bash arrays are 0-indexed
   pids=""
   endInd=$(($numbOfProc - 1))

   for t in $(seq 0 $endInd); do
      position=$(( $i + $t ))

      # Get srr, dgrp and instrument 
      currSRR=${SRRs[$position]}
      currDGRP=${DGRPs[$position]}
      currInstrument=${Intrument[$position]}

      echo 'Input arguments:'
      echo $refGen 
      echo $currSRR
      echo $currDGRP
      echo $currInstrument
      echo $fastqDir
      echo $trimmedDir
      echo $bwaDir
      echo $deduplDir
      echo $chrMdir

      trimMapDeduplExtrChrM $refGen $currSRR $currDGRP $currInstrument \
                            $fastqDir $trimmedDir $bwaDir $deduplDir $chrMdir \
                            1>$currSRR".trimToSc.out" \
                            2>$currSRR".trimToSc.err" &
      pids="$pids $!"
   done
   waitall $pids
done