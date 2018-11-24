#!/bin/bash

#------------------------------------------------------------------------------
# STEP 9: CHECK THE QUALITY OF THE MAPING ON DM6
#------------------------------------------------------------------------------

echo "Name NR_total NR_mapped NR_total_dedupl NR_mapped_dedupl NR_2L NR_2R NR_3L NR_3R NR_4 NR_mtDNA NR_X NR_Wolb chrMminCov chrMmedianCov chrMmaxCov" > align_to_dm6_stats.txt
allRuns=($(ls -d -1 6_mapped_dm6/* | grep align.srt.bam$ | sed 's@.*/@@g' | sed 's@_.*@@g'))

for ((i=0; i<=${#allRuns[@]}; i++)); do
    sample=${allRuns[i]}
    toPrint=($sample)
    toGrep=`echo $sample | cut -f1,2 -d'_'`
    toPrint+=(`head -200 6_mapped_dm6/$toGrep*_BWA_dm6.align.srt.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 6_mapped_dm6/$toGrep*_BWA_dm6.align.srt.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`samtools idxstats 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.bam | grep 'chr2L\|chr2R\|chr3L\|chr3R\|chr4\|chrX\|chrM\|Wolbachia' | grep -v random | cut -f3 | xargs`)
       
    # get min/max and median coverage per chrM
    toPrint+=(`samtools depth 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.bam -r chrM | head -14915 | cut -f3 | sort -nk 1 | head -1`)
    toPrint+=(`samtools depth 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.bam -r chrM | head -14915 | cut -f3 | sort -nk 1 | sed -n '7458p'`)
    toPrint+=(`samtools depth 7_dedupl_dm6/$toGrep*_BWA_dm6.align.srt.dedupl.bam -r chrM | head -14915 | cut -f3 | sort -nk 1 | tail -1`)
    
    printf '%s\n' "${toPrint[@]}" | xargs >> align_to_dm6_stats.txt

done


