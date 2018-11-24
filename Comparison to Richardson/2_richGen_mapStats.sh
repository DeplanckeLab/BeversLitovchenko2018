#!/bin/bash

#------------------------------------------------------------------------------
# STEP 2: CHECK THE QUALITY OF THE MAPING ON RICHADSON GENOME
#------------------------------------------------------------------------------

echo "Name NR_total NR_mapped NR_total_dedupl NR_mapped_dedupl NR_U NR_3L NR_Wolb" > align_to_rich_stats.txt
allRuns=($(ls -d -1 2_mapped_richardson/* | grep align.srt.bam$ | sed 's@.*/@@g' | sed 's@_richardson.align.srt.bam@@g'))

for ((i=0; i<=${#allRuns[@]}; i++)); do
    sample=${allRuns[i]}
    toPrint=($sample)
    toGrep=$toPrint
    toPrint+=(`head -200 2_mapped_richardson/$toGrep"_"*richardson.align.srt.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 2_mapped_richardson/$toGrep*_richardson.align.srt.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 3_dedupl_richardson/$toGrep*_richardson.align.srt.dedupl.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 3_dedupl_richardson/$toGrep*_richardson.align.srt.dedupl.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`samtools idxstats 3_dedupl_richardson/$toGrep*_richardson.align.srt.dedupl.bam | grep 'chrU\|chr3L\|Wolbachia' | grep -v random | cut -f3 | xargs`)

    printf '%s\n' "${toPrint[@]}" | xargs >> align_to_RichGen_stats.txt
    echo "Finished " $toPrint
done

exit 0;
