#!/bin/bash

#------------------------------------------------------------------------------
# STEP 4: CHECK THE QUALITY OF THE MAPING ON DM3
#------------------------------------------------------------------------------

Run1=($(ls -d -1 0_rawData/Run1/* | grep _R1_ | sed 's@.*/@@g' | sed 's/_R1.*//g'))
Run2=($(ls -d -1 0_rawData/Run2/* | grep _R1_ | sed 's@.*/@@g' | sed 's/_R1.*//g'))
Run3=($(ls -d -1 0_rawData/Run3/* | grep _R1_ | sed 's@.*/@@g' | sed 's/_R1.*//g'))

runNumbs=()
for i in "${Run1[@]}"; do
    runNumbs+=("1")
done

for i in "${Run2[@]}"; do
    runNumbs+=("2")
done

for i in "${Run3[@]}"; do
    runNumbs+=("3")
done


allRuns=("${Run1[@]}" "${Run2[@]}" "${Run3[@]}")
echo "Total number of samples" ${#allRuns[@]}

echo "Name Run_Numb NR_total NR_mapped NR_total_dedupl NR_mapped_dedupl NR_2L NR_2R NR_3L NR_3R NR_4 NR_U NR_X NR_mtDNA NR_Wolb" > align_to_dm3_stats.txt
for ((i=0; i<=${#allRuns[@]}; i++)); do
    sample=${allRuns[i]}
    toPrint=($sample)
    toPrint+=(${runNumbs[i]})
    toGrep=`echo $sample | cut -f1,2 -d'_'`

    toPrint+=(`head -200 2_mapped_dm3/$toGrep*_dm3.align.srt.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 2_mapped_dm3/$toGrep*_dm3.align.srt.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 3_dedupl_dm3/$toGrep*_dm3.align.srt.dedupl.flagstat | grep "in total" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`head -200 3_dedupl_dm3/$toGrep*_dm3.align.srt.dedupl.flagstat | grep "mapped (" | sed 's/.*t://' | sed 's/ .*//'`)
    toPrint+=(`samtools idxstats 3_dedupl_dm3/$toGrep*_dm3.align.srt.dedupl.bam | cut -f 1,3 | grep -v Het | grep -v extra | cut -f2 | xargs | cut -d' ' -f1-9`)
    printf '%s\n' "${toPrint[@]}" | xargs >> align_to_dm3_stats.txt
    
done

