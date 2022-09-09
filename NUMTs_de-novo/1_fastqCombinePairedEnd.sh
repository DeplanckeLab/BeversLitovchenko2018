#!/bin/bash
#SBATCH --job-name="fix"
#SBATCH --error="fix_%A_%a.err"
#SBATCH --output="fix_%A_%a.out"
#SBATCH --array=0-8 # Job Array
#SBATCH --nodes 1 # We can only ask for one node for free accounts
#SBATCH --ntasks=1 # Total number of MPI tasks
#SBATCH --time=48:00:00 # 6 hour time limit for free accounts
#SBATCH --ntasks-per-node=1 # Number of MPI task per compute node (<=8)
#SBATCH --cpus-per-task=1 # Number of OpenMP threads
#SBATCH --workdir /scratch/litovche/SRP000694-err/
#SBATCH --mem-per-cpu=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maria.litovchenko@epfl.ch

# DESCRIPTION -----------------------------------------------------------------
# Shell script wrapper for 1_fastqCombinePairedEnd.py to run on HPC

# INPUTS ----------------------------------------------------------------------
fastqDir=fastq/ # path to directory with fastqs
# a path to tab-separated table with 2 columns: first - SRR ID of a sample and
# second is instrument on which that sample was sequenced. File paths should 
# follow the pattern: fastqDir/intrument/SRR_1.fastq.gz
inTable='inventory.txt'

# READING IN ------------------------------------------------------------------
SRRs=('zero')
Instrument=('zero')

# read in the input table
while read col1 col2
do
    SRRs+=($col1)
    Instrument+=($col2)
done < $inTable

currSRR=${SRRs[${SLURM_ARRAY_TASK_ID}]}
intrumentDir=${Instrument[${SLURM_ARRAY_TASK_ID}]}

# SYNC FASTQs -----------------------------------------------------------------
python 1_fastqCombinePairedEnd.py\
       $fastqDir$intrumentDir"/"$currSRR"_1.fastq.gz" \
       $fastqDir$intrumentDir"/"$currSRR"_2.fastq.gz"
rm $fastqDir$intrumentDir"/"$currSRR"_1.fastq.gz" \
   $fastqDir$intrumentDir"/"$currSRR"_1.fastq.gz"