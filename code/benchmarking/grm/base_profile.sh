#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:00:00
#$ -l mem=1G
#$ -t 1-27
#$ -o logs/grm_benchmarking.out
#$ -e logs/grm_benchmarking.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.11
source tsenv/bin/activate

echo "Job started at: `date`"

python3 scripts/benchmarking/grm/01_base_profile.py $SGE_TASK_ID
#-censor

echo "Job finished with exit code $? at: `date`"
