#!/bin/bash -l
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=12:00:00
#$ -l mem=4G
#$ -t 2-100
#$ -o logs/balsac_grm.out
#$ -e logs/balsac_grm.err
#$ -M b.lehmann@ucl.ac.uk
#$ -m beas

# load software and environments
module load python/3.9
source tsrelatedness/bin/activate
#pip install --no-index --upgrade pip
#pip install --no-index -r $HOME/projects/ctb-sgravel/python_environments/pedsim_requirements.txt

dir=/home/ucakble/Projects/genome_simulations
ts_path=$dir/tree_sequences
pedigree_name=$dir/data/balsac_pedigree.csv

echo "Job started at: `date`"

python3 code/compute_grm.py \
-d $dir \
-o $ts_path \
-p $pedigree_name \
-chr 3 \
-rep $SGE_TASK_ID
#-censor

echo "Job finished with exit code $? at: `date`"
