#!/usr/bin/bash

#$ -P bf528
#$ -cwd
#$ -j y
#$ -pe mpi_16_tasks_per_node 16

echo "Running job $JOB_ID"
echo "Started: $(date +%F)"
echo "Running in directory: $PWD"

#set environment
module load anaconda3
source activate py38

#full path to call executable
/projectnb2/bf528/students/rkafrawi/project_3/Run_data/SRR3879604/data_curator_efficient_RK3.py

echo "Job finished: $(date +%F)"
