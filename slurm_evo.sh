#!/bin/bash
#
# filename: slurm.sh
#
#
#SBATCH -J sdcsim
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -p long
#

echo 'Your job is running on node(s):'
echo $SLURM_JOB_NODELIST
echo 'Cores per node:'
echo $SLURM_TASKS_PER_NODE

# load module enviroment
. modules.sh

# load python
ml Python
ml Makefile

python evo.py
