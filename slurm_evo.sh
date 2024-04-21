#!/bin/bash
#
# filename: slurm.sh
#
#
#SBATCH -J sdcsim
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -p long
#

echo 'Your job is running on node(s):'
echo $SLURM_JOB_NODELIST
echo 'Cores per node:'
echo $SLURM_TASKS_PER_NODE

# load module enviroment
. modules.sh

for pmut in 1 2 4 8 16 32 64 100
do
    for mgenes in 1 2 3 4
    do
        for drand in 1 5 10 100 500
        do
            echo "python evo.py --pop 8 --pmut $pmut --mgenes $mgenes --gen 300 --drand $drand --profile"
            python evo.py --pop 8 --pmut $pmut --mgenes $mgenes --gen 300 --drand $drand --profile --tour 2
        done
    done
done
