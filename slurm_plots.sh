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

#make plots

samples=1
iterations=400

#python3 stats.py ./build/sdcsim ./sdcasm/rule110.sdcasm --strands 101  --strands_range 1 --strands_step 1 -t 10000 --time_range 20000 --time_step 50 -T 25 -r 1 -s $samples -i $iterations -c dp100.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 21   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp20t.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp40t.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 11   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp10t.csv

python3 stats.py ./build/sdcsim ./sdcasm/rule28_000000.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp000000.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_111111.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp111111.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_101010.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp101010.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_110010.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp110010.csv

python3 stats.py ./build/sdcsim ./sdcasm/rule28_2.sdcasm  --strands 101   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp2.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_4.sdcasm  --strands 101   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp4.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_8.sdcasm  --strands 101   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp8.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_16.sdcasm --strands 101   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp16.csv
python3 stats.py ./build/sdcsim ./sdcasm/rule28_32.sdcasm --strands 101   --strands_range 1 --strands_step 1 -t 99999 --time_range 1 --time_step 1 -T 50 -r 100 -s $samples -i $iterations -c dp32.csv


#python3 stats.py ./build/sdcsim ./sdcasm/rule110.sdcasm --strands 101  --strands_range 1 --strands_step 1 -t 10000 --time_range 20000 --time_step 50 -T 25 -r 1 -s $samples -i $iterations -c dp100.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 21   --strands_range 1 --strands_step 1 -t 99999 --time_range 199998 --time_step 8000 -T 23 -r 1 -s $samples -i $iterations -c dp20.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 41   --strands_range 1 --strands_step 1 -t 99999 --time_range 199998 --time_step 8000 -T 23 -r 1 -s $samples -i $iterations -c dp40.csv
#python3 stats.py ./build/sdcsim ./sdcasm/rule28.sdcasm --strands 11   --strands_range 1 --strands_step 1 -t 99999 --time_range 199998 --time_step 8000 -T 23 -r 1 -s $samples -i $iterations -c dp10.csv

