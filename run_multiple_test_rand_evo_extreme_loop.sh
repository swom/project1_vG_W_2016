#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_rand_evo_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_test_rand_evo_extreme_loop.sh
#
# Peregrine directives:
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=test_run_multiple_rand_evo
#SBATCH --output=rand_evo_multiple_test.log

module load Qt5
export CC=g++
export CXX=g++
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

for j in $(seq 3 5)
do
echo "seq : " $j 
  for i in $(seq 1 3)
do
  echo "seed : "$i
  sbatch run_test_rand_evo_extr.sh $i $j
done
done

