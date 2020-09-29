#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_rand_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_rand_loop.sh
#
# Peregrine directives:
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run
#SBATCH --output=run_rand.log


# simulation_logic_only has this command-line interface:
#
# simulation_logic_only --rand s[seed] f[change frequency] a[amplitude of change] n[number of random conditions]
module load Qt5
export CC=g++
export CXX=g++
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

change_freqs=(0, 1, 10)
amplitude=(1.5, 2, 2.5, 3)
for i in $(seq 1 40)
do
  for j in "${change_freqs[@]}"
do
  for z in "${amplitude[@]}"
do
  echo $i
  echo $j
  echo $z
  sbatch run_rand_loop.sh $i $j $z
	done
  done
done 