#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_rand_evo_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_rand_evo_loop.sh
#
# Peregrine directives:
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_rand_evo
#SBATCH --output=run_rand_evo.log

module load Qt5
export CC=g++
export CXX=g++
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

change_freqs=(0)
amplitude=(3)
for i in $(seq 1 100)
do
  for j in "${change_freqs[@]}"
do
	for z in "${amplitude[@]}"
do
		for k in $(seq 1 50)
do
  echo $i
  echo $j
  echo $z
  echo $k

  sbatch run_rand_evo.sh $i $j $z $k
		done
	done
  done
done 