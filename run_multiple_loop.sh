#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_loop.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=test
#SBATCH --output=test.log

change_freqs=(0, 1, 10)
amplitude=(1.5, 2, 2.5, 3)
for i in $(seq 1 20)
do
  for j in "${change_freqs[@]}"
do
  for z in "${amplitude[@]}"
do
  echo $i
  echo $j
  echo $z
	done
  done
done 