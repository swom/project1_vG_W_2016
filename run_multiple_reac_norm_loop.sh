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
#SBATCH --job-name=reac_norm_m
#SBATCH --output=reac_norm_m.log

module load Qt5
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

change_freqs=(0, 1, 10)
for i in $(seq 1 40)
do
  for j in "${change_freqs[@]}"
do
  echo $i
  echo $j
sbatch run_reac_norm_loop.sh $i $j 
  done
done 