#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run
#SBATCH --output=run_sim.log


# simulation_logic_only has this command-line interface:
#
# simulation_logic_only s[seed] f[change_frequency]
module load Qt5
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

change_freqs=(0, 1, 10)
for i in $(seq 1 20)
do
  for j in "${change_freqs[@]}"
do
  echo $i
  echo $j
  sbatch run_sim_loop.sh $i $j
  done
done 