#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_switches.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_switches.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_multiple_switches
#SBATCH --output=run_multiple_switches.log


# simulation_logic_only has this command-linae interface:
#
# simulation_logic_only --[function] s[seed] f[change_frequency]
module load Qt5
export CC=g++
export CXX=g++
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

change_freqs=(0,1,0.5)

for i in (seq 1 100)
do
  for j in "${change_freqs[@]}"
do
  echo $i
  echo $j
  sbatch run_sim_loop.sh $i $j
  done
done 
