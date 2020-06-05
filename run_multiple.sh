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
#SBATCH --output=run.log


# simulation_logic_only has this command-line interface:
#
# simulation_logic_only [seed] [output filename]

for i in $(seq 1 2)
do
  echo $i
  sbatch run.sh $i 
done 
