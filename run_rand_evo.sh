#!/bin/bash
# Script to run the resulting population
# of a simualtion run against a set of random conditions
#
# Usage, locally:
#
#   ./run_rand.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_rand.sh
#
# Peregrine directives:
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=rand_evo_%j
#SBATCH --output=rand_evo_%j.log

echo "original seed" $1
echo "sequence_index" $2
./simulation_logic_only --rand_evo s$1 f0 a3 n50 cs50 si$2

