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
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=rand_evo
#SBATCH --output=rand_evo_%j.log

echo "seed: "$1
echo "freq: "$2
echo "ampl: "$3
echo "condition: "$4
./simulation_logic_only --rand_evo s$1 f$2 a$3 n50 sr0 nr$4

