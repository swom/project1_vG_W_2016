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
#SBATCH --job-name=rand_best
#SBATCH --output=rand_best_%j.log


# simulation_logic_only has this command-line interface:
#
# simulation_logic_only --rand_best s[seed] f[change frequency] a[amplitude of change] n[number of random conditions]
echo "seed: "$1
echo "freq: "$2
echo "ampl: "$3
./simulation_logic_only --rand_best s$1 f$2 a$3

