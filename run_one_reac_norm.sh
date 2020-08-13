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
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=reac_norm
#SBATCH --output=reac_norm%j.log

module load Qt5
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 

echo "seed: 1"
echo "freq: 0"
sbatch ./simulation_logic_only --reac_norm s1 f0


