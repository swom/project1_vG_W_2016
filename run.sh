#!/bin/bash
# Script to run the simulation
#
# Usage, locally:
#
#   ./run.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run.sh
#
# Peregrine directives:
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run
#SBATCH --output=run_%j.log
module load Qt5
module load gompic/2019b
qmake simulation_logic_only.pro
make 
echo "seed: "$1
./simulation_logic_only s$1 

