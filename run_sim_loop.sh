#!/bin/bash
# Script to run the simulation
#
# Usage, locally:
#
#   ./run_sim.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_sim.sh
#
# Peregrine directives:
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=sim
#SBATCH --output=sim_%j.log

echo "seed: "$1
./simulation_logic_only --sim s$1 f$2

