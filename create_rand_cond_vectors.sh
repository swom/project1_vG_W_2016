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
#SBATCH --job-name= create_rand_cond
#SBATCH --output= create_rand_cond.log


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

amplitude=(1.5, 2, 2.5, 3) 

 for i in "${amplitude[@]}"
do
  echo $i
  ./simulation_logic_only --create_rand_cond_vec a$i