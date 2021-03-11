#!/bin/bash
# Script to run multiple simulations
#
# Usage, locally:
#
#   ./run_multiple_rand_evo_loop.sh
#
# Usage, on Peregrine:
#
#   sbatch ./run_multiple_rand_evo_loop.sh
#
# Peregrine directives:
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --job-name=run_multiple_rand_evo
#SBATCH --output=rand_evo_multiple.log

module load Qt5
export CC=g++
export CXX=g++
module load gompic/2019b
module load gompi
make clean
qmake simulation_logic_only.pro
make 


  for i in $(seq 1 30)
do
  for k in $(seq 0 49)
do
  echo $i
  echo $k
  sbatch run_rand_evo.sh $i $k
done
done

#move funders_success files directly in data partition
# &&
#copy sim_demographic files to data partition but keep also in home
#thus the program will check if it has run the same simulation twice

watch -n 1200 "mv rand_evo_a?*fund* ../../../data/p288427/ || true && cp rand_evo_a?*sim* ../../../data/p288427/"

