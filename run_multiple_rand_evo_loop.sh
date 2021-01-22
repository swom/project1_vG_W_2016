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
#SBATCH --time=1:00:00
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


  for i in $(seq 31 60)
do
  for k in $(seq 1 50)
do
  echo $i
  echo $k
  sbatch run_rand_evo.sh $i $k
done
done
