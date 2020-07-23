#!/bin/bash

change_freqs=(0,1,10)
for i in $(seq 1 20)
do
  for j in "${change_freqs[@]}"
do
  echo $i
  echo $j
  sbatch run_sim.sh $i $j
  done
done 
