#!/bin/bash

change_freqs=(0, 1, 10)
amplitude=(1.5, 2, 2.5, 3)
for i in $(seq 1 20)
echo 1
do
  for j in "${change_freqs[@]}"
do
  for z in "${amplitude[@]}"
do
  echo $i
  echo $j
  echo $z
	done
  done
done 