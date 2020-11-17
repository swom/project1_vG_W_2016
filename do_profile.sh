#!/bin/bash
module load Qt5
make clean
export CXX=g++; export CC=gcc; qmake DEFINES+=IS_ON_TRAVIS DEFINES+=LOGIC_ONLY simulation_logic_only.pro
make release
./simulation_logic_only --profile
gprof simulation_logic_only > gprof.log
head gprof.log -n 300

