#!/bin/bash

W_TYPE=$1
R_TYPE=$2

rm simulationFTL
gcc -g -fsanitize=address -o simulationFTL main.c bloomfilter.c sha256.c -lm && ./simulationFTL $W_TYPE $R_TYPE
