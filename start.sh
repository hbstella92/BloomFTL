#!/bin/bash

rm simulationFTL
gcc -g -fsanitize=address -o simulationFTL main.c bloomfilter.c sha256.c zpipe.c -lm -lz && ./simulationFTL
