#!/bin/bash

rm main
gcc -o main main.c bloomfilter.c sha256.c -lm
./main
