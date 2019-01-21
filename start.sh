#!/bin/bash

gcc -o main main.c bloomfilter.c sha256.c -lm
./main
