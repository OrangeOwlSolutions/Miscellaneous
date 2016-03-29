#!/bin/bash

icpc -lpthread -lmkl_sequential -lmkl_core -lmkl_intel_lp64 -lrt -openmp svdMKLSingle.cpp TimingCPU.cpp InputOutput.cpp -O3 -o svdMKLSingle

declare -i numExecutions=30

for M in 2 3 4 5 6 7 8
do 
    for ((K = 16 ; K <= 1048576; K = K * 2))
    do
    echo $K $M $numExecutions
    ./svdMKLSingle $M $M $K $numExecutions >> timings.txt
    done
done
