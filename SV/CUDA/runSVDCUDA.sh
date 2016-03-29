#!/bin/bash

declare -i Nrows_=6

/usr/local/cuda-7.5/bin/nvcc -O3 --compile --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35  -x cu -DNrows_ -o "src/kernel.o" "../src/kernel.cu"

#declare -i numExecutions=30

#for M in 2 3 4 5 6 7 8
#do 
#    for ((K = 16; K <= 1048576; K = K * 2))
#    do
#    echo $K $M $numExecutions
#    ./svdMKLDouble $M $M $K $numExecutions >> timings.txt
#    done
#done
