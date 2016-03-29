#!/bin/bash

#for M in 2 3 4 5 6 7 8
for M in 4 5 6 7 8 9 10
do 
    declare -i Nrows__=$M
    declare -i Ncols__=$M-2
    for ((K = 16; K <= 1048576; K = K * 2))
    do
       echo $Nrows__ $Ncols__ $K 

       declare -i numMatrices__=$K

      /usr/local/cuda-7.5/bin/nvcc -O3 -gencode arch=compute_35,code=sm_35  -odir "src" -M -D Nrows_=$Nrows__ -D Ncols_=$Ncols__ -D numMatrices_=$numMatrices__ -o "src/kernel.d" "../src/kernel.cu"

       /usr/local/cuda-7.5/bin/nvcc -O3 --compile --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35  -x cu -D Nrows_=$Nrows__ -D Ncols_=$Ncols__ -D numMatrices_=$numMatrices__ -o "src/kernel.o" "../src/kernel.cu"

       /usr/local/cuda-7.5/bin/nvcc --cudart static --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35 -link -o  "SVD_v2"  ./src/InputOutput.o ./src/TimingCPU.o ./src/TimingGPU.o ./src/Utilities.o ./src/kernel.o   -lcublas

       ./SVD_v2
    done
done

