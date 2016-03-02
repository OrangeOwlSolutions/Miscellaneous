#!/bin/bash

mv timings.txt timings.old.txt

declare -i numExecutions=30

for M in 2 3 4 5 6 7 8
do 
    for ((K = 16; K <= 1048576; K = K * 2))
    do
    echo $K $M $numExecutions
    ./svdMKLDouble $M $M $K $numExecutions >> timings.txt
    done
done
