#!/bin/bash
rm outputs/*
iterations=4
threads=6
for ((i = 0 ; i < $iterations ; i++)); do
    # echo "iteration = $i"
    ./generate_matrices.o
    ./serial.o
    python3 conv.py
    for ((j = 2 ; j <= $threads ; j++)); do
        ./parallel.sh $j
        cat outputs/distributed_time.txt >> "outputs/distributed_time_it${iterations}_th${j}.txt"
        rm outputs/distributed_time.txt
    done
    ./compare_matrices.o
    rc=$?
    if [ $rc -ne 0 ]
    then
        break
    fi
done
mv outputs/serial_time.txt "outputs/serial_time_it${iterations}.txt"
