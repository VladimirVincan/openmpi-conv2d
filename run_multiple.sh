#!/bin/bash

echo "3 3 3 3 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_3x3_3x3/
mv outputs/* results/mat_3x3_3x3/

echo "10 10 10 10 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_10x10_10x10/
mv outputs/* results/mat_10x10_10x10/

echo "100 100 100 100 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_100x100_100x100/
mv outputs/* results/mat_100x100_100x100/

echo "1000 1000 1000 1000 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_1000x1000_1000x1000/
mv outputs/* results/mat_1000x1000_1000x1000/
