#!/bin/bash

mkdir outputs/
mkdir input/

echo "4 4 4 4 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_4x4_4x4/
mv outputs/* results/mat_4x4_4x4/

echo "8 8 8 8 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_8x8_8x8/
mv outputs/* results/mat_8x8_8x8/

echo "16 16 16 16 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_16x16_16x16/
mv outputs/* results/mat_16x16_16x16/

echo "32 32 32 32 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_32x32_32x32/
mv outputs/* results/mat_32x32_32x32/

echo "64 64 64 64 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_64x64_64x64/
mv outputs/* results/mat_64x64_64x64/

echo "128 128 128 128 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_128x128_128x128/
mv outputs/* results/mat_128x128_128x128/

echo "256 256 256 256 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_256x256_256x256/
mv outputs/* results/mat_256x256_256x256/

echo "512 512 512 512 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_512x512_512x512/
mv outputs/* results/mat_512x512_512x512/

echo "1024 1024 1024 1024 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_1024x1024_1024x1024/
mv outputs/* results/mat_1024x1024_1024x1024/

echo "2048 2048 2048 2048 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_2048x2048_2048x2048/
mv outputs/* results/mat_2048x2048_2048x2048/

echo "4096 4096 4096 4096 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_4096x4096_4096x4096/
mv outputs/* results/mat_4096x4096_4096x4096/

echo "8192 8192 8192 8192 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_8192x8192_8192x8192/
mv outputs/* results/mat_8192x8192_8192x8192/

echo "16384 16384 16384 16384 generic" > mat_size.txt
./run_single.sh
rm  outputs/*matrix.txt
mkdir results/mat_16384x16384_16384x16384/
mv outputs/* results/mat_16384x16384_16384x16384/

rm -r outputs/
rm -r input/
