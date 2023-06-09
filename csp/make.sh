#!/bin/bash

rm -r build
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
# cmake -DCMAKE_BUILD_TYPE=Release ..
make

# ./serial ../out_full.txt --nogen
# ./serial ../out_100_100_full.txt --nogen
# Run this 100 times
./openmp ../graphs/full_100.txt --nogen
# for i in {1..1000}; do
#     ./openmp ../out_100_100_full.txt --nogen
# done

# cd /global/homes/a/aarora/CS267/hw2-1/build
# ./serial -n 1000 -s 1 -o correct.parts.out
# ./openmp -n 1000 -s 1 -o openmp.parts.out

# module load python/3.9-anaconda-2021.11
# ../hw2-correctness/correctness-check.py openmp.parts.out correct.parts.out

# salloc -N 1 -C cpu -q interactive -t 01:00:00
# ./serial
# ./serial -n 10000
