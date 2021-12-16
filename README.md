# Run tests #
Run *./run_multiple.sh* in terminal. The results will be written in the *results/* folder. 

## Options ##
Inside *run_multiple.sh* determine matrix sizes.

Inside *run_single.sh* determine number of iterations and number of threads.

# Analyze results #
Run *python3 analyze_results.py* in terminal. Images will be shown in *results/* folder.

# Additional info #
## Matrix size ##
Determine matrix size inside *mat_size.txt*:
ah aw bh bw
generic|sobelx|sobely

Run *./generate_matrices.o* in terminal. The generated random matrices will appear in *input.txt*.
