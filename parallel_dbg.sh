# https://stackoverflow.com/questions/329259/how-do-i-debug-an-mpi-program
mpiexec -n 2 xterm -hold -e gdb ./parallel.o -ex run
