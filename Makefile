.PHONY : all clean debug
all: serial.o parallel.o generate_matrices.o compare_matrices.o regular_serial.o regular_distributed.o
serial.o: serial.cc matrix_operators.hpp mpi_operators.hpp
	g++ -o serial.o serial.cc
parallel.o: parallel.cc matrix_operators.hpp mpi_operators.hpp
	mpic++ -o parallel.o parallel.cc
regular_serial.o: regular_serial.cc matrix_operators.hpp
	g++ -o regular_serial.o regular_serial.cc
regular_distributed.o: regular_distributed.cc matrix_operators.hpp
	mpic++ -o regular_distributed.o regular_distributed.cc
generate_matrices.o: generate_matrices.cc matrix_operators.hpp
	g++ -o generate_matrices.o generate_matrices.cc
compare_matrices.o: compare_matrices.cc colors.hpp
	g++ -o compare_matrices.o compare_matrices.cc
debug: serial.o parallel.o regular_serial.o regular_distributed.o generate_matrices.o
	g++ -g -o serial.o serial.cc
	g++ -g -o regular_serial.o regular_serial.cc
	mpic++ -g -o parallel.o parallel.cc
	mpic++ -g -o regular_distributed.o regular_distributed.cc
clean:
	rm serial.o parallel.o regular_serial.o regular_distributed.o generate_matrices.o
