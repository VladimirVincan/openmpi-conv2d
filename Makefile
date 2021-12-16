.PHONY : all clean debug
all: serial.o parallel.o generate_matrices.o compare_matrices.o
serial.o: serial.cc matrix_operators.hpp mpi_operators.hpp
	g++ -o serial.o serial.cc
parallel.o: parallel.cc matrix_operators.hpp mpi_operators.hpp
	mpic++ -o parallel.o parallel.cc
generate_matrices.o: generate_matrices.cc matrix_operators.hpp
	g++ -o generate_matrices.o generate_matrices.cc
compare_matrices.o: compare_matrices.cc colors.hpp
	g++ -o compare_matrices.o compare_matrices.cc
debug: serial.o parallel.o
	g++ -g -o serial.o serial.cc
	mpic++ -g -o parallel.o parallel.cc
clean:
	rm serial.o parallel.o generate_matrices.o