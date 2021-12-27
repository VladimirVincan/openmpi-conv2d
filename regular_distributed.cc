#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

#include "matrix_operators.hpp"
#include "mpi_operators.hpp"

using namespace std;

void conv2(int const &ah, int const &aw, int const &bh, int const &bw, int &ch,
           int &cw, double const *const *const a, double const *const *const b,
           double **&c, int const &prank, int const &csize) {
  // initialize height & width
  cw = aw + bw - 1;
  ch = ah + bh - 1;

  // sendcounts, displs
  int *sendcounts = create_elementwise_sendcounts(csize, ch * cw);
  int *displs = create_displs(sendcounts, csize);

  // allocate memory
  if (!prank) alloc_matrix(ch, cw, c);
  double *c_local = new double[sendcounts[prank]]();

  // convolution
  for (int c_it = 0; c_it < sendcounts[prank]; ++c_it) {
    int ci = (displs[prank] + c_it) / cw;
    int cj = (displs[prank] + c_it) % cw;
    int aw_start = std::max(0, cj - bw + 1);
    // int bw_start = std::max(0, cj - aw + 1);
    int aw_end = std::min(aw, cj + 1);
    int bw_end = std::min(bw, cj + 1);
    int ah_start = std::max(0, ci - bh + 1);
    // int bh_start = std::max(0, ci - ah + 1);
    int ah_end = std::min(ah, ci + 1);
    int bh_end = std::min(bh, ci + 1);
    int wlen = aw_end - aw_start;  // same as bw_end - bw_start
    int hlen = ah_end - ah_start;  // same as bh_end - bh_start
    for (int i = 0; i < hlen; ++i)
      for (int j = 0; j < wlen; ++j)
        c_local[c_it] +=
            a[ah_start + i][aw_start + j] * b[bh_end - i - 1][bw_end - j - 1];
  }

  // gatherv
  if (!prank)
    MPI_Gatherv(c_local, sendcounts[prank], MPI_DOUBLE, c[0], sendcounts, displs,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(c_local, sendcounts[prank], MPI_DOUBLE, NULL, sendcounts,
                displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void read(char const *const file_name, int &ah, int &aw, int &bh, int &bw,
          double **&matrix, double **&kernel) {
  ifstream f(file_name);
  f >> ah >> aw >> bh >> bw;

  alloc_matrix(ah, aw, matrix);
  alloc_matrix(bh, bw, kernel);

  for (int i = 0; i < ah; ++i)
    for (int j = 0; j < aw; ++j) f >> matrix[i][j];

  for (int i = 0; i < bh; ++i)
    for (int j = 0; j < bw; ++j) f >> kernel[i][j];

  f.close();
}

void write(char const *const matrix_file_name, char const *const time_file_name,
           double const &mind, int const &ch, int const &cw,
           double const *const *const conv) {
  // print_matrix(ch, cw, conv);
  fprint_matrix(matrix_file_name, ch, cw, conv);

  // // cout << "Elapsed time distributed (" << csize
  // //      << " processes) = " << long(mind * 1000000000) << " [ns] " << endl;
  ofstream f(time_file_name, std::ios_base::app);
  f << long(mind * 1000000000) << std::endl;
  f.close();
}

int main(int argc, char *argv[]) {
  // init mpi
  int initialized, csize, prank;
  MPI_Initialized(&initialized);
  if (!initialized) MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &csize);
  MPI_Comm_rank(MPI_COMM_WORLD, &prank);

  double **matrix, **kernel, **conv;
  int ah, aw, bh, bw, ch, cw;

  // read sizes and matrices
  if (!prank) {
    char const file_name[] = "input/input.txt";
    read(file_name, ah, aw, bh, bw, matrix, kernel);
  }

  // start time measurement
  double start_time = MPI_Wtime();

  MPI_Bcast(&ah, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&aw, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bh, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&bw, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // alloc matrices
  if (prank) {
    alloc_matrix(ah, aw, matrix);
    alloc_matrix(bh, bw, kernel);
  }

  // distribute matrices
  MPI_Bcast(matrix[0], ah * aw, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(kernel[0], bh * bw, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  conv2(ah, aw, bh, bw, ch, cw, matrix, kernel, conv, prank, csize);

  // dealloc matrices
  if (prank) {
    dealloc_matrix(ah, matrix);
    dealloc_matrix(bh, kernel);
  }

  // end time measurement
  double end_time = MPI_Wtime();
  double total_time = end_time - start_time;
  double mind;
  MPI_Reduce(&total_time, &mind, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!prank) {
    char const matrix_file_name[] = "outputs/regular_distributed_matrix.txt";
    char const time_file_name[] = "outputs/regular_distributed_time.txt";
    write(matrix_file_name, time_file_name, mind, ch, cw, conv);

    dealloc_matrix(ah, matrix);
    dealloc_matrix(bh, kernel);
    dealloc_matrix(ch, conv);
  }

  // finalize mpi
  int finalized;
  MPI_Finalized(&finalized);
  if (!finalized) MPI_Finalize();

  return 0;
}
