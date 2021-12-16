#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>

// #include "colors.hpp"
#include "matrix_operators.hpp"
#include "mpi_operators.hpp"

using namespace std;

void initTwiddle(double *const &wCOS, double *const &wSIN, int const &length,
                 int const &prank, int const &csize) {
  int *sendcounts = create_elementwise_sendcounts(csize, length / 2);
  int *displs = create_displs(sendcounts, csize);
  double *wCOS_v = new double[sendcounts[prank]];
  double *wSIN_v = new double[sendcounts[prank]];

  for (int i = displs[prank]; i < displs[prank] + sendcounts[prank]; ++i) {
    wCOS_v[i - displs[prank]] = cos(-2.0 * PI * i / length);
    wSIN_v[i - displs[prank]] = sin(-2.0 * PI * i / length);
  }

  MPI_Allgatherv(wCOS_v, sendcounts[prank], MPI_DOUBLE, wCOS, sendcounts,
                 displs, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(wSIN_v, sendcounts[prank], MPI_DOUBLE, wSIN, sendcounts,
                 displs, MPI_DOUBLE, MPI_COMM_WORLD);

  delete[] wCOS_v;
  delete[] wSIN_v;
  delete[] sendcounts;
  delete[] displs;
}

void butterfly(double const topRE_i, double const topIM_i,
               double const bottomRE_i, double const bottomIM_i,
               double *const topRE_o, double *const topIM_o,
               double *const bottomRE_o, double *const bottomIM_o, int const &k,
               double const *const wCOS, double const *const wSIN) {
  double botRE = bottomRE_i * wCOS[k] - bottomIM_i * wSIN[k];
  double botIM = bottomRE_i * wSIN[k] + bottomIM_i * wCOS[k];

  *topRE_o = topRE_i + botRE;
  *topIM_o = topIM_i + botIM;
  *bottomRE_o = topRE_i - botRE;
  *bottomIM_o = topIM_i - botIM;
}

void fft(int const &size, int const &log2_size, double *const dataRE,
         double *const dataIM, double const *const wCOS,
         double const *const wSIN) {
  for (int i = 0; i < log2_size; ++i) {
    int m = 1 << (i + 1);
    int m2 = 1 << i;  // "half of m"
    for (int j = 0; j < m2; ++j) {
      for (int k = j; k < size; k += m) {
        butterfly(dataRE[k], dataIM[k], dataRE[k + m2], dataIM[k + m2],
                  dataRE + k, dataIM + k, dataRE + k + m2, dataIM + k + m2,
                  j << (log2_size - i - 1), wCOS, wSIN);
      }
    }
  }
}

int bit_reversal(int val, int const &log2_val) {
  int reversed = 0;
  for (int i = 0; i < log2_val; ++i) {
    reversed <<= 1;
    reversed |= (val & 1);
    val >>= 1;
  }
  return reversed;
}

void fft2_a_and_b(int const &height, int const &log2h, int const &width,
                  int const &log2w, double *const *const aRE,
                  double *const *const aIM, double *const *const bRE,
                  double *const *const bIM, double *const *const cRE,
                  double *const *const cIM, int const &prank,
                  int const &csize) {
  // allocate direct and transposed sendcounts and displs
  int *sendcounts_D = create_sendcounts(csize, height, width);
  int *sendcounts_T = create_sendcounts(csize, width, height);
  int *displs_D = NULL, *displs_T = NULL;
  if (!prank) {
    displs_D = create_displs(sendcounts_D, csize);
    displs_T = create_displs(sendcounts_T, csize);
  }

  // allocate temporary direct and transposed matrices A
  double **aRE_tD, **aIM_tD, **aRE_tT, **aIM_tT;
  alloc_matrix(sendcounts_D[prank] / width, width, aRE_tD);
  alloc_matrix(sendcounts_D[prank] / width, width, aIM_tD);
  alloc_matrix(sendcounts_T[prank] / height, height, aRE_tT);
  alloc_matrix(sendcounts_T[prank] / height, height, aIM_tT);

  // allocate temporary direct and transposed matrices B
  double **bRE_tD, **bIM_tD, **bRE_tT, **bIM_tT, **cRE_tT, **cIM_tT;
  alloc_matrix(sendcounts_D[prank] / width, width, bRE_tD);
  alloc_matrix(sendcounts_D[prank] / width, width, bIM_tD);
  alloc_matrix(sendcounts_T[prank] / height, height, bRE_tT);
  alloc_matrix(sendcounts_T[prank] / height, height, bIM_tT);
  alloc_matrix(sendcounts_T[prank] / height, height, cRE_tT);
  alloc_matrix(sendcounts_T[prank] / height, height, cIM_tT);

  // allocate temporary row/column
  double *arrRE = new double[max(height, width)]();
  double *arrIM = new double[max(height, width)]();

  // allocate twiddle
  double *wCOS = new double[max(height, width) / 2]();
  double *wSIN = new double[max(height, width) / 2]();

  // scatterv direct matrices A
  if (!prank)
    MPI_Scatterv(aRE[0], sendcounts_D, displs_D, MPI_DOUBLE, aRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, aRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(aIM[0], sendcounts_D, displs_D, MPI_DOUBLE, aIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, aIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // scatterv direct matrices A
  if (!prank)
    MPI_Scatterv(bRE[0], sendcounts_D, displs_D, MPI_DOUBLE, bRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, bRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(bIM[0], sendcounts_D, displs_D, MPI_DOUBLE, bIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, bIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // horizontal FFT A
  initTwiddle(wCOS, wSIN, width, prank, csize);
  for (int i = 0; i < sendcounts_D[prank] / width; ++i) {
    for (int j = 0; j < width; ++j) {
      // potentially make a bit_reversal array which could be scattered
      int reversed = bit_reversal(j, log2w);
      arrRE[j] = aRE_tD[i][reversed];
      arrIM[j] = aIM_tD[i][reversed];
    }

    fft(width, log2w, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < width; ++j) {
      aRE_tD[i][j] = arrRE[j];
      aIM_tD[i][j] = arrIM[j];
    }
  }

  // horizontal FFT B
  for (int i = 0; i < sendcounts_D[prank] / width; ++i) {
    for (int j = 0; j < width; ++j) {
      // potentially make a bit_reversal array which could be scattered
      int reversed = bit_reversal(j, log2w);
      arrRE[j] = bRE_tD[i][reversed];
      arrIM[j] = bIM_tD[i][reversed];
    }

    fft(width, log2w, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < width; ++j) {
      bRE_tD[i][j] = arrRE[j];
      bIM_tD[i][j] = arrIM[j];
    }
  }

  // gatherv direct matrices A
  if (!prank)
    MPI_Gatherv(aRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, aRE[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(aRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL, sendcounts_D,
                displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Gatherv(aIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, aIM[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(aIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL, sendcounts_D,
                displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // gatherv direct matrices B
  if (!prank)
    MPI_Gatherv(bRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, bRE[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(bRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL, sendcounts_D,
                displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Gatherv(bIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, bIM[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(bIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL, sendcounts_D,
                displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // transpose matrix A
  double **aRE_T = NULL, **aIM_T = NULL;
  if (!prank) {
    alloc_matrix(height, width, aRE_T);
    alloc_matrix(height, width, aIM_T);
    transpose(height, width, aRE, aRE_T);
    transpose(height, width, aIM, aIM_T);
  }

  // transpose matrix B
  double **bRE_T = NULL, **bIM_T = NULL;
  if (!prank) {
    alloc_matrix(height, width, bRE_T);
    alloc_matrix(height, width, bIM_T);
    transpose(height, width, bRE, bRE_T);
    transpose(height, width, bIM, bIM_T);
  }

  // scatterv transposed matrices A
  if (!prank)
    MPI_Scatterv(aRE_T[0], sendcounts_T, displs_T, MPI_DOUBLE, aRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, aRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(aIM_T[0], sendcounts_T, displs_T, MPI_DOUBLE, aIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, aIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // scatterv transposed matrices B
  if (!prank)
    MPI_Scatterv(bRE_T[0], sendcounts_T, displs_T, MPI_DOUBLE, bRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, bRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(bIM_T[0], sendcounts_T, displs_T, MPI_DOUBLE, bIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, bIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // vertical FFT A
  if (height != width) initTwiddle(wCOS, wSIN, height, prank, csize);
  for (int i = 0; i < sendcounts_T[prank] / height; ++i) {
    for (int j = 0; j < height; ++j) {
      int reversed = bit_reversal(j, log2h);
      arrRE[j] = aRE_tT[i][reversed];
      arrIM[j] = aIM_tT[i][reversed];
    }
    fft(height, log2h, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < height; ++j) {
      aRE_tT[i][j] = arrRE[j];
      aIM_tT[i][j] = arrIM[j];
    }
  }

  // vertical FFT B
  for (int i = 0; i < sendcounts_T[prank] / height; ++i) {
    for (int j = 0; j < height; ++j) {
      int reversed = bit_reversal(j, log2h);
      arrRE[j] = bRE_tT[i][reversed];
      arrIM[j] = bIM_tT[i][reversed];
    }
    fft(height, log2h, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < height; ++j) {
      bRE_tT[i][j] = arrRE[j];
      bIM_tT[i][j] = arrIM[j];
    }
  }

  // elementwise multiplication
  for (int i = 0; i < sendcounts_T[prank] / height; ++i) {
    for (int j = 0; j < height; ++j) {
      cRE_tT[i][j] = aRE_tT[i][j] * bRE_tT[i][j] - aIM_tT[i][j] * bIM_tT[i][j];
      cIM_tT[i][j] =
          -(aRE_tT[i][j] * bIM_tT[i][j] + aIM_tT[i][j] * bRE_tT[i][j]);
    }
  }

  // gatherv transposed matrices C_tT (we are not interested in A and B)
  if (!prank)
    MPI_Gatherv(cRE_tT[0], sendcounts_T[prank], MPI_DOUBLE, aRE_T[0],
                sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(cRE_tT[0], sendcounts_T[prank], MPI_DOUBLE, NULL, sendcounts_T,
                displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Gatherv(cIM_tT[0], sendcounts_T[prank], MPI_DOUBLE, aIM_T[0],
                sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(cIM_tT[0], sendcounts_T[prank], MPI_DOUBLE, NULL, sendcounts_T,
                displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // transpose matrices again A -> C (we could write C_tT -> C but that requires
  // more memory)
  if (!prank) {
    transpose(width, height, aRE_T, cRE);
    transpose(width, height, aIM_T, cIM);
  }

  // deallocate memory
  delete[] arrRE;
  delete[] arrIM;
  delete[] wCOS;
  delete[] wSIN;

  if (!prank) {
    dealloc_matrix(width, aRE_T);
    dealloc_matrix(width, aIM_T);

    dealloc_matrix(width, bRE_T);
    dealloc_matrix(width, bIM_T);
  }

  dealloc_matrix(sendcounts_D[prank] / width, aRE_tD);
  dealloc_matrix(sendcounts_D[prank] / width, aIM_tD);
  dealloc_matrix(sendcounts_T[prank] / height, aRE_tT);
  dealloc_matrix(sendcounts_T[prank] / height, aIM_tT);

  dealloc_matrix(sendcounts_D[prank] / width, bRE_tD);
  dealloc_matrix(sendcounts_D[prank] / width, bIM_tD);
  dealloc_matrix(sendcounts_T[prank] / height, bRE_tT);
  dealloc_matrix(sendcounts_T[prank] / height, bIM_tT);

  dealloc_matrix(sendcounts_T[prank] / height, cRE_tT);
  dealloc_matrix(sendcounts_T[prank] / height, cIM_tT);

  delete[] sendcounts_D;
  delete[] sendcounts_T;
  if (!prank) {
    delete[] displs_D;
    delete[] displs_T;
  }
}

void fft2_c(int const &height, int const &log2h, int const &width,
            int const &log2w, double *const *const matRE,
            double *const *const matIM, double *const *const c,
            int const &prank, int const &csize) {
  // allocate direct and transposed sendcounts and displs
  int *sendcounts_D = create_sendcounts(csize, height, width);
  int *sendcounts_T = create_sendcounts(csize, width, height);
  int *displs_D = NULL, *displs_T = NULL;
  if (!prank) {
    displs_D = create_displs(sendcounts_D, csize);
    displs_T = create_displs(sendcounts_T, csize);
  }

  // allocate temporary direct and transposed matrices
  double **matRE_tD, **matIM_tD, **matRE_tT, **matIM_tT;
  alloc_matrix(sendcounts_D[prank] / width, width, matRE_tD);
  alloc_matrix(sendcounts_D[prank] / width, width, matIM_tD);
  alloc_matrix(sendcounts_T[prank] / height, height, matRE_tT);
  alloc_matrix(sendcounts_T[prank] / height, height, matIM_tT);

  // allocate temporary row/column
  double *arrRE = new double[max(height, width)]();
  double *arrIM = new double[max(height, width)]();

  // allocate twiddle
  double *wCOS = new double[max(height, width) / 2]();
  double *wSIN = new double[max(height, width) / 2]();

  // scatterv direct matrices
  if (!prank)
    MPI_Scatterv(matRE[0], sendcounts_D, displs_D, MPI_DOUBLE, matRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, matRE_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(matIM[0], sendcounts_D, displs_D, MPI_DOUBLE, matIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_D, displs_D, MPI_DOUBLE, matIM_tD[0],
                 sendcounts_D[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // horizontal FFT
  initTwiddle(wCOS, wSIN, width, prank, csize);
  for (int i = 0; i < sendcounts_D[prank] / width; ++i) {
    for (int j = 0; j < width; ++j) {
      // potentially make a bit_reversal array which could be scattered
      int reversed = bit_reversal(j, log2w);
      arrRE[j] = matRE_tD[i][reversed];
      arrIM[j] = matIM_tD[i][reversed];
    }

    fft(width, log2w, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < width; ++j) {
      matRE_tD[i][j] = arrRE[j];
      matIM_tD[i][j] = arrIM[j];
    }
  }

  // gatherv direct matrices
  if (!prank)
    MPI_Gatherv(matRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, matRE[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(matRE_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL,
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Gatherv(matIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, matIM[0],
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(matIM_tD[0], sendcounts_D[prank], MPI_DOUBLE, NULL,
                sendcounts_D, displs_D, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // transpose matrix
  double **matRE_T = NULL, **matIM_T = NULL;
  if (!prank) {
    alloc_matrix(height, width, matRE_T);
    alloc_matrix(height, width, matIM_T);
    transpose(height, width, matRE, matRE_T);
    transpose(height, width, matIM, matIM_T);
  }

  // scatterv transposed matrices
  if (!prank)
    MPI_Scatterv(matRE_T[0], sendcounts_T, displs_T, MPI_DOUBLE, matRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, matRE_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (!prank)
    MPI_Scatterv(matIM_T[0], sendcounts_T, displs_T, MPI_DOUBLE, matIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Scatterv(NULL, sendcounts_T, displs_T, MPI_DOUBLE, matIM_tT[0],
                 sendcounts_T[prank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // vertical FFT
  int size2d = height * width;
  if (height != width) initTwiddle(wCOS, wSIN, height, prank, csize);
  for (int i = 0; i < sendcounts_T[prank] / height; ++i) {
    for (int j = 0; j < height; ++j) {
      int reversed = bit_reversal(j, log2h);
      arrRE[j] = matRE_tT[i][reversed];
      arrIM[j] = matIM_tT[i][reversed];
    }
    fft(height, log2h, arrRE, arrIM, wCOS, wSIN);

    // elementwise division
    for (int j = 0; j < height; ++j) arrRE[j] /= size2d;

    for (int j = 0; j < height; ++j) {
      matRE_tT[i][j] = arrRE[j];
      matIM_tT[i][j] = arrIM[j];
    }
  }

  // gatherv transposed matrices
  if (!prank)
    MPI_Gatherv(matRE_tT[0], sendcounts_T[prank], MPI_DOUBLE, matRE_T[0],
                sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  else
    MPI_Gatherv(matRE_tT[0], sendcounts_T[prank], MPI_DOUBLE, NULL,
                sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // if (!prank)
  //   MPI_Gatherv(matIM_tT[0], sendcounts_T[prank], MPI_DOUBLE, matIM_T[0],
  //               sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  // else
  //   MPI_Gatherv(matIM_tT[0], sendcounts_T[prank], MPI_DOUBLE, NULL,
  //               sendcounts_T, displs_T, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // transpose matrices again
  if (!prank) {
    transpose(width, height, matRE_T, c);
    // transpose(width, height, matIM_T, matIM);
  }

  // deallocate memory
  delete[] arrRE;
  delete[] arrIM;
  delete[] wCOS;
  delete[] wSIN;

  if (!prank) {
    dealloc_matrix(width, matRE_T);
    dealloc_matrix(width, matIM_T);
  }
  dealloc_matrix(sendcounts_D[prank] / width, matRE_tD);
  dealloc_matrix(sendcounts_D[prank] / width, matIM_tD);
  dealloc_matrix(sendcounts_T[prank] / height, matRE_tT);
  dealloc_matrix(sendcounts_T[prank] / height, matIM_tT);
  delete[] sendcounts_D;
  delete[] sendcounts_T;
  if (!prank) {
    delete[] displs_D;
    delete[] displs_T;
  }
}

void debug(int const &ch, int const &cw, double const *const *const aRE,
           double const *const *const aIM, double const *const *const bRE,
           double const *const *const bIM, int const &prank) {
  if (!prank) {
    cout << "aRE" << endl;
    print_matrix(ch, cw, aRE);
    cout << "aIM" << endl;
    print_matrix(ch, cw, aIM);
    cout << "bRE" << endl;
    print_matrix(ch, cw, bRE);
    cout << "bIM" << endl;
    print_matrix(ch, cw, bIM);
  }
}

void conv2(int const &ah, int const &aw, int const &bh, int const &bw, int &ch,
           int &cw, double const *const *const a, double const *const *const b,
           double **&c, int const &prank, int const &csize) {
  // allocate memory
  double **aRE, **aIM, **bRE, **bIM, **cRE, **cIM;

  // initialize height & width
  int log2w = int(ceil(log2(aw + bw - 1)));
  int log2h = int(ceil(log2(ah + bh - 1)));
  cw = (1 << log2w);
  ch = (1 << log2h);

  if (!prank) {
    alloc_matrix(ch, cw, aIM);
    alloc_matrix(ch, cw, aRE);
    alloc_matrix(ch, cw, bRE);
    alloc_matrix(ch, cw, bIM);
    alloc_matrix(ch, cw, cRE);
    alloc_matrix(ch, cw, cIM);
    alloc_matrix(ch, cw, c);
    copy_matrix(ah, aw, a, aRE);
    copy_matrix(bh, bw, b, bRE);
  }

  // direct fft2 && elementwise multiplication
  fft2_a_and_b(ch, log2h, cw, log2w, aRE, aIM, bRE, bIM, cRE, cIM, prank,
               csize);

  // inverse fft2 && elementwise division
  fft2_c(ch, log2h, cw, log2w, cRE, cIM, c, prank, csize);

  // deallocate memory
  if (!prank) {
    dealloc_matrix(ch, aRE);
    dealloc_matrix(ch, aIM);
    dealloc_matrix(ch, bRE);
    dealloc_matrix(ch, bIM);
    dealloc_matrix(ch, cRE);
    dealloc_matrix(ch, cIM);
  }

  ch = ah + bh - 1;
  cw = aw + bw - 1;
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

  // cout << "Elapsed time distributed (" << csize
  //      << " processes) = " << long(mind * 1000000000) << " [ns] " << endl;
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

  conv2(ah, aw, bh, bw, ch, cw, matrix, kernel, conv, prank, csize);

  // end time measurement
  double end_time = MPI_Wtime();
  double total_time = end_time - start_time;
  double mind;
  MPI_Reduce(&total_time, &mind, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (!prank) {
    char const matrix_file_name[] = "outputs/distributed_matrix.txt";
    char const time_file_name[] = "outputs/distributed_time.txt";
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
