#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "matrix_operators.hpp"

using namespace std;

void initTwiddle(double **const wCOS, double **const wSIN, int const size) {
  for (int i = 0; i < size / 2; ++i) {
    (*wCOS)[i] = cos(-2.0 * PI * i / size);
    (*wSIN)[i] = sin(-2.0 * PI * i / size);
  }
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

void fft2(int const &height, int const &log2h, int const &width,
          int const &log2w, double *const *const matRE,
          double *const *const matIM) {
  // allocate temorary row/column
  double *arrRE = new double[max(height, width)]();
  double *arrIM = new double[max(height, width)]();

  // allocate twiddle
  double *wCOS = new double[max(height, width) / 2]();
  double *wSIN = new double[max(height, width) / 2]();

  // horizontal FFT
  initTwiddle(&wCOS, &wSIN, width);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int reversed = bit_reversal(j, log2w);
      arrRE[j] = matRE[i][reversed];
      arrIM[j] = matIM[i][reversed];
    }
    fft(width, log2w, arrRE, arrIM, wCOS, wSIN);

    for (int j = 0; j < width; ++j) {
      matRE[i][j] = arrRE[j];
      matIM[i][j] = arrIM[j];
    }
  }

  // vertical FFT
  if (height != width) initTwiddle(&wCOS, &wSIN, height);
  for (int j = 0; j < width; ++j) {
    for (int i = 0; i < height; ++i) {
      int reversed = bit_reversal(i, log2h);
      arrRE[i] = matRE[reversed][j], arrIM[i] = matIM[reversed][j];
    }
    fft(height, log2h, arrRE, arrIM, wCOS, wSIN);

    for (int i = 0; i < height; ++i) {
      matRE[i][j] = arrRE[i];
      matIM[i][j] = arrIM[i];
    }
  }

  // deallocate memory
  delete[] arrRE;
  delete[] arrIM;
  delete[] wCOS;
  delete[] wSIN;
}

void debug(int const &ch, int const &cw, double const *const *const aRE,
           double const *const *const aIM, double const *const *const bRE,
           double const *const *const bIM) {
  cout << "aRE" << endl;
  print_matrix(ch, cw, aRE);
  cout << "aIM" << endl;
  print_matrix(ch, cw, aIM);
  cout << "bRE" << endl;
  print_matrix(ch, cw, bRE);
  cout << "bIM" << endl;
  print_matrix(ch, cw, bIM);
}

void conv2(int const &ah, int const &aw, int const &bh, int const &bw, int &ch,
           int &cw, double const *const *const a, double const *const *const b,
           double **&c) {
  // initialize height & width
  int log2w = int(ceil(log2(aw + bw - 1)));
  int log2h = int(ceil(log2(ah + bh - 1)));
  cw = (1 << log2w);
  ch = (1 << log2h);

  // allocate memory
  double **aRE, **aIM, **bRE, **bIM, **cRE, **cIM;
  alloc_matrix(ch, cw, aRE);
  alloc_matrix(ch, cw, bRE);
  alloc_matrix(ch, cw, aIM);
  alloc_matrix(ch, cw, bIM);
  alloc_matrix(ch, cw, cRE);
  alloc_matrix(ch, cw, cIM);
  alloc_matrix(ch, cw, c);

  copy_matrix(ah, aw, a, aRE);
  copy_matrix(bh, bw, b, bRE);

  // direct fft2
  fft2(ch, log2h, cw, log2w, aRE, aIM);
  fft2(ch, log2h, cw, log2w, bRE, bIM);

  // convolution -> elementwise multiplication
  for (int i = 0; i < ch; ++i)
    for (int j = 0; j < cw; ++j) {
      cRE[i][j] = aRE[i][j] * bRE[i][j] - aIM[i][j] * bIM[i][j];
      cIM[i][j] = -(aRE[i][j] * bIM[i][j] + aIM[i][j] * bRE[i][j]);
    }

  // inverse fft2
  fft2(ch, log2h, cw, log2w, cRE, cIM);
  for (int i = 0; i < ch; ++i)
    for (int j = 0; j < cw; ++j) c[i][j] = cRE[i][j] / (cw * ch);

  // deallocate memory
  dealloc_matrix(ch, aRE);
  dealloc_matrix(ch, aIM);
  dealloc_matrix(ch, bRE);
  dealloc_matrix(ch, bIM);
  dealloc_matrix(ch, cRE);
  dealloc_matrix(ch, cIM);

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
           std::chrono::steady_clock::time_point const begin,
           std::chrono::steady_clock::time_point const end, int const &ch,
           int const &cw, double const *const *const conv) {
  // cout << "ch = " << ch << " cw = " << cw << endl;
  // print_matrix(ch, cw, conv);
  fprint_matrix(matrix_file_name, ch, cw, conv);

  // std::cout << "Elapsed time serial = " <<
  // std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
  // << "[ns]" << std::endl;
  ofstream f(time_file_name, std::ios_base::app);
  f << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()
    << std::endl;
  f.close();
}

int main() {
  double **matrix, **kernel, **conv;
  int ah, aw, bh, bw, ch, cw;
  const char file_name[] = "input/input.txt";
  read(file_name, ah, aw, bh, bw, matrix, kernel);

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  conv2(ah, aw, bh, bw, ch, cw, matrix, kernel, conv);
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  char const matrix_file_name[] = "outputs/serial_matrix.txt";
  char const time_file_name[] = "outputs/serial_time.txt";
  write(matrix_file_name, time_file_name, begin, end, ch, cw, conv);

  dealloc_matrix(ah, matrix);
  dealloc_matrix(bh, kernel);
  dealloc_matrix(ch, conv);

  return 0;
}
