#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

#include "matrix_operators.hpp"

using namespace std;

void conv2(int const &ah, int const &aw, int const &bh, int const &bw, int &ch,
           int &cw, double const *const *const a, double const *const *const b,
           double **&c) {
  // initialize height & width
  cw = aw + bw - 1;
  ch = ah + bh - 1;

  // allocate memory
  alloc_matrix(ch, cw, c);

  // convolution
  for (int ci = 0; ci < ch; ++ci)
    for (int cj = 0; cj < cw; ++cj) {
      // clang-format off
      int aw_start = std::max(0, cj - bw + 1),
          // bw_start = std::max(0, cj - aw + 1),
          aw_end = std::min(aw, cj + 1),
          bw_end = std::min(bw, cj + 1),
          ah_start = std::max(0, ci - bh + 1),
          // bh_start = std::max(0, ci - ah + 1),
          ah_end = std::min(ah, ci + 1),
          bh_end = std::min(bh, ci + 1),
          wlen = aw_end - aw_start, // same as bw_end - bw_start
          hlen = ah_end - ah_start; // same as bh_end - bh_start
      // clang-format on
      for (int i = 0; i < hlen; ++i)
        for (int j = 0; j < wlen; ++j)
          c[ci][cj] +=
            a[ah_start + i][aw_start + j] * b[bh_end - i - 1][bw_end - j - 1];
    }
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
  print_matrix(ch, cw, conv);
  // fprint_matrix(matrix_file_name, ch, cw, conv);

  // ofstream f(time_file_name, std::ios_base::app);
  // f << std::chrono::duration_cast<std::chrono::nanoseconds>(end -
  // begin).count()
  //   << std::endl;
  // f.close();
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

  char const matrix_file_name[] = "outputs/regular_serial_matrix.txt";
  char const time_file_name[] = "outputs/regular_serial_time.txt";
  write(matrix_file_name, time_file_name, begin, end, ch, cw, conv);

  dealloc_matrix(ah, matrix);
  dealloc_matrix(bh, kernel);
  dealloc_matrix(ch, conv);

  return 0;
}
