#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>

#include "matrix_operators.hpp"

using namespace std;

void read(char const *const file_name, int &ah, int &aw, int &bh, int &bw,
          char *&kernel_name) {
  ifstream f(file_name);
  f >> ah >> aw >> bh >> bw >> kernel_name;
  f.close();
}

void create_matrices(int const &ah, int const &aw, int const &bh, int const &bw,
                     char *kernel_name, double **&a, double **&b) {
  alloc_matrix(ah, aw, a);
  alloc_matrix(bh, bw, b);
  for (int i = 0; i < ah; ++i)
    for (int j = 0; j < aw; ++j) a[i][j] = double(rand() % 255);

  if (!strcmp(kernel_name, "sobelx")) {
    b[0][0] = 1;
    b[0][2] = -1;
    b[1][0] = 2;
    b[1][2] = -2;
    b[2][0] = 1;
    b[2][2] = -1;
  } else if (!strcmp(kernel_name, "sobely")) {
    b[0][0] = 1;
    b[0][1] = 2;
    b[0][2] = 1;
    b[2][0] = -1;
    b[2][1] = -2;
    b[2][2] = -1;
  } else
    for (int i = 0; i < bh; ++i)
      for (int j = 0; j < bw; ++j) b[i][j] = double(rand() % 255);
}

void write(char const *const file_name, int const &ah, int const &aw,
           int const &bh, int const &bw, double const *const *const &a,
           double const *const *const b) {
  ofstream f(file_name);
  f << ah << " " << aw << " " << bh << " " << bw << endl;
  for (int i = 0; i < ah; ++i) {
    for (int j = 0; j < aw; ++j) f << int(a[i][j]) << " ";
    f << endl;
  }

  for (int i = 0; i < bh; ++i) {
    for (int j = 0; j < bw; ++j) f << int(b[i][j]) << " ";
    f << endl;
  }
  f.close();
}

int main() {
  int ah, aw, bh, bw;
  double **a, **b;
  srand(time(NULL));
  char *kernel_name = new char[10];

  read("mat_size.txt", ah, aw, bh, bw, kernel_name);
  create_matrices(ah, aw, bh, bw, kernel_name, a, b);
  write("input/input.txt", ah, aw, bh, bw, a, b);

  delete [] kernel_name;
  dealloc_matrix(ah, a);
  dealloc_matrix(bh, b);
  return 0;
}
