#pragma once
#include <cmath>
#include <cstdio>
// #include "colors.hpp"
#define PI 3.141592653589793238462643383279502884197

// https://stackoverflow.com/questions/7546620/operator-new-initializes-memory-to-zero
// https://stackoverflow.com/questions/1403150/how-do-you-dynamically-allocate-a-matrix

// matrix is one memory block
void inline alloc_matrix(int const height, int const width, double **&matrix) {
  matrix = new double *[height];
  if (height) {
    matrix[0] = new double[height * width]();
    for (int i = 1; i < height; ++i) matrix[i] = matrix[0] + i * width;
  }
}
void inline dealloc_matrix(int const height, double **&matrix) {
  if (height) delete[] matrix[0];
  delete[] matrix;
}

// row is one memory block
void inline alloc_matrix_separated(int const height, int const width,
                                   double **&matrix) {
  matrix = new double *[height];
  for (int i = 0; i < height; i++) matrix[i] = new double[width]();
}
void inline dealloc_matrix_separated(int const height, double **const matrix) {
  for (int i = 0; i < height; i++) delete[] matrix[i];
  delete[] matrix;
}

void inline transpose(int const &height, int const &width,
                      double const *const *const matrix,
                      double *const *const &transposed) {
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j) transposed[j][i] = matrix[i][j];
}

// https://isocpp.org/wiki/faq/const-correctness#constptrptr-conversion
void inline copy_matrix(int const &height, int const &width,
                        double const *const *const src,
                        double *const *const dest) {
  for (int i = 0; i < height; i++)
    for (int j = 0; j < width; j++) dest[i][j] = src[i][j];
}

void inline print_matrix(int const height, int const width,
                         double const *const *const matrix) {
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) printf("%.2f ", (matrix)[i][j]);
    printf("\n");
  }
}
void inline fprint_matrix(char const *const file_name, int const height,
                          int const width, double const *const *const matrix) {
  FILE *f = fopen(file_name, "w");
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) fprintf(f, "%.2f ", (matrix)[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}
void inline fprint_matrix_append(char const *const file_name, int const height,
                                 int const width,
                                 double const *const *const matrix) {
  FILE *f = fopen(file_name, "a");
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) fprintf(f, "%.2f ", (matrix)[i][j]);
    fprintf(f, "\n");
  }
  fclose(f);
}

void inline empty_file(char const *const file_name) {
  fclose(fopen(file_name, "w"));
}

void inline print_array(int const &length, double const *const array) {
  // printf(RESET);
  for (int i = 0; i < length; ++i) printf("%.2f ", array[i]);
  printf("\n");
}

void inline print_array(int const &length, int const *const array) {
  // printf(RESET);
  for (int i = 0; i < length; ++i) printf("%d ", array[i]);
  printf("\n");
}

// void inline print_array(int const length, double const *const array, char
// const *const color) {
//   printf(color);
//   for (int i = 0; i < length; ++i) printf("%.2f ", array[i]);
//   printf("\n");
//   printf(RESET);
// }

// void inline print_array(int const length, int const *const array, char const
// *const color) {
//   printf(color);
//   for (int i = 0; i < length; ++i) printf("%d ", array[i]);
//   printf("\n");
//   printf(RESET);
// }
