#pragma once
#include <mpi.h>

int *create_sendcounts(int const &csize, int const &height, int const &width) {
  int *sendcounts = new int[csize]();
  for (int i = 0; i < csize; ++i) sendcounts[i] = width * (height / csize);
  for (int i = 0; i < height % csize; ++i) sendcounts[i] += width;
  return sendcounts;
}

int *create_recvcounts(int const &csize, int const &dim) {
  int *recvcounts = new int[dim]();
  for (int i = 0; i < csize; ++i) recvcounts[i] = dim / csize;
  for (int i = 0; i < dim % csize; ++i) recvcounts[i] += 1;
  return recvcounts;
}

int *create_displs(int const *sendcounts, int const &csize) {
  int *displs = new int[csize]();
  for (int i = 1; i < csize; ++i) displs[i] = displs[i - 1] + sendcounts[i - 1];
  return displs;
}

int *create_elementwise_sendcounts(int const &csize, int const &length) {
  int *sendcounts = new int[csize]();
  for (int i = 0; i < csize; ++i) sendcounts[i] = length / csize;
  for (int i = 0; i < length % csize; ++i) ++sendcounts[i];
  return sendcounts;
}

int *create_elementwise_recvcounts(int const &csize, int const &length) {
  int *recvcounts = new int[csize]();
  for (int i = 0; i < csize; ++i) recvcounts[i] = length / csize;
  for (int i = 0; i < length % csize; ++i) ++recvcounts[i];
  return recvcounts;
}
