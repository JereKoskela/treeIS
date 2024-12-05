#ifndef MATRIX
#define MATRIX

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

class Matrix {

public:
  Matrix(const int nr, const int nc, const double entry)
      : nrow(nr), ncol(nc), entries(nr * nc, entry) {}

  Matrix(const int n) : nrow(n), ncol(n), entries(n * n, 0.0) {}

  void copy(const Matrix &other) {
    nrow = other.nrow;
    ncol = other.ncol;
    entries = other.entries;
  }

  double getEntry(const int row, const int col) const {
    return entries[row * ncol + col];
  }

  void setEntry(const int row, const int col, const double entry) {
    entries[row * ncol + col] = entry;
    return;
  }

  std::vector<double> getRow(const int row) const {
    std::vector<double> ret(ncol);
    for (int col = 0; col < ncol; col++) {
      ret[col] = entries[row * ncol + col];
    }
    return ret;
  }

  int nrow;
  int ncol;
  std::vector<double> entries;
};

#endif
