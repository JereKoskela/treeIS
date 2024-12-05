#ifndef MATRIX
#define MATRIX

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Matrix {

  Matrix(std::string filename) : nrow(0), ncol(), entries() {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    while (getline(file, line)) {
      nrow++;
      std::stringstream iss;
      iss << line;
      while (getline(iss, token, ' ')) {
        entries.push_back(atof(token.c_str()));
      }
      if (nrow == 1) {
        ncol = int(entries.size());
      }
    }
  }

  Matrix(const int nr, const int nc, const double entry)
      : nrow(nr), ncol(nc), entries(nr * nc, entry) {}

  Matrix(const int n) : nrow(n), ncol(n), entries(n * n, 0.0) {}

  void copy(const Matrix &other) {
    nrow = other.nrow;
    ncol = other.ncol;
    entries = other.entries;
    return;
  }

  void resize(const int nr, const int nc) {
    nrow = nr;
    ncol = nc;
    entries.resize(nr * nc);
    return;
  }

  double getEntry(const int row, const int col) const {
    return entries[row * ncol + col];
  }

  void setEntry(const int row, const int col, const double entry) {
    entries[row * ncol + col] = entry;
    return;
  }

  void removeRow(const int row) {
    entries.erase(entries.begin() + row * ncol,
                  entries.begin() + (row + 1) * ncol);
    nrow--;
    return;
  }

  void removeColumn(const int col) {
    for (int i = nrow - 1; i > -1; i--) {
      entries.erase(entries.begin() + i * ncol + col);
    }
    ncol--;
    return;
  }

  void print() const {
    for (int i = 0; i < nrow; i++) {
      for (int j = 0; j < ncol; j++) {
        std::cout << entries[i * ncol + j] << " ";
      }
      std::cout << std::endl;
    }
    return;
  }

  int nrow, ncol;
  std::vector<double> entries;
};

#endif
