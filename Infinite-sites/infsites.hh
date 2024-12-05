#ifndef InfSites
#define InfSites

#include "matrix.hh"
#include <cassert>
#include <fstream>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_spmatrix.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct infSites {

  infSites(std::string filename, const double th)
      : theta(th), total_prob(0), sample_size(0), row_counts(), col_counts(),
        col_sums(), num_sites(), row_prop(), data(1) {
    std::ifstream file;
    file.open(filename);
    std::string line, token;
    gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, 1);
    int l = 0;
    while (getline(file, line)) {
      std::vector<int> tmp_vec;
      std::stringstream iss;
      iss << line;
      while (getline(iss, token, ' ')) {
        tmp_vec.push_back(atoi(token.c_str()));
      }
      row_counts.push_back(tmp_vec.back());
      row_prop.push_back(0);
      sample_size += row_counts.back();
      if (l == 0) {
        col_sums.resize(tmp_vec.size() - 1, 0);
        col_counts.resize(tmp_vec.size() - 1, 0);
        num_sites.resize(tmp_vec.size() - 1, 1);
      }
      for (unsigned int i = 0; i < tmp_vec.size() - 1; i++) {
        if (tmp_vec[i] == 1) {
          gsl_spmatrix_set(tmp, l, i, 1);
          col_sums[i] += 1;
          col_counts[i] += row_counts.back();
        } else {
          gsl_spmatrix_set(tmp, l, i, 0);
        }
      }
      l++;
    }
    data.resize(int(tmp->size1), int(tmp->size2));
    for (int i = 0; i < int(tmp->size1); i++) {
      for (int j = 0; j < int(tmp->size2); j++) {
        data.setEntry(i, j, gsl_spmatrix_get(tmp, i, j));
      }
    }
    gsl_spmatrix_free(tmp);
    int match, i, r;
    for (int j = data.ncol; j > -1; j--) {
      match = 0;
      r = 0;
      while (match == 0 && r < j) {
        i = 0;
        match = 1;
        while (match == 1 && i < data.nrow) {
          if (data.getEntry(i, j) != data.getEntry(i, r)) {
            match = 0;
          } else {
            i++;
          }
        }
        if (match == 0) {
          r++;
        }
      }
      if (match == 1) {
        num_sites[r]++;
        data.removeColumn(j);
        col_sums.erase(col_sums.begin() + j);
        col_counts.erase(col_counts.begin() + j);
        num_sites.erase(num_sites.begin() + j);
      }
    }
  }

  infSites(const int n, const double m)
      : theta(m), total_prob(0), sample_size(n), row_counts(), col_counts(),
        col_sums(), num_sites(), row_prop(), data(1) {}

  infSites(const infSites &other)
      : theta(other.theta), total_prob(other.total_prob),
        sample_size(other.sample_size), row_counts(other.row_counts),
        col_counts(other.col_counts), col_sums(other.col_sums),
        num_sites(other.num_sites), row_prop(other.row_prop), data(other.data) {
  }

  void assign(const infSites &other) {
    theta = other.theta;
    total_prob = other.total_prob;
    sample_size = other.sample_size;
    row_counts = other.row_counts;
    col_counts = other.col_counts;
    col_sums = other.col_sums;
    num_sites = other.num_sites;
    row_prop = other.row_prop;
    data.copy(other.data);
  }

  int find_singleton(const int i) {
    for (int j = 0; j < data.ncol; j++) {
      if (data.getEntry(i, j) == 1 && col_sums[j] == 1) {
        return 1;
      }
    }
    return 0;
  }

  int singleton_count(const int i) {
    int ret = 0;
    for (int j = 0; j < data.ncol; j++) {
      if (data.getEntry(i, j) == 1 && col_sums[j] == 1 && row_counts[i] == 1) {
        ret++;
      }
    }
    return ret;
  }

  int sample_singleton(const int i, gsl_rng *gen) {
    int ret = -1;
    std::vector<double> w(data.ncol, 0);
    double w_sum = 0;
    std::vector<int> singletons;
    for (int j = 0; j < data.ncol; j++) {
      if (data.getEntry(i, j) == 1 && col_sums[j] == 1) {
        w[j] = num_sites[j];
        w_sum += w[j];
      }
    }
    if (w_sum > 0) {
      double coin = gsl_rng_uniform(gen);
      int ind = 0;
      double ub = w[0] / w_sum;
      while (ub < coin) {
        ind++;
        ub += w[ind] / w_sum;
      }
      ret = ind;
    }
    return ret;
  }

  double proposal(const Matrix &prop, const int i, const int j) {
    if (row_counts[i] == sample_size) {
      return 1;
    }
    double ret = 0;
    if (row_counts[i] > 1 || find_singleton(i) == 1) {
      ret = 1;
      if (sample_size > 2 || row_counts[i] > 1) {
        ret = prop.getEntry(sample_size - 2, col_counts[j] - 1);
      }
      if (data.getEntry(i, j) == 1) {
        ret *= num_sites[j] * row_counts[i] / double(col_counts[j]);
      } else {
        ret = num_sites[j] * (1 - ret) * row_counts[i] /
              double(sample_size - col_counts[j]);
      }
    }
    return ret;
  }

  int sample_row(gsl_rng *gen) {
    int ret = 0;
    double coin = gsl_rng_uniform(gen);
    double ub = row_prop[0] / total_prob;
    while (coin > ub) {
      ret++;
      ub += row_prop[ret] / total_prob;
    }
    return ret;
  }

  void recompute_proposal(const Matrix &prop) {
    if (data.nrow == 1) {
      row_prop[0] = 1;
      total_prob = 1;
    } else {
      total_prob = 0;
      for (int i = 0; i < data.nrow; i++) {
        row_prop[i] = 0;
        if (row_counts[i] > 1 || find_singleton(i) == 1) {
          for (int j = 0; j < data.ncol; j++) {
            row_prop[i] += proposal(prop, i, j);
          }
          total_prob += row_prop[i];
        }
      }
    }
    return;
  }

  int match(const int i, const int j) const {
    double diff = 0;
    int r = 0;
    int c;
    while (r < data.nrow) {
      if (r != i) {
        diff = 0;
        c = 0;
        while (diff == 0 && c < data.ncol) {
          if (c != j) {
            diff += fabs(data.getEntry(i, c) - data.getEntry(r, c));
          } else {
            diff += data.getEntry(r, c);
          }
          c++;
        }
        if (diff == 0) {
          return r;
        }
      }
      r++;
    }
    return r;
  }

  int match(const int i) const {
    double diff = 0;
    int r = 0;
    int c;
    while (r < data.nrow) {
      if (r != i) {
        diff = 0;
        c = 0;
        while (diff == 0 && c < data.ncol) {
          diff += fabs(data.getEntry(i, c) - data.getEntry(r, c));
          c++;
        }
        if (diff == 0) {
          return r;
        }
      }
      r++;
    }
    return r;
  }

  int remove_mutation(const int i, const int j) {
    int ret = i;
    if (num_sites[j] > 1) {
      num_sites[j]--;
    } else {
      data.setEntry(i, j, 0);
      data.removeColumn(j);
      col_sums.erase(col_sums.begin() + j);
      col_counts.erase(col_counts.begin() + j);
      num_sites.erase(num_sites.begin() + j);
      int r = match(i);
      if (r < data.nrow) {
        row_counts[r] += row_counts[i];
        for (int c = 0; c < data.ncol; c++) {
          col_sums[c] -= data.getEntry(i, c);
        }
        row_counts.erase(row_counts.begin() + i);
        data.removeRow(i);
        if (i < r) {
          r--;
        }
        ret = r;
      }
    }
    return ret;
  }

  int remove_mutation(const int i, const int j, const Matrix &prop) {
    int ret = i;
    for (int r = 0; r < data.nrow; r++) {
      total_prob -= proposal(prop, r, j);
      row_prop[r] -= proposal(prop, r, j);
    }
    if (num_sites[j] > 1) {
      num_sites[j]--;
      for (int r = 0; r < data.nrow; r++) {
        row_prop[r] += proposal(prop, r, j);
        total_prob += proposal(prop, r, j);
      }
    } else {
      data.setEntry(i, j, 0);
      data.removeColumn(j);
      col_sums.erase(col_sums.begin() + j);
      col_counts.erase(col_counts.begin() + j);
      num_sites.erase(num_sites.begin() + j);
      int r = match(i);
      if (r < data.nrow) {
        total_prob -= row_prop[r] + row_prop[i];
        row_counts[r] += row_counts[i];
        for (int c = 0; c < data.ncol; c++) {
          col_sums[c] -= data.getEntry(i, c);
        }
        row_counts.erase(row_counts.begin() + i);
        row_prop.erase(row_prop.begin() + i);
        data.removeRow(i);
        if (i < r) {
          r--;
        }
        if (data.ncol == 0) {
          row_prop[r] = 1;
        } else {
          row_prop[r] = 0;
          for (int c = 0; c < data.ncol; c++) {
            row_prop[r] += proposal(prop, r, c);
          }
        }
        total_prob += row_prop[r];
        ret = r;
      } else {
        total_prob -= row_prop[i];
        row_prop[i] = 0;
      }
    }
    return ret;
  }

  void generate_data(gsl_rng *gen) {
    int n = 2;
    gsl_spmatrix *tmp = gsl_spmatrix_alloc(1, 1);
    row_counts.push_back(2);
    int child = 0;
    int ub = row_counts[0];
    double coin = 0;
    int next_segregating_site = 0;
    while (n < sample_size + 1) {
      child = 0;
      ub = row_counts[0];
      coin = gsl_rng_uniform(gen);
      while (coin > double(ub) / n) {
        child++;
        ub += row_counts[child];
      }
      if (gsl_rng_uniform(gen) < theta / (n - 1 + theta)) {
        if (row_counts[child] == 1) {
          gsl_spmatrix_set(tmp, child, next_segregating_site, 1);
        } else {
          row_counts[child]--;
          row_counts.push_back(1);
          int new_child = row_counts.size() - 1;
          for (int i = 0; i < next_segregating_site; i++) {
            gsl_spmatrix_set(tmp, new_child, i,
                             gsl_spmatrix_get(tmp, child, i));
          }
          gsl_spmatrix_set(tmp, new_child, next_segregating_site, 1);
        }
        col_sums.push_back(1);
        next_segregating_site++;
      } else {
        if (n == sample_size) {
          break;
        }
        for (int i = 0; i < next_segregating_site; i++) {
          col_sums[i] += gsl_spmatrix_get(tmp, child, i);
        }
        row_counts[child]++;
        n++;
      }
    }
    data.resize(int(tmp->size1), int(tmp->size2));
    for (int i = 0; i < int(tmp->size1); i++) {
      for (int j = 0; j < int(tmp->size2); j++) {
        data.setEntry(i, j, gsl_spmatrix_get(tmp, i, j));
      }
    }
    gsl_spmatrix_free(tmp);
    return;
  }

  int sample_parent(gsl_rng *gen) {
    int ret = 0;
    double coin = gsl_rng_uniform(gen);
    double ub = double(row_counts[0]) / sample_size;
    while (coin > ub) {
      ret++;
      ub += double(row_counts[ret]) / sample_size;
    }
    return ret;
  }

  void kill(const int row, const int n) {
    row_counts[row] -= n;
    sample_size -= n;
    if (row_counts[row] == 0) {
      for (int i = 0; i < data.ncol; i++) {
        col_sums[i] -= data.getEntry(row, i);
      }
      data.removeRow(row);
      row_counts.erase(row_counts.begin() + row);
      for (int i = data.ncol - 1; i > -1; i--) {
        if (col_sums[i] == 0 || col_sums[i] == data.nrow) {
          data.removeColumn(i);
          col_sums.erase(col_sums.begin() + i);
        }
      }
    }
    return;
  }

  void print() const {
    for (int i = 0; i < data.nrow; i++) {
      for (int j = 0; j < data.ncol; j++) {
        std::cout << data.getEntry(i, j) << " ";
      }
      std::cout << row_counts[i] << std::endl;
    }
    for (int i = 0; i < data.ncol; i++) {
      std::cout << col_sums[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < data.ncol; i++) {
      std::cout << col_counts[i] << " ";
    }
    std::cout << std::endl;
    return;
  }

  void print(char *fileName) const {
    std::ofstream outputFile;
    outputFile.open(fileName);
    for (int i = 0; i < data.nrow; i++) {
      for (int j = 0; j < data.ncol; j++) {
        outputFile << data.getEntry(i, j) << " ";
      }
      outputFile << row_counts[i] << std::endl;
    }
    outputFile.close();
    return;
  }

  double theta, total_prob;
  int sample_size;
  std::vector<int> row_counts, col_counts, col_sums, num_sites;
  std::vector<double> row_prop;
  Matrix data;
};

#endif
