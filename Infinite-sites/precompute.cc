#include "matrix.hh"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_sf.h>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cout << "Call " << argv[0] << " <sample_size> <theta>" << std::endl;
    return 1;
  }
  int sample_size = atoi(argv[1]);
  double theta = atof(argv[2]);

  Matrix prop(sample_size - 1, sample_size - 1, 0.0);
  double num;
  double den;
  int n, d;
  for (int i = 0; i <= sample_size - 2; i++) {
    n = i + 2;
    d = 1;
    num = 1 / (n - 1 + theta);
    den = 0;
    for (int k = 2; k <= n; k++) {
      den += (k - 1) / ((k - 1 + theta) * (n - 1));
    }
    prop.setEntry(i, 0, num / den);
    for (int j = 1; j <= n - 2; j++) {
      num = 0;
      den = 0;
      d = j + 1;
      for (int k = 2; k <= n - d + 1; k++) {
        num += (d - 1) *
               exp(gsl_sf_lnchoose(n - d - 1, k - 2) -
                   gsl_sf_lnchoose(n - 1, k - 1)) /
               ((n - k) * (k - 1 + theta));
        den += exp(gsl_sf_lnchoose(n - d - 1, k - 2) -
                   gsl_sf_lnchoose(n - 1, k - 1)) /
               (k - 1 + theta);
      }
      prop.setEntry(i, j, num / den);
    }
  }
  prop.print();

  return 0;
}
