#ifndef APPROX_CSD
#define APPROX_CSD

#include "finite_sites.hh"
#include "fs_sample.hh"
#include "matrix.hh"
#include <gsl/gsl_sf_gamma.h>
#include <iostream>

class CSD {
public:
  CSD(const double theta_, const int sample_size, const FiniteSites &mutation)
      : theta(theta_), quad_points(4), row_length(sample_size - 1),
        M(4 * (sample_size - 1), 2) {
    quad_points[0] = 0.322548;
    quad_points[1] = 1.74576;
    quad_points[2] = 4.53662;
    quad_points[3] = 9.39507;

    double rate, temp;
    for (unsigned int point = 0; point < quad_points.size(); point++) {
      for (int n = 0; n < sample_size - 1; n++) {
        rate = (double(n) + 1.0) / 2.0;
        for (int a = 0; a < 2; a++) {
          for (int b = 0; b < 2; b++) {
            temp = 0.0;
            for (int m = 0; m < 20; m++) {
              if (((m % 2 == 0) && (a == b)) || ((m % 2 != 0) && (a != b))) {
                temp +=
                    exp(double(m) * (log(theta / double(2 * mutation.site_no)) +
                                     log(quad_points[point]) - log(rate)) -
                        gsl_sf_lnfact(m) -
                        exp(log(theta / double(2 * mutation.site_no)) +
                            log(quad_points[point]) - log(rate)));
              }
            }
            M[point * row_length + n].setEntry(a, b, temp);
          }
        }
      }
    }
  }

  double pi_hat(const unsigned long long int type, const FSSample &sample,
                const FiniteSites &mutation) const {
    double ret = 0.0;
    const std::vector<int> type_vec = mutation.int_to_vec(type);
    std::vector<int> temp_type(mutation.site_no);
    for (int type_index = 0; type_index < sample.sample_types; type_index++) {
      temp_type = mutation.int_to_vec(sample.types[type_index]);
      ret += double(sample.counts[type_index]) *
             quadrature(temp_type, type_vec, sample.sample_size) /
             double(sample.sample_size);
    }
    return ret;
  }

  double quadrature(const std::vector<int> &a, const std::vector<int> &b,
                    const int sample_size) const {
    double ret = 0.0, temp;
    for (unsigned int i = 0; i < quad_points.size(); i++) {
      temp = 1.0;
      for (unsigned int j = 0; j < a.size(); j++) {
        temp *= M[i * row_length + sample_size - 1].getEntry(a[j], b[j]);
      }
      ret +=
          temp * quad_points[i] /
          pow(((double)(quad_points.size()) + 1.0) *
                  (-pow(quad_points[i], 5.0) + 25.0 * pow(quad_points[i], 4.0) -
                   200.0 * pow(quad_points[i], 3.0) +
                   600.0 * pow(quad_points[i], 2.0) - 600.0 * quad_points[i] +
                   120.0) /
                  120.0,
              2.0);
    }
    return ret;
  }

  const double theta;
  std::vector<double> quad_points;
  const int row_length;
  std::vector<Matrix> M;
};

#endif
