#ifndef FIN_SITES
#define FIN_SITES

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>

class FiniteSites {
public:
  FiniteSites(const int s_n) : site_no(s_n) {}

  std::vector<int> int_to_vec(const unsigned long long int type) const {
    std::vector<int> ret(site_no);
    unsigned long long int total = 0;
    for (int i = 0; i < site_no; i++) {
      ret[i] = floor((type - total) / pow(2, site_no - i - 1));
      if (ret[i] > 0) {
        total += ret[i] * pow(2, site_no - i - 1);
      }
    }
    return ret;
  }

  unsigned long long int vec_to_int(const std::vector<int> &sites) const {
    unsigned long long int ret = 0;
    for (int i = 0; i < site_no; i++) {
      ret += sites[i] * pow(2, site_no - i - 1);
    }
    return ret;
  }

  std::vector<unsigned long long int>
  neighbours(const unsigned long long int parent_type) const {
    std::vector<unsigned long long int> ret(site_no);
    std::vector<int> parent_vec = int_to_vec(parent_type);
    for (int i = 0; i < site_no; i++) {
      parent_vec[i] = 1 - parent_vec[i];
      ret[i] = vec_to_int(parent_vec);
      parent_vec[i] = 1 - parent_vec[i];
    }
    return ret;
  }

  unsigned long long int mutate(const unsigned long long int parent_type,
                                gsl_rng *generator) const {
    std::vector<int> parent_vec = int_to_vec(parent_type);
    int mutant = floor(gsl_rng_uniform(generator) * site_no);
    parent_vec[mutant] = 1 - parent_vec[mutant];
    return vec_to_int(parent_vec);
  }

  unsigned long long int random_type(gsl_rng *generator) const {
    std::vector<int> type_vec(site_no);
    for (unsigned int i = 0; i < type_vec.size(); i++) {
      type_vec[i] = floor(2.0 * gsl_rng_uniform(generator));
    }
    return vec_to_int(type_vec);
  }

  double stationary_law(const unsigned long long int type) const {
    return 1.0 / pow(2.0, (double)(site_no));
  }

  int site_no;
};

#endif
