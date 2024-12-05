#include "infsites.hh"
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <vector>

double log_sum_exp_moment(const std::vector<double> &v, const int j,
                          const int u) {
  double m = v[0];
  for (int i = 1; i < u; i++) {
    if (v[i] > m) {
      m = v[i];
    }
  }
  double ret = 0;
  for (int i = 0; i < u; i++) {
    ret += exp(j * (v[i] - m));
  }
  return j * m + log(ret) - log(u);
}

double log_chan_lai(const std::vector<double> &v, const std::vector<int> &eve,
                    const int min_particle_no) {
  double mu = log_sum_exp_moment(v, 1, int(v.size()));
  std::vector<double> maxima(min_particle_no, -DBL_MAX);
  std::vector<int> counts(min_particle_no, 0);
  for (unsigned int i = 0; i < v.size(); i++) {
    maxima[eve[i]] = fmax(maxima[eve[i]], v[i]);
    counts[eve[i]]++;
  }
  std::vector<double> r(min_particle_no, 0);
  for (unsigned int i = 0; i < v.size(); i++) {
    r[eve[i]] += exp(v[i] - maxima[eve[i]]);
  }
  double r_max = -DBL_MAX;
  double local_mu;
  for (int i = 0; i < min_particle_no; i++) {
    if (counts[i] > 0) {
      local_mu = mu + log(counts[i]);
      r[i] = maxima[i] + log(r[i]);
      if (r[i] > local_mu) {
        r[i] = 2 * (r[i] + log(1 - exp(local_mu - r[i])));
      } else {
        r[i] = 2 * (local_mu + log(1 - exp(r[i] - local_mu)));
      }
      r_max = fmax(r_max, r[i]);
    }
  }
  double ret = 0;
  for (int i = 0; i < min_particle_no; i++) {
    if (counts[i] > 0) {
      ret += exp(r[i] - r_max);
    }
  }
  ret = r_max + log(ret) - log(v.size());
  return ret;
}

double ess(const std::vector<double> &w) {
  double ret = w.size() * exp(2 * log_sum_exp_moment(w, 1, int(w.size())) -
                              log_sum_exp_moment(w, 2, int(w.size())));
  return ret;
}

void resample(std::vector<infSites> &samples, std::vector<double> &w,
              std::vector<int> &eve, gsl_rng *gen) {
  int n = samples.size();
  double u = gsl_rng_uniform(gen);
  double log_w_sum = log_sum_exp_moment(w, 1, int(w.size())) + log(n);
  int parent = 0;
  double ub = exp(w[0] - log_w_sum);
  std::vector<infSites> new_samples(w.size(), samples[0]);
  std::vector<int> new_eve(w.size(), 0);
  for (int i = 0; i < n; i++) {
    while ((i + u) / n > ub) {
      parent++;
      ub += exp(w[parent] - log_w_sum);
    }
    new_samples[i].assign(samples[parent]);
    new_eve[i] = eve[parent];
  }
  for (int i = 0; i < n; i++) {
    samples[i].assign(new_samples[i]);
    w[i] = log_w_sum - log(n);
    eve[i] = new_eve[i];
  }
  return;
}

double step(infSites &infsites, gsl_rng *gen, const double theta) {
  int nc = 0;
  for (unsigned int i = 0; i < infsites.row_counts.size(); i++) {
    if (infsites.row_counts[i] > 1) {
      nc += infsites.row_counts[i];
    } else if (infsites.find_singleton(i) == 1) {
      nc += 1;
    }
  }
  double coin = gsl_rng_uniform(gen);
  int ind = -1;
  double ub = 0;
  do {
    ind++;
    if (infsites.row_counts[ind] > 1) {
      ub += double(infsites.row_counts[ind]) / nc;
    } else if (infsites.find_singleton(ind) == 1) {
      ub += 1.0 / nc;
    }
  } while (coin > ub);
  double ret = 0;
  if (infsites.row_counts[ind] > 1) {
    ret =
        ((infsites.row_counts[ind] - 1) / (infsites.sample_size - 1 + theta)) /
        (double(infsites.row_counts[ind]) / nc);
    infsites.row_counts[ind]--;
    for (int j = 0; j < infsites.data.ncol; j++) {
      infsites.col_counts[j] -= infsites.data.getEntry(ind, j);
    }
    infsites.sample_size--;
  } else {
    int j = infsites.sample_singleton(ind, gen);
    j = infsites.remove_mutation(ind, j);
    ret = nc * (double(infsites.row_counts[j]) / infsites.sample_size) *
          (theta / (infsites.sample_size - 1 + theta));
  }
  return log(ret);
}

void smc(const infSites &infsites, const double theta,
         const int min_particle_no, const int max_particle_no, const int output,
         const double resampling_frac, gsl_rng *gen) {
  infSites temp_sample(infsites);
  std::vector<infSites> temp_samples(max_particle_no, temp_sample);
  std::vector<double> likelihoods(max_particle_no, 0.0);
  int boundary = infsites.sample_size - 1;
  int res;
  double c = 0.1;
  int target = floor(pow(infsites.sample_size,
                         pow(c, 1 / (theta * log(infsites.sample_size)))));
  std::vector<int> eve(max_particle_no, 0);
  for (int i = 1; i < min_particle_no; i++) {
    eve[i] = i;
  }
  int particle_no = min_particle_no;
  for (int k = 0; k < 2; k++) {
    while (boundary > target) {
      res = 0;
      if (ess(likelihoods) < particle_no * resampling_frac) {
        res = 1;
        resample(temp_samples, likelihoods, eve, gen);
      }
      for (int i = 0; i < particle_no; i++) {
        while (temp_samples[i].sample_size > boundary) {
          likelihoods[i] += step(temp_samples[i], gen, theta);
        }
      }
      if (output == 1) {
        std::cout << double(boundary) / infsites.sample_size << " "
                  << exp(log_sum_exp_moment(likelihoods, 2, particle_no) -
                         2 * log_sum_exp_moment(likelihoods, 1, particle_no))
                  << " " << exp(log_sum_exp_moment(likelihoods, 1, particle_no))
                  << " " << res << std::endl;
      }
      boundary--;
    }
    if (k == 0) {
      int ratio = max_particle_no / min_particle_no;
      int ind = max_particle_no - 1;
      for (int i = min_particle_no - 1; i > -1; i--) {
        for (int j = 0; j < ratio; j++) {
          temp_samples[ind].assign(temp_samples[i]);
          likelihoods[ind] = likelihoods[i];
          eve[ind] = eve[i];
          ind--;
        }
      }
      target = 0;
      particle_no = max_particle_no;
    }
  }
  res = 0;
  if (output == 1) {
    std::cout << double(boundary) / infsites.sample_size << " "
              << exp(log_sum_exp_moment(likelihoods, 2, max_particle_no) -
                     2 * log_sum_exp_moment(likelihoods, 1, max_particle_no))
              << " " << exp(log_sum_exp_moment(likelihoods, 1, max_particle_no))
              << " " << res << std::endl;
  }
  if (output == 0) {
    std::cout << theta << " "
              << log_sum_exp_moment(likelihoods, 1, max_particle_no) << " "
              << log_chan_lai(likelihoods, eve, min_particle_no) << std::endl;
  }
  return;
}

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cout << "Call " << argv[0]
              << " <input_file> <theta> <min_particles> <max_particles> "
                 "<output> <resampling_fraction>"
              << std::endl;
    std::cout << "<output> = 0 => log mean and standard error of weights"
              << std::endl;
    std::cout << "<output> = 1 => running relative second moment" << std::endl;
    return 1;
  }
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  double theta = atof(argv[2]);
  infSites infsites(argv[1], theta);
  int min_particles = atoi(argv[3]);
  int max_particles = atoi(argv[4]);
  int output = atoi(argv[5]);
  const double resampling_frac = atof(argv[6]);

  smc(infsites, theta, min_particles, max_particles, output, resampling_frac,
      gen);

  gsl_rng_free(gen);
  return 0;
}
