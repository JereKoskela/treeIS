#include "csd.hh"
#include "finite_sites.hh"
#include "fs_sample.hh"
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <sstream>
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

void resample(std::vector<FSSample> &samples, std::vector<double> &w,
              std::vector<int> &eve, gsl_rng *gen) {
  int n = samples.size();
  double u = gsl_rng_uniform(gen);
  double log_w_sum = log_sum_exp_moment(w, 1, int(w.size())) + log(n);
  int parent = 0;
  double ub = exp(w[0] - log_w_sum);
  std::vector<FSSample> new_samples(w.size(), samples[0]);
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

double step(FSSample &sample, gsl_rng *generator, const FiniteSites &mutation,
            const CSD &acsd) {
  int parent_index = 0;
  double coin = gsl_rng_uniform(generator);
  double ub = double(sample.counts[0]) / sample.sample_size;
  while (coin > ub) {
    parent_index++;
    ub += double(sample.counts[parent_index]) / sample.sample_size;
  }
  double coal_prob = sample.counts[parent_index] - 1.0;
  std::vector<double> mut_probs(mutation.site_no);
  double mutation_prob = 0.0;
  FSSample temp_sample(sample);
  temp_sample.kill(temp_sample.types[parent_index], 1);
  std::vector<unsigned long long int> neighbours =
      mutation.neighbours(sample.types[parent_index]);
  for (int i = 0; i < mutation.site_no; i++) {
    mut_probs[i] = acsd.theta *
                   acsd.pi_hat(neighbours[i], temp_sample, mutation) /
                   mutation.site_no;
    mutation_prob += mut_probs[i];
  }
  double num, denom;
  if (gsl_rng_uniform(generator) <
      mutation_prob / (mutation_prob + coal_prob)) {
    int mutant_index = 0;
    coin = gsl_rng_uniform(generator);
    ub = mut_probs[0] / mutation_prob;
    while (coin > ub) {
      mutant_index++;
      ub += mut_probs[mutant_index] / mutation_prob;
    }
    denom = sample.counts[parent_index] * acsd.theta *
            acsd.pi_hat(neighbours[mutant_index], temp_sample, mutation) /
            (acsd.pi_hat(sample.types[parent_index], temp_sample, mutation) *
             sample.sample_size * (sample.sample_size + acsd.theta - 1) *
             mutation.site_no);
    num = acsd.theta * (sample.abundance(neighbours[mutant_index]) + 1) /
          (sample.sample_size * (sample.sample_size - 1 + acsd.theta) *
           mutation.site_no);
    sample.kill(sample.types[parent_index], 1);
    sample.add(neighbours[mutant_index], 1);
  } else {
    denom = 2 * gsl_sf_choose(sample.counts[parent_index], 2) /
            (acsd.pi_hat(sample.types[parent_index], temp_sample, mutation) *
             sample.sample_size * (sample.sample_size + acsd.theta - 1));
    num = (sample.counts[parent_index] - 1) /
          (sample.sample_size - 1 + acsd.theta);
    sample.kill(sample.types[parent_index], 1);
  }
  return log(num) - log(denom);
}

void smc(const FSSample &sample, const double theta,
         const FiniteSites &mutation, const int min_particle_no,
         const int max_particle_no, const int output,
         const double resampling_frac, gsl_rng *generator) {
  CSD acsd(theta, sample.sample_size, mutation);
  FSSample temp_sample(sample);
  std::vector<FSSample> temp_samples(max_particle_no, temp_sample);
  std::vector<double> likelihoods(max_particle_no, 0.0);
  std::vector<int> eve(max_particle_no, 0);
  for (int i = 1; i < min_particle_no; i++) {
    eve[i] = i;
  }
  int boundary = sample.sample_size - 1;
  int res;
  double c = 0.1;
  int target = floor(
      pow(sample.sample_size, pow(c, 1 / (theta * log(sample.sample_size)))));
  int particle_no = min_particle_no;
  for (int k = 0; k < 2; k++) {
    while (boundary > target) {
      res = 0;
      if (ess(likelihoods) < particle_no * resampling_frac) {
        res = 1;
        resample(temp_samples, likelihoods, eve, generator);
      }
      for (int i = 0; i < particle_no; i++) {
        while (temp_samples[i].sample_size > boundary) {
          likelihoods[i] += step(temp_samples[i], generator, mutation, acsd);
        }
      }
      if (output == 1) {
        std::cout << double(boundary) / sample.sample_size << " "
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
  for (int i = 0; i < max_particle_no; i++) {
    likelihoods[i] += log(mutation.stationary_law(temp_samples[i].types[0]));
  }
  if (output == 1) {
    std::cout << double(boundary) / sample.sample_size << " "
              << exp(log_sum_exp_moment(likelihoods, 2, max_particle_no) -
                     2 * log_sum_exp_moment(likelihoods, 1, max_particle_no))
              << " " << exp(log_sum_exp_moment(likelihoods, 1, max_particle_no))
              << " " << res << std::endl;
  }
  if (output == 0) {
    std::cout << acsd.theta * 2 << " "
              << log_sum_exp_moment(likelihoods, 1, max_particle_no) << " "
              << log_chan_lai(likelihoods, eve, min_particle_no) << std::endl;
  }
  return;
}

int main(int argc, char **argv) {
  if (argc != 8) {
    std::cout << "Call " << argv[0]
              << " <sample_size> <site_no> <theta> <min_particles> "
                 "<max_particles> <output> <resampling_fraction>"
              << std::endl;
    std::cout << "<output> = 0 => log mean and standard error of weights"
              << std::endl;
    std::cout << "<output> = 1 => running relative second moment" << std::endl;
    return 1;
  }
  const int sample_size = atoi(argv[1]);
  const int site_no = atoi(argv[2]);
  const double theta = atof(argv[3]);
  const int min_particle_no = atoi(argv[4]);
  const int max_particle_no = atoi(argv[5]);
  const int output = atoi(argv[6]);
  const double resampling_frac = atof(argv[7]);
  const FiniteSites mutation(site_no);
  std::stringstream fileStream;
  fileStream << "../Samples/finitealleles-" << sample_size << ".dat";
  std::string fileName = fileStream.str();
  FSSample sample(const_cast<char *>(fileName.c_str()));
  gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(gen, time(NULL));
  smc(sample, theta, mutation, min_particle_no, max_particle_no, output,
      resampling_frac, gen);
  gsl_rng_free(gen);

  return 1;
}
