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

double log_sum_exp_moment(const std::vector<double> &v, const int j) {
  double m = v[0];
  for (unsigned int i = 1; i < v.size(); i++) {
    if (v[i] > m) {
      m = v[i];
    }
  }
  double ret = 0;
  for (unsigned int i = 0; i < v.size(); i++) {
    ret += exp(j * (v[i] - m));
  }
  return j * m + log(ret) - log(v.size());
}

double ess(const std::vector<double> &w) {
  double ret =
      w.size() * exp(2 * log_sum_exp_moment(w, 1) - log_sum_exp_moment(w, 2));
  return ret;
}

void resample(std::vector<FSSample> &samples, std::vector<double> &w,
              std::vector<int> &b, gsl_rng *gen) {
  int n = samples.size();
  double u = gsl_rng_uniform(gen);
  double log_w_sum = log_sum_exp_moment(w, 1) + log(n);
  int parent = 0;
  double ub = exp(w[0] - log_w_sum);
  std::vector<FSSample> new_samples(w.size(), samples[0]);
  std::vector<int> new_b(b.size(), 0);
  for (int i = 0; i < n; i++) {
    while ((i + u) / n > ub) {
      parent++;
      ub += exp(w[parent] - log_w_sum);
    }
    new_samples[i].assign(samples[parent]);
    new_b[i] = b[parent];
  }
  for (int i = 0; i < n; i++) {
    samples[i].assign(new_samples[i]);
    b[i] = new_b[i];
    w[i] = log_w_sum - log(n);
  }
  return;
}

double step(FSSample &sample, gsl_rng *generator, const FiniteSites &mutation,
            const double theta, int &b) {
  std::vector<double> coal_probs(sample.sample_types, 0);
  std::vector<double> mut_probs(sample.sample_types, 0);
  std::vector<unsigned long long int> neighbours(mutation.site_no);
  double total_prob = 0;
  for (int i = 0; i < sample.sample_types; i++) {
    if (sample.counts[i] > 1) {
      coal_probs[i] = (sample.counts[i] - 1) / (sample.sample_size - 1 + theta);
    }
    neighbours = mutation.neighbours(sample.types[i]);
    for (int j = 0; j < mutation.site_no; j++) {
      mut_probs[i] += (sample.abundance(neighbours[j]) + 1) * theta /
                      (sample.sample_size * (sample.sample_size + theta - 1) *
                       mutation.site_no);
    }
    total_prob += coal_probs[i] + mut_probs[i];
  }
  int parent_index = 0;
  double coin = gsl_rng_uniform(generator);
  double ub = (coal_probs[0] + mut_probs[0]) / total_prob;
  while (coin > ub) {
    parent_index++;
    ub += (coal_probs[parent_index] + mut_probs[parent_index]) / total_prob;
  }
  if (gsl_rng_uniform(generator) <
      mut_probs[parent_index] /
          (mut_probs[parent_index] + coal_probs[parent_index])) {
    b++;
    neighbours = mutation.neighbours(sample.types[parent_index]);
    coin = gsl_rng_uniform(generator);
    int mutant_index = 0;
    ub = (sample.abundance(neighbours[0]) + 1) * theta /
         (sample.sample_size * (sample.sample_size + theta - 1) *
          mutation.site_no * mut_probs[parent_index]);
    while (coin > ub) {
      mutant_index++;
      ub += (sample.abundance(neighbours[mutant_index]) + 1) * theta /
            (sample.sample_size * (sample.sample_size + theta - 1) *
             mutation.site_no * mut_probs[parent_index]);
    }
    sample.kill(sample.types[parent_index], 1);
    sample.add(neighbours[mutant_index], 1);
  } else {
    sample.kill(sample.types[parent_index], 1);
  }
  return log(total_prob);
}

void smc(const FSSample &sample, const double theta,
         const FiniteSites &mutation, const int particle_no, const int bound,
         gsl_rng *generator) {
  FSSample temp_sample(sample);
  std::vector<FSSample> temp_samples(particle_no, temp_sample);
  std::vector<double> likelihoods(particle_no, 0.0);
  std::vector<int> b(particle_no, 0);
  int boundary = sample.sample_size - 1;
  int res;
  while (boundary > 0) {
    res = 0;
    if (ess(likelihoods) < 1) {
      res = 1;
      resample(temp_samples, likelihoods, b, generator);
    }
    for (int i = 0; i < particle_no; i++) {
      while (temp_samples[i].sample_size > boundary && b[i] <= bound) {
        likelihoods[i] +=
            step(temp_samples[i], generator, mutation, theta, b[i]);
      }
      if (b[i] > bound) {
        likelihoods[i] = -DBL_MAX / 2;
      }
    }
    std::cout << double(boundary) / sample.sample_size << " "
              << exp(log_sum_exp_moment(likelihoods, 2) -
                     2 * log_sum_exp_moment(likelihoods, 1))
              << " " << exp(log_sum_exp_moment(likelihoods, 1)) << " " << res
              << std::endl;
    boundary--;
  }
  for (int i = 0; i < particle_no; i++) {
    likelihoods[i] += log(mutation.stationary_law(temp_samples[i].types[0]));
  }
  std::cout << double(boundary) / sample.sample_size << " "
            << exp(log_sum_exp_moment(likelihoods, 2) -
                   2 * log_sum_exp_moment(likelihoods, 1))
            << " " << exp(log_sum_exp_moment(likelihoods, 1)) << " " << res
            << std::endl;
  return;
}

int main(int argc, char **argv) {
  if (argc != 6) {
    std::cout << "Call " << argv[0]
              << " <sample_size> <particle_no> <site_no> <theta> <bound>"
              << std::endl;
    return 1;
  }
  const int sample_size = atoi(argv[1]);
  const int particle_no = atoi(argv[2]);
  const int site_no = atoi(argv[3]);
  const double theta = atof(argv[4]);
  const int b = atoi(argv[5]);
  const FiniteSites mutation(site_no);
  std::stringstream fileStream;
  fileStream << "../Samples/finitealleles-" << sample_size << ".dat";
  std::string fileName = fileStream.str();
  FSSample sample(const_cast<char *>(fileName.c_str()));
  gsl_rng *mersenne_twister = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(mersenne_twister, time(NULL));
  smc(sample, theta, mutation, particle_no, b, mersenne_twister);
  gsl_rng_free(mersenne_twister);

  return 1;
}
