#ifndef FS_SAMPLE
#define FS_SAMPLE

#include "finite_sites.hh"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <vector>

class FSSample {
public:
  FSSample(const int count, FiniteSites &mutation, gsl_rng *generator)
      : counts(1, count), types(1), sample_size(count), sample_types(1) {
    std::vector<int> type(mutation.site_no);
    for (int i = 0; i < mutation.site_no; i++) {
      type[i] = floor(gsl_rng_uniform(generator) * 2.0);
    }
    types[0] = mutation.vec_to_int(type);
  }

  FSSample(char *filename)
      : counts(0), types(0), sample_size(0), sample_types(0) {
    std::ifstream dataFile;
    dataFile.open(filename);
    if (!dataFile.good()) {
      std::cout << "Failed to open data file. Aborting" << std::endl;
      abort();
    }
    char *stopstring;
    int base = 10;
    std::string line, token;
    std::stringstream iss;
    getline(dataFile, line);
    while (!dataFile.eof()) {
      iss.str(std::string());
      iss.clear();
      iss << line;
      getline(iss, token, ' ');
      counts.push_back(atoi(token.c_str()));
      getline(iss, token, ' ');
      types.push_back(strtoull(token.c_str(), &stopstring, base));
      sample_size += counts.back();
      sample_types++;
      getline(dataFile, line);
    }
  }

  FSSample(const FSSample &other)
      : counts(other.counts), types(other.types),
        sample_size(other.sample_size), sample_types(other.sample_types) {}

  void add(const unsigned long long int type, const int no) {
    for (int i = 0; i < sample_types; i++) {
      if (types[i] == type) {
        counts[i] += no;
        sample_size += no;
        return;
      }
    }
    counts.push_back(no);
    sample_size += no;
    types.push_back(type);
    sample_types++;
    return;
  }

  void kill(const unsigned long long int type, const int no) {
    for (unsigned int i = 0; i < types.size(); i++) {
      if (types[i] == type) {
        counts[i] -= no;
        sample_size -= no;
        if (counts[i] == 0) {
          counts.erase(counts.begin() + i);
          types.erase(types.begin() + i);
          sample_types--;
        }
        return;
      }
    }
    std::cout << "Type " << type << " not found in sample. Aborting."
              << std::endl;
    abort();
  }

  void print(char *fileName) const {
    std::ofstream outputFile;
    outputFile.open(fileName);
    if (!outputFile.good()) {
      std::cout << "Could not open output file " << fileName << std::endl;
      return;
    }
    for (int i = 0; i < sample_types; i++) {
      outputFile << counts[i] << " " << types[i] << std::endl;
    }

    outputFile.close();
    return;
  }

  int abundance(const unsigned long long int type) const {
    for (int i = 0; i < sample_types; i++) {
      if (types[i] == type) {
        return counts[i];
      }
    }
    return 0;
  }

  unsigned long long int sample_parent(gsl_rng *generator) const {
    int index = 0;
    double ub = double(counts[0]) / sample_size;
    double coin = gsl_rng_uniform(generator);
    while (coin > ub) {
      index++;
      ub += double(counts[index]) / sample_size;
    }
    return types[index];
  }

  void assign(const FSSample &other) {
    counts = other.counts;
    types = other.types;
    sample_size = other.sample_size;
    sample_types = other.sample_types;
  }

  std::vector<int> counts;
  std::vector<unsigned long long int> types;
  int sample_size, sample_types;
};

#endif
