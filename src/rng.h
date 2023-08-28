/**** Random number generator wraper class ****/

#ifndef RNG_H
#define RNG_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include <Rmath.h>

class RNG {

private:
  Rcpp::RNGScope rng_scope;

public:
  RNG() {}
  ~RNG() {}

  double runif(double min = 0.0, double max = 1.0) {return R::runif(min, max);}

  double rnorm(double mu = 0.0, double sigma = 1.0) {return R::rnorm(mu, sigma);}

  // gamma distribution in shape-rate parameterization
  double rgamma(double shape, double rate) {return R::rgamma(shape, 1.0 / rate);}

  // discrete uniform distribution on nonnegative integers {min, min+1, ..., max}
  unsigned int rdunif(unsigned int min, unsigned int max) {
    double sample = std::ceil(R::runif(static_cast<double>(min) - 1.0, static_cast<double>(max)));
    return static_cast<unsigned int>(sample);
  }

  // discrete distribution on {0, 1, ..., probs.size()-1}
  // @param probs: probs[i] is the probability for i
  unsigned int rdiscrete(std::vector<double> probs) {
    // get cdf
    for (unsigned int i = 0; i < probs.size(); ++i) {
      if (probs[i] < 0)
        throw std::invalid_argument("Probabality cannot be negative");
      if (i > 0) probs[i] += probs[i - 1];
    }

    double u = runif(0, probs.back());
    // binary search
    auto it = std::lower_bound(probs.begin(), probs.end(), u);
    return std::distance(probs.begin(), it);
  }

  // other distributions
};

#endif
