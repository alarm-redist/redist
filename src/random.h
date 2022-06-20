#ifndef RANDOM_H
#define RANDOM_H

#include <RcppArmadillo.h>
#include <vector>
#include <random>

using namespace arma;

// INITIALIZE MT RNG
extern std::random_device rd;
extern std::mt19937 generator;
extern std::uniform_real_distribution<double> unif;

/*
 * Generate a uniform random integer in [0, max).
 */
int rint(int max);

/*
 * Generate a random integer in [0, max) according to weights.
 */
int rint(int max, vec cum_wgts);

/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
// [[Rcpp::export]]
arma::ivec resample_lowvar(arma::vec wgts);


#endif
