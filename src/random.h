#ifndef RANDOM_H
#define RANDOM_H

#include <RcppArmadillo.h>
#include <vector>
#include <cstdint>
#include <random>

using namespace arma;

/*
 * Set RNG seed
 */
void seed_rng(int seed);

/*
 * Generate a uniform random integer in [0, max). Very slightly biased.
 */
int r_int(uint32_t max);

/*
 * Generate a uniform random double in [0, 1). Very slightly biased.
 */
double r_unif();

/*
 * Generate a random integer in [0, max) according to weights.
 */
int r_int_wgt(int max, vec cum_wgts);

/*
 * Generate a random integer within a stratum with some probability p
 */
int r_int_mixstrat(int max, int stratum, double p, vec cum_wgts);

/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
// [[Rcpp::export]]
arma::ivec resample_lowvar(arma::vec wgts);


#endif
