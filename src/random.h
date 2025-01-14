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
 * Generate a random integer in [0, max) according to cumulative normalized weights.
 */
int r_int_wgt(int max, vec cum_wgts);

//' Generate a random index of `unnormalized_wgts` with probability proportional to its weight
//'
//' Takes a vector of strictly positive weights and returns an index with probability 
//' proportional to its weight. In other words, it selects index `i` with probability
//' proporitional to `unnormalized_wgts[i]` 
//' (or exactly `unnormalized_wgts[i]/sum(unnormalized_wgts)`). This does not support
//' inputs where some of the weights are zero. This has positive probability of 
//' returning indices that have weight zero. 
//'
//'
//' @param unnormalized_wgts An arma vector of positive numbers
//'
//' @details no Modifications to inputs made
//'
//' @returns An integer in [0, `unnormalized_wgts.size()`)
//'
//' @keyword internal
//' @noRd
int r_int_unnormalized_wgt(const vec &unnormalized_wgts);

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
