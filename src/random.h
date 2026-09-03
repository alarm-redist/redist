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
 * Get / restore the internal RNG state, for checkpoint / resume.
 *
 * The state is packed into a length-6 integer vector by bit-casting:
 *   [0:1] = state_sr    (one uint64_t)
 *   [2:5] = state_xo[4] (four uint32_t)
 * The integers are opaque (they may read as negative or as NA); only their
 * bit patterns matter, and they round-trip exactly through save / load.
 */
Rcpp::IntegerVector rng_state_get();
void rng_state_set(const Rcpp::IntegerVector& state);

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
