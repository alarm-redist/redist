#ifndef SMC_BASE_H
#define SMC_BASE_H

#include <stdlib.h>
#include <vector>
#include <random>
#include <iostream>
#include <limits>
#include <RcppArmadillo.h>
#include "redist_types.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
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
 * Get the index of the k-th smallest element of x
 */
// TESTED
int select_k(std::vector<double> x, int k);

#endif
