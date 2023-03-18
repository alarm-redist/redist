#ifndef SMC_BASE_H
#define SMC_BASE_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include <vector>
#include <limits>
#include <RcppArmadillo.h>
#include "redist_types.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

#include "random.h"

/*
 * Get the index of the k-th smallest element of x
 */
// TESTED
int select_k(std::vector<double> x, int k);

/*
 * Make a progress bar configuration with format string `fmt`
 */
List cli_config(bool clear = false,
                const char * fmt = "{cli::pb_bar} {cli::pb_percent} | ETA:{cli::pb_eta}");

#endif
