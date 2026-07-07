#ifndef UTILS_H
#define UTILS_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include "redist_types.h"
#include <RcppArmadillo.h>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]


#include "random.h"



/*
 * Make a progress bar configuration with format string `fmt`
 */
Rcpp::List cli_config(bool clear = false,
                const char *fmt = "{cli::pb_bar} {cli::pb_percent} | ETA:{cli::pb_eta}");

#endif
