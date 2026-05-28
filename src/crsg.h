#ifndef CRSG_H
#define CRSG_H

#include "constraint_calc_helper.h"
#include "distance_helpers.h"
#include "make_swaps_helper.h"
#include "redist_types.h"
#include "shatter_search.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

List crsg(List adj_list, NumericVector population, NumericVector area, NumericVector x_center,
          NumericVector y_center, int Ndistrict, double target_pop, double thresh, int maxiter);

#endif