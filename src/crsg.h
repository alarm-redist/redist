#ifndef CRSG_H
#define CRSG_H

#include <Rcpp.h>
#include "shatter_search.h"
#include "distance_helpers.h"
using namespace Rcpp;

List crsg(List adj_list,
          NumericVector population,
          NumericVector area,
          NumericVector x_center,
          NumericVector y_center,
          int Ndistrict,
          double target_pop,
          double thresh,
          int maxiter);

#endif