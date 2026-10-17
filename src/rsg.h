#ifndef RSG_H
#define RSG_H

#include <Rcpp.h>

Rcpp::List rsg(Rcpp::List adj_list, Rcpp::NumericVector population, int Ndistrict,
               double target_pop, double thresh, int maxiter);

#endif
