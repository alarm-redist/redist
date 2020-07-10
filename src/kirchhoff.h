#ifndef KIRCHHOFF_H
#define KIRCHHOFF_H

#include "smc_base.h"

/*
 * Compute the log number of spanning trees which could generate a given set of maps.
 * `districts` should have each column be a map
 */
// TESTED
// [[Rcpp::export]]
NumericVector log_st_map(const List &g, const IntegerMatrix &districts,
                         const IntegerVector &counties, int n_distr);

/*
 * Compute the log number of spanning trees for `district` intersect `county`
 */
// TESTED
double log_st_distr(const List &g, const IntegerMatrix &districts,
                    const IntegerVector &counties, int idx,
                    int district, int county);

/*
 * Compute the log number of spanning trees for the contracted graph
 */
// TESTED
double log_st_contr(const List &g, const IntegerMatrix &districts,
                    const IntegerVector &counties, int n_cty,
                    int idx, int district);

#endif
