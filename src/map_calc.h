#include "smc_base.h"
#include "tree_op.h"

#ifndef MAP_CALC_H
#define MAP_CALC_H

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_boundary(const List &g, const IntegerMatrix::Column &districts,
                    int distr_root, int distr_other);

/*
 * Compute the deviation from the equal population constraint.
 */
NumericMatrix pop_dev(const IntegerMatrix &districts,
                      const NumericVector &pop, int n_distr);

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// [[Rcpp::export]]
NumericVector max_dev(const IntegerMatrix &districts,
                      const NumericVector &pop, int n_distr);

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
void tree_dev(Tree &ust, int root, vec res,
              const IntegerVector &pop, double total_pop, double target);

/*
 * Compute the isoperimetric quotient for maps, given the area of each precinct
 * and a matrix of pariwise boundary lengths (with external boundary lengths
 * along the diagonal)
 */
// [[Rcpp::export]]
NumericMatrix calc_isoper_quo(const IntegerMatrix &districts,
                              const NumericVector &areas,
                              const NumericMatrix &perims, int n_distr);


#endif