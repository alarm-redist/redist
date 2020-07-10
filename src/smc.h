#ifndef SMC_H
#define SMC_H

#include <math.h>
#include "smc_base.h"
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "kirchhoff.h"

/*
 * Main entry point.
 *
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is within `tol`
 */
// [[Rcpp::export]]
IntegerMatrix smc_plans(int N, List g, const IntegerVector &counties,
                        const IntegerVector &pop, int n_distr, double tol,
                        double gamma, NumericVector &log_prob, double thresh,
                        double alpha, int infl=5, int verbosity=1);

/*
 * Split off a piece from each map in `districts`, keeping deviation within `tol`
 */
void split_maps(List g, const IntegerVector &counties, Multigraph &cg,
                const IntegerVector &pop, IntegerMatrix &districts,
                NumericVector &lp, NumericVector &pop_left, int N, int n_distr,
                int dist_ctr, double distr_pop, double tol, double gamma,
                int k, int verbosity);

/*
 * Resample partially-drawn maps according to their weights.
 */
void resample_maps(int N_drawn, int N_sample, double alpha, IntegerMatrix &districts,
                   NumericVector &lp, NumericVector &pop_left);

/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map(List g, const IntegerVector &counties, Multigraph &cg,
                 IntegerMatrix::Column districts, int dist_ctr,
                 const IntegerVector &pop, double total_pop,
                 double &lower, double upper, double target, int k);

/*
 * Cut spanning subtree into two pieces of roughly equal population
 */
double cut_districts(Tree &ust, int k, int root, IntegerMatrix::Column districts,
                     int dist_ctr, const IntegerVector &pop, double total_pop,
                     double lower, double upper, double target);

/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_parameters(List g, int &k, double &prob, int N_adapt, int valid,
                      const NumericVector &lp, double thresh,
                      double tol, const IntegerMatrix &districts,
                      const IntegerVector &counties, Multigraph &cg,
                      const IntegerVector &pop, const NumericVector &pop_left,
                      double target);

/*
 * Partition `x` and its indices `idxs` between `right` and `left` by `pivot`
 */
// TESTED
void partition_vec(std::vector<double> &x, std::vector<int> &idxs, int left,
                   int right, int &pivot);

/*
 * Get the index of the k-th smallest element of x
 */
// TESTED
int select_k(std::vector<double> x, int k);

#endif
