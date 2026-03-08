#ifndef MMSS_H
#define MMSS_H

#include "merge_split.h"

/*
 * Main entry point for Multi-Merge-Split MCMC.
 *
 * Sample `N` redistricting plans on map `g`, merging and re-splitting
 * `l_merge` districts at each step.  When l_merge=2 this reduces to
 * standard merge-split.
 */
Rcpp::List mmss_plans(int N, List l, const arma::uvec init,
                     const arma::uvec &counties, const arma::uvec &pop,
                     int n_distr, double target, double lower, double upper,
                     double rho, List constraints, List control,
                     int k, int thin, int l_merge, int verbosity);

/*
 * Select `l` connected districts from the district graph via random BFS.
 * Returns a vector of 1-indexed district labels.
 * `log_prob` is set to log(probability of this selection) for MH correction.
 */
std::vector<int> select_l_districts(int n_distr, const Graph &dist_g,
                                    int l, double &log_prob);

#endif
