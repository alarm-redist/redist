#ifndef MERGESPLIT_H
#define MERGESPLIT_H

#include "smc_base.h"

#include <string>
#include <cli/progress.h>

// [[Rcpp::depends(redistmetrics)]]

#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include <kirchhoff_inline.h>
#include "mcmc_gibbs.h"

/*
 * Main entry point.
 *
 * USING MCMMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
Rcpp::List ms_plans(int N, List l, const arma::uvec init, const arma::uvec &counties,
                    const arma::uvec &pop, int n_distr, double target, double lower,
                    double upper, double rho, List constraints, List control,
                    int k, int thin, int verbosity);


/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map_ms(const Graph &g, const uvec &counties, Multigraph &cg,
                    subview_col<uword> districts, int distr_1, int distr_2,
                    const uvec &pop, double lower, double upper, double target,
                    int k);

/*
 * Cut district into two pieces of roughly equal population
 */
// TESTED
bool cut_districts_ms(Tree &ust, int k, int root, subview_col<uword> &districts,
                      int distr_1, int distr_2, const uvec &pop, double total_pop,
                      double lower, double upper, double target);

/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_ms_parameters(const Graph &g, int n_distr, int &k, double thresh,
                         double tol, const uvec &plan, const uvec &counties,
                         Multigraph &cg, const uvec &pop, double target);

/*
 * Select a pair of neighboring districts i, j
 */
void select_pair(int n_distr, const Graph &dist_g, int &i, int &j);

#endif
