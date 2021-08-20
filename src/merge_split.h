#ifndef MERGESPLIT_H
#define MERGESPLIT_H

#include "smc_base.h"
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "kirchhoff.h"

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
                    double upper, double rho, double beta_sq, const arma::uvec &current,
                    int n_current, double beta_vra, double tgt_min, double tgt_other,
                    double pow_vra, const arma::uvec &min_pop,
                    double beta_vra_hinge, const arma::vec &tgts_min,
                    double beta_inc, const arma::uvec &incumbents,
                    double beta_splits, double beta_fractures,
                    double thresh, int k, int verbosity);

/*
 * Add specific constraint weights & return the cumulative weight vector
 */
double calc_gibbs_tgt(const subview_col<uword> &plan, int n_distr, int V,
                      int distr_1, int distr_2, const uvec &pop, double beta_sq,
                      const uvec &current, int n_current,
                      double beta_vra, double tgt_min, double tgt_other,
                      double pow_vra, const uvec &min_pop,
                      double beta_vra_hinge, const vec &tgts_min,
                      double beta_inc, const uvec &incumbents,
                      double beta_counties, double beta_fractures,
                      const uvec &counties, int n_cty);


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
void select_pair(int n, const Graph &g, const uvec &plan, int &i, int &j);

#endif
