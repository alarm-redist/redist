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
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
arma::umat smc_plans(int N, List l, const arma::uvec &counties,
                     const arma::uvec &pop, int n_distr, double target,
                     double lower, double upper, double rho,
                     double beta_sq, const arma::uvec &current, int n_current,
                     double beta_vra, double tgt_min, double tgt_other,
                     double pow_vra, const arma::uvec &min_pop,
                     double beta_vra_hinge, const arma::vec &tgts_min,
                     double beta_inc, const arma::uvec &incumbents,
                     arma::vec &lp, double thresh,
                     double alpha, double pop_temper=0.1, int verbosity=1);

/*
 * Split off a piece from each map in `districts`,
 * keeping deviation between `lower` and `upper`
 */
void split_maps(const Graph &g, const uvec &counties, Multigraph &cg,
                const uvec &pop, umat &districts, vec &cum_wgt, vec &lp,
                vec &pop_left, vec &log_temper, double pop_temper, int n_distr,
                int dist_ctr, double lower, double upper, double target,
                double rho, int k, int verbosity);


/*
 * Add specific constraint weights & return the cumulative weight vector
 */
vec get_wgts(const umat &districts, int n_distr, int distr_ctr,
             double alpha, vec &lp, const uvec &pop,
             double beta_sq, const uvec &current, int n_current,
             double beta_vra, double tgt_min, double tgt_other,
             double pow_vra, const uvec &min_pop,
             double beta_vra_hinge, const vec &tgts_min,
             double beta_inc, const uvec &incumbents, int verbosity);

/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map(const Graph &g, const uvec &counties, Multigraph &cg,
                 subview_col<uword> districts, int dist_ctr, const uvec &pop,
                 double total_pop, double &lower, double upper, double target, int k);

/*
 * Cut spanning subtree into two pieces of roughly equal population
 */
double cut_districts(Tree &ust, int k, int root, subview_col<uword> &districts,
                     int dist_ctr, const uvec &pop, double total_pop,
                     double lower, double upper, double target);

/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_parameters(const Graph &g, int &k, const vec &lp, double thresh,
                      double tol, const umat &districts, const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      const vec &pop_left, double target, int verbosity);

#endif
