#pragma once
#ifndef SMC_H
#define SMC_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "labeling.h"

/*
 * Main entry point.
 *
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
List smc_plans(int N, List l, const arma::uvec &counties, const arma::uvec &pop,
               int n_distr, double target, double lower, double upper, double rho,
               arma::umat districts, int n_drawn, int n_steps,
               List constraints, List control, int verbosity=1);

/*
 * Split off a piece from each map in `districts`,
 * keeping deviation between `lower` and `upper`
 */
void split_maps(const Graph &g, const uvec &counties, Multigraph &cg,
                const uvec &pop, umat &districts, vec &cum_wgt, vec &lp,
                vec &pop_left, vec &log_temper, double pop_temper,
                double &accept_rate, int n_distr, int dist_ctr,
                std::vector<Graph> &dist_grs, vec &log_labels,
                umat &ancestors, const std::vector<int> &lags,
                bool adjust_labels, double est_label_mult, int &n_unique,
                double lower, double upper, double target,
                double rho, int k, bool check_both,
                RcppThread::ThreadPool &pool, int verbosity);


/*
 * Add specific constraint weights & return the cumulative weight vector
 */
vec get_wgts(const umat &districts, int n_distr, int distr_ctr, bool final,
             double alpha, vec &lp, double &neff,
             const uvec &pop, double parity, const Graph g,
             List constraints, int verbosity);

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
void adapt_parameters(const Graph &g, int &k, int last_k, const vec &lp, double thresh,
                      double tol, const umat &districts, const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      const vec &pop_left, double target, int verbosity);

#endif
