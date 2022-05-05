#include <algorithm>
#include <set>
#include <RcppThread.h>
#include "smc_base.h"
#include "tree_op.h"

#ifndef MAP_CALC_H
#define MAP_CALC_H

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_boundary(const Graph &g, const subview_col<uword> &districts,
                    int distr_root, int distr_other);

/*
 * Compute the status quo penalty for district `distr`
 */
double eval_sq_entropy(const subview_col<uword> &districts, const uvec &current,
                       int distr, const uvec &pop, int n_distr, int n_current, int V);

/*
 * Compute the new, hinge VRA penalty for district `distr`
 */
double eval_grp_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop);

/*
 * Compute the new, hinge VRA penalty for district `distr`
 */
double eval_grp_inv_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop);

/*
 * Compute the old VRA penalty for district `distr`
 */
double eval_grp_pow(const subview_col<uword> &districts, int distr,
                    const uvec &grp_pop, const uvec &total_pop,
                    double tgt_grp, double tgt_other, double pow);

/*
 * Compute the incumbent-preserving penalty for district `distr`
 */
double eval_inc(const subview_col<uword> &districts, int distr, const uvec &incumbents);

/*
 * Compute the county split penalty for district `distr`
 */
double eval_splits(const subview_col<uword> &districts, int distr,
                   const uvec &counties, int n_cty, bool smc);

/*
 * Compute the county fracture penalty for district `distr`
 */
double eval_multisplits(const subview_col<uword> &districts, int distr,
                        const uvec &counties, int n_cty, bool smc);

/*
 * Compute the county split penalty for district `distr`
 */
double eval_total_splits(const subview_col<uword> &districts, int distr,
                   const uvec &counties, int n_cty);

/*
 * Compute the Polsby Popper penalty for district `distr`
 */
double eval_polsby(const subview_col<uword> &districts, int distr,
            const ivec &from,
            const ivec &to,
            const vec &area,
            const vec &perimeter);

/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const subview_col<uword> &districts, int distr,
                     const uvec &total_pop, mat ssdmat, double denominator);

/*
 * Compute the population penalty for district `distr`
 */
double eval_pop_dev(const subview_col<uword> &districts, int distr,
                       const uvec &total_pop, double parity);

/*
 * Compute the segregation penalty for district `distr`
 */
double eval_segregation(const subview_col<uword> &districts, int distr,
                        const uvec &grp_pop, const uvec &total_pop);

/*
 * Compute the qps penalty for district `distr`
 */
double eval_qps(const subview_col<uword> &districts, int distr,
                const uvec &total_pop, const uvec &cities, int n_city,
                int nd);

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_log_st(const subview_col<uword> &districts, const Graph g,
                   arma::uvec counties, int ndists);

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_er(const subview_col<uword> &districts, const Graph g, int ndists);



/*
 * Compute the cooccurence matrix for a set of precincts indexed by `idxs`,
 * given a collection of plans
 */
// [[Rcpp::export]]
arma::mat prec_cooccur(arma::umat m, arma::uvec idxs, int ncores=0);

/*
 * Compute the percentage of `group` in each district. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
NumericMatrix group_pct(arma::umat m, arma::vec group_pop, arma::vec total_pop, int n_distr);

/*
 * Compute the deviation from the equal population constraint.
 */
// [[Rcpp::export]]
NumericMatrix pop_tally(IntegerMatrix districts, arma::vec pop, int n_distr);

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// [[Rcpp::export]]
NumericVector max_dev(const IntegerMatrix districts, const arma::vec pop, int n_distr);

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
std::vector<double> tree_dev(Tree &ust, int root, const uvec &pop,
                             double total_pop, double target);


#endif
