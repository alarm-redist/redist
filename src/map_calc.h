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
 * 
 * If county constraints are on then it won't count any boundaries in invalid counties
 */
double log_graph_boundary(const Graph &g, const subview_col<uword> &region_ids,
                    int const region1_id, int const region2_id, 
                    int const num_counties, arma::uvec counties);



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
                   const uvec &counties, int n_cty, bool smc);

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
NumericMatrix group_pct(IntegerMatrix const &plans_mat, 
    arma::vec const &group_pop, arma::vec const &total_pop, 
    int const n_distr, int const num_threads = 0);

/*
 * Tally a variable by district.
 */
// [[Rcpp::export]]
NumericMatrix pop_tally(IntegerMatrix const &districts, arma::vec const &pop, int const n_distr,
                        int const num_threads = 0);


/*
 * Infer the sizes of the regions 
 */
// [[Rcpp::export]]
Rcpp::IntegerMatrix infer_region_sizes(
    Rcpp::IntegerMatrix const &region_pops,
    double const lower, double const upper,
    int const total_seats,
    int const num_threads = 0
);

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// [[Rcpp::export]]
NumericVector max_dev(const IntegerMatrix &districts, const arma::vec &pop, int const n_distr,
                      int const num_threads = 0);

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
std::vector<double> tree_dev(Tree &ust, int root, const uvec &pop,
                             double const total_pop, double const target);


// computes log number of spanning trees on region intersect county
// In either a region or a merged region 
double compute_log_region_and_county_spanning_tree(
    Graph const &g, const uvec &counties, int const county,
    PlanVector const &region_ids,
    int const region1_id, int const potential_region2_id = -42
);


/*
 * Compute the log number of spanning trees for the contracted (ie county level) graph
 */
double compute_log_county_level_spanning_tree(
    Graph const &g, const uvec &counties, int const n_cty,
    PlanVector const &region_ids,
    int const region1_id, int const potential_region2_id = -42
);

// [[Rcpp::export]]
Rcpp::NumericVector order_district_stats(
    Rcpp::NumericVector const &district_stats, 
    int const ndists,
    int const num_threads
);


/**************************
 * Parallel Versions of redistmetric functions for working with very large plans 
 * Can probably remove in production version
 ****************************/

/*
 * Compute the number of edges removed
 * 
 * Parallel version lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
// [[Rcpp::export]]
NumericVector parallel_n_removed(const Graph &g, const IntegerMatrix &districts, int const n_distr, 
    int const num_threads);


// lifted from 
// https://github.com/alarm-redist/redistmetrics/blob/5f7b36d8a7f9c7bc3c9098a7b6c6aa561d9074c7/src/partisan.cpp#L72
// [[Rcpp::export(rng = false)]]
NumericVector parallel_effgap(
    NumericMatrix const &dcounts, NumericMatrix const &rcounts, 
    int const totvote, int const num_threads
);


// lifted from 
// https://github.com/alarm-redist/redistmetrics/blob/5f7b36d8a7f9c7bc3c9098a7b6c6aa561d9074c7/src/partisan.cpp#L72
// [[Rcpp::export(rng = false)]]
NumericMatrix parallel_agg_p2d(
    IntegerMatrix const &dm, NumericVector const &vote, 
    int const nd, int const num_threads);

// lifted from 
// https://github.com/alarm-redist/redistmetrics/blob/5f7b36d8a7f9c7bc3c9098a7b6c6aa561d9074c7/src/partisan.cpp#L72
// [[Rcpp::export(rng = false)]]
NumericVector parallel_biasatv(
    NumericMatrix const &dvs, 
    double const v, int const nd, int const num_threads);




/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/kirchhoff.cpp
 */
// [[Rcpp::export(rng = false)]]
NumericMatrix parallelDVS(
    NumericMatrix const &dcounts, NumericMatrix const &rcounts,
    int const num_threads);




/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/splits.cpp
 */
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector parallel_splits(
    const IntegerMatrix &dm, const IntegerVector &community,
    int const nd, int const max_split, int const num_threads,
    bool const skip_last);


/*
 * Parallel version, lifted directly from here
 * https://github.com/alarm-redist/redistmetrics/blob/main/src/compactness.cpp
 */
// [[Rcpp::export(rng = false)]]
NumericMatrix parallel_polsbypopper(IntegerVector const &from,
                           IntegerVector const &to,
                           NumericVector const &area,
                           NumericVector const &perimeter,
                           IntegerMatrix const &dm,
                           int const nd,
                           int const num_threads);

#endif
