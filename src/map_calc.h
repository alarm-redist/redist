#include <algorithm>
#include <set>
#include <RcppThread.h>
#include "smc_base.h"
#include "tree_op.h"


#ifndef MAP_CALC_H
#define MAP_CALC_H



/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const subview_col<uword> &districts, int distr,
                     const uvec &total_pop, mat ssdmat, double denominator);





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
 * Infer the number of seats of the regions 
 */
// [[Rcpp::export]]
Rcpp::IntegerMatrix infer_region_seats(
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



// computes log number of spanning trees on region intersect county
// In either a region or a merged region 
double compute_log_region_and_county_spanning_tree(
    Graph const &g, const uvec &counties, int const county,
    PlanVector const &region_ids,
    int const region1_id, int const region2_id
);


/*
 * Compute the log number of spanning trees for the contracted (ie county level) graph
 */
double compute_log_county_level_spanning_tree(
    Graph const &g, const uvec &counties, int const n_cty,
    PlanVector const &region_ids,
    int const region1_id, int const region2_id
);

// [[Rcpp::export]]
Rcpp::NumericVector order_district_stats(
    Rcpp::NumericVector const &district_stats, 
    int const ndists,
    int const num_threads
);



// [[Rcpp::export]]
Rcpp::DataFrame order_columns_by_district(
    Rcpp::DataFrame const &df,
    Rcpp::CharacterVector const &columns,
    int const ndists,
    int const num_threads = 0);

#endif
