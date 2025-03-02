#include <algorithm>
#include <set>
#include <RcppThread.h>
#include "smc_base.h"
#include "tree_op.h"
#include "base_plan_type.h"

#ifndef MAP_CALC_H
#define MAP_CALC_H

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_graph_boundary(const Graph &g, const subview_col<uword> &districts,
                    int const region1_id, int const region2_id);


/*  
 * @title Count the Number of Valid Adjacent Region Pairs
 * 
 * Counts and returns the number of valid adjacent region pairs in a graph. 
 * We define a valid of pair of adjacent regions to be two regions that 
 *      1. Share at least one edge between them in `g`
 *      2. The regions sizes are valid to merge with eligibility determined
 *         by the `valid_merge_pairs` matrix
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return The number of valid adjacent region pairs
 */
int count_valid_adj_regions(
    Graph const &g, Plan const &plan,
    std::vector<bool> const &check_adj_to_regions,
    std::vector<std::vector<bool>> const &valid_merge_pairs
);


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
NumericMatrix group_pct(arma::umat m, arma::vec group_pop, arma::vec total_pop, int n_distr);

/*
 * Tally a variable by district.
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
                             double const total_pop, double const target);


#endif
