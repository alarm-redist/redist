#pragma once
#ifndef SPLITTING_H
#define SPLITTING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "gredist_types.h"
#include "base_plan_type.h"
#include "tree_splitter_types.h"





/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(const Graph &g, int &k, int const last_k, 
                      const std::vector<double> &unnormalized_weights, double thresh,
                      double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
                      const uvec &counties,
                      Multigraph &cg, const uvec &pop, 
                      int const min_region_cut_size, int const max_region_cut_size,
                      bool split_district_only,
                      double const target, int const verbosity);

//' Splits a multidistrict in all of the plans
//'
//' Using the procedure outlined in <PAPER HERE> this function attempts to split
//' a multidistrict in a previous steps plan until M successful splits have been made. This
//' is based on the `split_maps` function in smc.cpp
//'
//' @title Split all the maps
//'
//' @param g A graph (adjacency list) passed by reference
//' @param counties Vector of county labels of each vertex in `g`
//' @param cg County level multigraph
//' @param pop A vector of the population associated with each vertex in `g`
//' @param old_plans_vec A vector of plans from the previous step
//' @param new_plans_vec A vector which will be filled with plans that had a
//' multidistrict split to make them
//' @param parent_index_vec A vector used to track the index of the previous plan
//' sampled that was successfully split. The value of `parent_index_vec[i]` is the
//' index of the old plan from which the new plan `new_plans_vec[i]` was
//' successfully split from. In other words `new_plans_vec[i]` is equal to
//' `attempt_region_split(old_plans_vec[parent_index_vec[i]], ...)`
//' @param unnormalized_sampling_weights A vector of weights used to sample indices
//' of the `old_plans_vec`. The value of `unnormalized_sampling_weights[i]` is
//' the unnormalized probability that index i is selected
//' @param draw_tries_vec A vector used to keep track of how many plan split
//' attempts were made for index i. The value `draw_tries_vec[i]` represents how
//' many split attempts were made for the i-th new plan (including the successful
//' split). For example, `draw_tries_vec[i] = 1` means that the first split
//' attempt was successful.
//' @param parent_unsuccessful_tries_vec A vector used to keep track of how many times the
//' previous rounds plans were sampled and unsuccessfully split. The value
//' `parent_unsuccessful_tries_vec[i]` represents how many times `old_plans_vec[i]` was sampled
//' and then unsuccessfully split while creating all `M` of the new plans.
//' THIS MAY NOT BE THREAD SAFE
//' @param accept_rate The number of accepted splits over the total number of
//' attempted splits. This is equal to `sum(draw_tries_vec)/M`
//' @param n_unique_parent_indices The number of unique parent indices, ie the
//' number of previous plans that had at least one descendant amongst the new
//' plans. This is equal to `unique(parent_index_vec)`
//' @param ancestors Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lags Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param k_param The top edges to pick parameter for the region splitting
//' algorithm
//' @param split_district_only Whether or not to only allow for single district
//' splits. If set to `true` will only attempt to split off one district at a
//' time
//' @param pool A threadpool for multithreading
//' @param verbosity A parameter controlling the amount of detail printed out
//' during the algorithms running
//'
//' @details Modifications
//'    - The `new_plans_vec` is updated with all the newly split plans
//'    - The `old_plans_vec` is updated with all the newly split plans as well.
//'    Note that the reason both this and `new_plans_vec` are updated is because
//'    of the nature of the code you need both vectors and so both are passed by
//'    reference to save memory.
//'    - The `original_ancestor_vec` is updated to contain the indices of the
//'    original ancestors of the new plans
//'    - The `parent_index_vec` is updated to contain the indices of the parents of the
 //'    new plans
//'    - If two new valid regions are split then the new_region_ids is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'    - The `draw_tries_vec` is updated to contain the number of tries for each
//'    of the new plans
//'    - The `parent_unsuccessful_tries_vec` is updated to contain the number of unsuccessful
//'    samples of the old plans
//'    - The `accept_rate` is updated to contain the average acceptance rate for
//'    this iteration
//'    - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
//'    with the unique number of parents and original ancestors for all the new
//'    plans respectively
//'    - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
//'    I DO NOT KNOW WHAT IT MEANS
//'
//' @noRd
//' @keywords internal
void run_smc_step(
        const MapParams &map_params, 
        std::vector<std::unique_ptr<Plan>> &old_plans_vec, 
        std::vector<std::unique_ptr<Plan>> &new_plans_vec,
        std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
        Rcpp::IntegerMatrix::Column parent_index_vec,
        const arma::vec &normalized_cumulative_weights,
        Rcpp::IntegerMatrix::Column draw_tries_vec,
        Rcpp::IntegerMatrix::Column parent_unsuccessful_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        umat &ancestors, const std::vector<int> &lags,
        int const min_region_cut_size, int const max_region_cut_size, 
        bool const split_district_only,
        RcppThread::ThreadPool &pool,
        int verbosity, int diagnostic_level
);


#endif
