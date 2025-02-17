#pragma once
#ifndef WEIGHTS_H
#define WEIGHTS_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"


#include <string>
#include <random>
#include <unordered_map>
#include <set>
#include <utility> // for std::pair
#include <cmath>
#include "gredist_types.h"
#include "base_plan_type.h"
#include "splitting_schedule_types.h"
#include "tree_op.h"


//' Computes the effective sample size from log incremental weights
//'
//' Takes a vector of log incremental weights and computes the effective sample
//' size which is the sum of the weights squared divided by the sum of squared
//' weights
//'
//'
//' @title Compute Effective Sample Size
//'
//' @param log_wgt vector of log incremental weights
//'
//' @details No modifications to inputs made
//'
//' @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
//'
double compute_n_eff(const arma::subview_col<double> log_wgt);



//' Returns a vector of the triple (smaller region id, bigger region id, boundary len)
//' for all valid pairs of adjacent regions in the plan. (Either all adjacent regions if
//' doing generalized region splits or just adjacent to the remainder if only doing 
//' one district splits.)
//'
//'
//' @title Get All Valid Adjacent Regions and their Boundary Length
//'
//' @param g A graph (adjacency list) passed by reference
//' @param plan A plan object
//' @param split_district_only If true only gets regions adjacent to the remainder but if 
//' false then gets all adjacent regions in the plan
//'
//' @details No modifications to inputs made
//'
//' @return A vector of integer arrays of size 3 where the values are
//' (smaller region id, bigger region id, boundary len)
//'
std::vector<std::array<int, 3>> get_valid_adj_regions_and_boundary_lens_vec(
    Graph const &g, Plan const &plan,
    bool const split_district_only 
);

void get_all_adj_pairs(
    Graph const &g, std::vector<std::pair<int, int>> &adj_pairs_vec,
    arma::subview_col<arma::uword> const &vertex_region_ids,
    std::vector<bool> const valid_regions
);



// Returns a sampler over a vector of adjacent pairs where the probability 
// of a pair is decided according to `selection_type`
//'
//' Current supported options are
//'     - uniform - Every pair has equal probability
//'     - district_pair - double district pairs have weight 1000, one district is 10,
//'         and two multidistricts have 1/(1+sum of their dvals)
//'
//' @title Get Sampler over Adj Regions List
//'
//' @param plan A plan object
//' @param adj_pairs_and_boundary_lens A vector where each pair is 
//' (adj region1, adj region2, boundary length between 2 regions)
//' @param selection_type A string controlling the function to use
//' in assigning the unnormalized weight to each pair
//'
//' @details No modifications to inputs made
//'
//' @return A sampler where index i has probability proportional to the weight 
//' given to that pair 
//'
std::discrete_distribution<>  get_adj_pair_sampler(
    Plan const &plan,
    std::vector<std::array<int, 3>> const &adj_pairs_and_boundary_lens,
    std::string const &selection_type
);


double get_log_mh_ratio(
    const Graph &g, 
    const int region1_id, const int region2_id,
    const int old_region_boundary_length, const int new_region_boundary_length, 
    const double current_pair_merge_prob, const double new_pair_merge_prob, 
    Plan &current_plan, Plan &new_plan
);

//' Computes log unnormalized weights for vector of plans
//'
//' Using the procedure outlined in <PAPER HERE> this function computes the log
//' incremental weights and the unnormalized weights for a vector of plans (which
//' may or may not be the same depending on the parameters).
//'
//' @title Compute Log Unnormalized Weights
//'
//' @param pool A threadpool for multithreading
//' @param g A graph (adjacency list) passed by reference
//' @param plans_ptr_vec A vector of plans to compute the log unnormalized weights
//' of
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not. If `split_district_only` is true
//' then uses optimal weights from one-district split scheme.
//' @param log_incremental_weights A vector of the log incremental weights
//' computed for the plans. The value of `log_incremental_weights[i]` is
//' the log incremental weight for `plans_vec[i]`
//' @param unnormalized_sampling_weights A vector of the unnormalized sampling
//' weights to be used with sampling the `plans_vec` in the next iteration of the
//' algorithm. Depending on the other hyperparameters this may or may not be the
//' same as `exp(log_incremental_weights)`
//' @param target Target population of a single district
//' @param pop_temper <DETAILS NEEDED>
//'
//' @details Modifications
//'    - The `log_incremental_weights` is updated to contain the incremental
//'    weights of the plans
//'    - The `unnormalized_sampling_weights` is updated to contain the unnormalized
//'    sampling weights of the plans for the next round
void get_all_plans_log_optimal_weights(
        RcppThread::ThreadPool &pool,
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
        bool split_district_only,
        arma::subview_col<double> log_incremental_weights,
        std::vector<double> &unnormalized_sampling_weights,
        double pop_temper
);



void get_all_plans_uniform_adj_weights(
        RcppThread::ThreadPool &pool,
        const SplittingSchedule &splitting_schedule,
        const Graph &g, std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
        bool split_district_only,
        arma::subview_col<double> log_incremental_weights,
        std::vector<double> &unnormalized_sampling_weights,
        double target, double pop_temper
);



double get_log_retroactive_splitting_prob_for_joined_tree(
    MapParams const &map_params,
    Plan const &plan, Tree &ust, const TreeSplitter &edge_splitter,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_root, const int region2_root,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
);



double compute_optimal_forest_log_incremental_weight(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        const Plan &plan, const TreeSplitter &edge_splitter,
        bool split_district_only,
        const double pop_temper);



void get_all_forest_plans_log_optimal_weights(
        RcppThread::ThreadPool &pool,
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
        const std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
        bool split_district_only,
        arma::subview_col<double> log_incremental_weights,
        std::vector<double> &unnormalized_sampling_weights,
        double pop_temper, int verbosity
);



double compute_optimal_log_incremental_weight(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    const Plan &plan, bool split_district_only,
    const double pop_temper);

#endif
