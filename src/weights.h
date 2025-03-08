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
#include "scoring.h"
#include "map_calc.h"


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
    std::vector<std::tuple<int, int, double>> const &adj_pairs_and_boundary_lens,
    std::string const &selection_type
);


double get_log_mh_ratio(
    const Graph &g, 
    const int region1_id, const int region2_id,
    const int old_region_boundary_length, const int new_region_boundary_length, 
    const double current_pair_merge_prob, const double new_pair_merge_prob, 
    Plan &current_plan, Plan &new_plan
);





void get_all_plans_uniform_adj_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    const std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    std::vector<double> &unnormalized_sampling_weights,
    int verbosity
);



double get_log_retroactive_splitting_prob_for_joined_tree(
    MapParams const &map_params,
    Plan const &plan, const TreeSplitter &edge_splitter,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_root, const int region2_root,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
);





void compute_all_plans_log_optimal_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    const std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    std::vector<double> &unnormalized_sampling_weights,
    int verbosity
);


#endif
