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
#include "redist_types.h"
#include "base_plan_type.h"
#include "splitting_schedule_types.h"
#include "tree_op.h"
#include "scoring.h"
#include "map_calc.h"
#include "graph_ops.h"


/* Computes Compute Effective Sample Size from log incremental weights
 *
 *
 * Takes a vector of log incremental weights and computes the effective sample
 * size which is the sum of the weights squared divided by the sum of squared
 * weights
 * 
 * @param log_wgt vector of log incremental weights
 * 
 * @details No modifications to inputs made
 * 
 * @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
 */
double compute_n_eff(const arma::subview_col<double> log_wgt);



void compute_all_plans_log_simple_incremental_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    TreeSplitter const &tree_splitter, 
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    int verbosity
);

double compute_log_optimal_incremental_weights(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, const TreeSplitter &edge_splitter,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double const rho,
    bool compute_log_splitting_prob, bool is_final_plan
);

void compute_all_plans_log_optimal_incremental_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    TreeSplitter const &tree_splitter,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    int verbosity
);


#endif
