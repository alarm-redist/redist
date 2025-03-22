#pragma once
#ifndef MERGING_H
#define MERGING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <set>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "gredist_types.h"
#include "splitting.h"
#include "weights.h"
#include "ust_sampler.h"




std::tuple<bool, bool, double> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &new_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    std::string const merge_prob_type, bool save_edge_selection_prob,
    std::vector<std::pair<int,int>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts,
    double const rho, bool const is_final
);

int run_merge_split_steps(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &dummy_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    std::string const merge_prob_type,
    double const rho, bool const is_final, 
    int num_steps_to_run
);

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec, 
    TreeSplitter const &tree_splitter,
    std::string const merge_prob_type, 
    double const rho, bool const is_final, 
    int const nsteps_to_run,
    Rcpp::IntegerMatrix::Column success_count_vec,
    int verbosity
);

#endif
