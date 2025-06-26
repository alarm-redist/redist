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
#include "redist_types.h"
#include "splitting.h"
#include "weights.h"
#include "ust_sampler.h"



arma::vec get_adj_pair_unnormalized_weights(
    Plan const &plan,
    std::vector<std::pair<RegionID, RegionID>> const &valid_region_adj_pairs,
    std::string const &selection_type
);


std::tuple<bool, bool, double, int> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &new_plan,
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph,
    std::string const merge_prob_type, bool save_edge_selection_prob,
    std::vector<std::pair<RegionID, RegionID>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts,
    double const rho, bool const is_final
);

int run_merge_split_steps(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &dummy_plan,
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph,
    std::string const merge_prob_type,
    double const rho, bool const is_final, 
    int num_steps_to_run,
    std::vector<int> &tree_sizes, std::vector<int> &successful_tree_sizes
);



#endif
