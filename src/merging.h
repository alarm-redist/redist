#pragma once
#ifndef MERGING_H
#define MERGING_H

#include "advanced_types.h"

#include <RcppArmadillo.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class Plan;
class PlanMultigraph;
class ScoringFunction;
class SplittingSchedule;
class RNGState;
class USTSampler;
class TreeSplitter;
class WeightCache;

// Simple struct for tracking granular time
struct GranularMCMCTimes {
    double wilson_time = 0.0;
    double backfill_wilson_time = 0.0;
    double get_valid_mergepairs = 0.0;
    double selecting_merge_pair = 0.0;
    double hard_constraint_time = 0.0;
    double eff_boundary_length = 0.0;
    double region_scores = 0.0;
    double plan_scores = 0.0;
    double tau_terms = 0.0;
    double plan_copying = 0.0;
};



arma::vec get_adj_pair_unnormalized_weights(
    Plan const &plan, std::vector<std::pair<RegionID, RegionID>> const &valid_region_adj_pairs,
    std::string const &selection_type);

std::tuple<bool, bool, double, int> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, RNGState &rng_state,
    SamplingSpace const sampling_space, Plan &plan, Plan &new_plan, USTSampler &ust_sampler,
    TreeSplitter &tree_splitter, PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph, std::string const merge_prob_type,
    bool save_edge_selection_prob, std::vector<std::pair<RegionID, RegionID>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts, double const rho, bool const is_final, bool const do_mh,
    bool const using_caching, WeightCache *weight_cache,
    GranularMCMCTimes &granular_times);

int run_merge_split_steps(MapParams const &map_params,
                          const SplittingSchedule &splitting_schedule,
                          ScoringFunction const &scoring_function, RNGState &rng_state,
                          SamplingSpace const sampling_space, Plan &plan, Plan &dummy_plan,
                          USTSampler &ust_sampler, TreeSplitter &tree_splitter,
                          PlanMultigraph &current_plan_multigraph,
                          PlanMultigraph &proposed_plan_multigraph,
                          std::string const merge_prob_type, double const rho,
                          bool const is_final, int num_steps_to_run,
                          std::vector<int> &tree_sizes, std::vector<int> &successful_tree_sizes,
                          bool const using_caching, WeightCache *weight_cache, 
                          GranularMCMCTimes &granular_times);

#endif
