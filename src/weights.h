#pragma once
#ifndef WEIGHTS_H
#define WEIGHTS_H




#include "advanced_types.h" // SamplingSpace


class Plan;
class PlanMultigraph;
class SplittingSchedule;
class ScoringFunction;
class USTSampler;
class TreeSplitter;
class WeightCache;

// Simple struct for tracking granular time
struct GranularWeightTimes {
    double get_valid_pairs = 0.0;
    double splitting_prob = 0.0;
    double region_scores = 0.0;
    double plan_scores = 0.0;
    double tau_terms = 0.0;
};

double compute_log_optimal_incremental_weights(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, USTSampler &ust_sampler,
    TreeSplitter &edge_splitter, SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double const rho,
    double const whole_map_compactness_term, bool compute_log_splitting_prob,
    double const multidistrict_selection_alpha,
    bool is_final_split, bool const using_caching, WeightCache *weight_cache,
    GranularWeightTimes &granular_times);
    

double compute_simple_log_incremental_weight(Plan const &plan, PlanMultigraph &plan_multigraph,
                                             const SplittingSchedule &splitting_schedule,
                                             USTSampler &ust_sampler,
                                             TreeSplitter &edge_splitter,
                                             SamplingSpace const sampling_space,
                                             ScoringFunction const &scoring_function,
                                             double rho, bool compute_log_splitting_prob,
                                             double const multidistrict_selection_alpha,
                                             bool is_final_split, bool const using_caching,
                                             WeightCache *weight_cache,
                                             GranularWeightTimes &granular_times);

#endif
