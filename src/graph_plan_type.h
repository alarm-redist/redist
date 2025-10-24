#pragma once
#ifndef GRAPH_PLAN_TYPE_H
#define GRAPH_PLAN_TYPE_H

#include "base_plan_type.h"
#include "tree_op.h"
#include "weights.h"

class GraphPlan : public Plan {
    // private member variable
    

public:
    using Plan::Plan;
    // implementation of the pure virtual function
    

    void update_vertex_and_plan_specific_info_from_cut(
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler, EdgeCut const cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool const add_region
    ) override;


    std::vector<std::tuple<RegionID, RegionID, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        ScoringFunction const &scoring_function, bool const is_final_split,
        USTSampler &ust_sampler, TreeSplitter &tree_splitter
    ) const override;


    double get_log_eff_boundary_len(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        USTSampler &ust_sampler, TreeSplitter &tree_splitter, 
        ScoringFunction const &scoring_function,
        const int region1_id, int const region2_id
    ) const override;
};


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    RNGState &rng_state,
    int &k, int const last_k, 
    const arma::vec &unnormalized_weights, double thresh,
    double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
    bool split_district_only,
    int const verbosity);



/*
 * Choose k and multiplier for efficient, accurate sampling
 */
int estimate_mergesplit_cut_k(
    Plan const &plan, PlanMultigraph const &plan_multigraph,  
    SplittingSchedule const &splitting_schedule,
    double const thresh, double const tol,
    RNGState &rng_state
);


#endif