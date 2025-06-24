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
        TreeSplitter const &tree_splitter
    ) const override;


    double get_log_eff_boundary_len(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const override;
};

#endif