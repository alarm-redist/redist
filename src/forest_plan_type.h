#pragma once
#ifndef FOREST_PLAN_TYPE_H
#define FOREST_PLAN_TYPE_H

#include "base_plan_type.h"
#include "weights.h"


class ForestPlan : public Plan {


public:
    // constructor for a blank plan 
    ForestPlan(int const total_seats,
        int const total_pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added
   ):
    Plan(total_seats, total_pop, 
        this_plan_region_ids, this_plan_region_sizes, this_plan_region_pops, this_plan_order_added
    ){forest_graph.resize(region_ids.size());};

   // constructor for partial plan (more than 1 region)
   ForestPlan(int const ndists, int const num_regions,
        const arma::uvec &pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added,
        MapParams const &map_params, 
        Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore,
        RNGState &rng_state, const Rcpp::List &initial_forest_adj_list = {}
    );
    

    // We now need to keep track of trees as undirected graphs
    VertexGraph get_forest_adj() override;

    void update_vertex_and_plan_specific_info_from_cut (
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler, EdgeCut const cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool const add_region
    ) override;

    double get_log_eff_boundary_len(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const override;


    std::vector<std::tuple<RegionID, RegionID, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter
    ) const override;

};

#endif