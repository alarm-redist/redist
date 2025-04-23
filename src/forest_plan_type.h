#pragma once
#ifndef FOREST_PLAN_TYPE_H
#define FOREST_PLAN_TYPE_H

#include "base_plan_type.h"
#include "weights.h"


class ForestPlan : public Plan {


public:
    // constructor for a blank plan 
    ForestPlan(int const V, int const ndists, int const nsim_number,
        int const total_pop,
        AllPlansVector &all_plans_vec, 
        AllRegionSizesVector &all_region_sizes_vec,
        std::vector<int> &all_region_pops_vec,
        std::vector<int> &all_region_order_added_vec
   ):
    Plan(V, ndists, nsim_number, total_pop, 
       all_plans_vec, all_region_sizes_vec, all_region_pops_vec, all_region_order_added_vec
    ){forest_graph.resize(V);};

   // constructor for partial plan (more than 1 region)
   ForestPlan(int const V, int const ndists, int const num_regions,
        int const nsim_number, const arma::uvec &pop,
        AllPlansVector &all_plans_vec, 
        AllRegionSizesVector &all_region_sizes_vec,
        std::vector<int> &all_region_pops_vec,
        std::vector<int> &all_region_order_added_vec,
        const Rcpp::List &initial_forest_adj_list
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
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const override;


    std::vector<std::tuple<int, int, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter,
        std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
    ) const override;

};

#endif