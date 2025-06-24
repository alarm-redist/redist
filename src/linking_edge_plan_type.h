#pragma once
#ifndef LINKING_PLAN_TYPE_H
#define LINKING_PLAN_TYPE_H

#include "base_plan_type.h"
#include "weights.h"

class LinkingEdgePlan : public Plan {

    public:
    // constructor for a blank plan 
    LinkingEdgePlan(int const total_seats, 
        int const total_pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added
   );

   // constructor for partial plan (more than 1 region)
   LinkingEdgePlan(    
    int const ndists, int const num_regions,
    const arma::uvec &pop,
    PlanVector &this_plan_region_ids, 
    RegionSizes &this_plan_region_sizes,
    IntPlanAttribute &this_plan_region_pops,
    IntPlanAttribute &this_plan_order_added,
    const std::vector<std::array<double, 3>> &linking_edges,
    const Rcpp::List &initial_forest_adj_list
    );


        VertexGraph get_forest_adj() override;

        void Rprint(bool verbose = false) const override;

        void update_vertex_and_plan_specific_info_from_cut(
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

        // For a plan and edge across two regions it returns the log probability 
        // that edge would be split in the merged tree across the two regions
        double get_regions_log_splitting_prob(
            TreeSplitter const &tree_splitter, USTSampler &ust_sampler,
            const int region1_root, const int region2_root
        ) const;

        // This builds a plan multigraph and 
        // This just returns pairs with a linking edge between them and valid size 
        std::pair<bool, std::vector<std::pair<RegionID,RegionID>>> attempt_to_get_valid_mergesplit_pairs(
            PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
        ) const override;

        // Get a vector of all valid adj region pairs for the backwards kernel
        std::vector<std::pair<RegionID,RegionID>> get_valid_smc_merge_regions(
            PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
        ) const override;

        std::vector<std::array<double, 3>> get_linking_edges() override{
            std::vector<std::array<double, 3>> output;
            output.reserve(linking_edges.size());
            for(auto const& an_edge: linking_edges){
                output.push_back({
                    static_cast<double>(std::get<0>(an_edge)),
                    static_cast<double>(std::get<1>(an_edge)),
                    std::get<2>(an_edge)
                });
            }
            return output;
        };

};

#endif