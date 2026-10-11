#pragma once
#ifndef FOREST_PLAN_TYPE_H
#define FOREST_PLAN_TYPE_H

#include "base_plan_type.h"

class USTSampler;
class PlanMultigraph;
class ScoringFunction;
class TreeSplitter;
class SplittingSchedule;
class RNGState;

class ForestPlan : public Plan {

  public:
    // constructor for a blank plan
    ForestPlan(int const total_seats, int const total_pop, PlanVector &this_plan_region_ids,
               RegionSizes &this_plan_region_sizes, IntPlanAttribute &this_plan_region_pops,
               IntPlanAttribute &this_plan_order_added, PlanEdgeBits &this_plan_forest_edge_bits)
        : Plan(total_seats, total_pop, this_plan_region_ids, this_plan_region_sizes,
               this_plan_region_pops, this_plan_order_added, this_plan_forest_edge_bits) {
    };

    // constructor for partial plan (more than 1 region). If
    // `initial_forest_words` is passed the forest is copied from those packed
    // words (unchecked), otherwise a tree is drawn on each region at random.
    ForestPlan(int const ndists, int const num_regions, const std::vector<unsigned int> &pop,
                       PlanVector &this_plan_region_ids, RegionSizes &this_plan_region_sizes,
                       IntPlanAttribute &this_plan_region_pops,
                       IntPlanAttribute &this_plan_order_added,
                       PlanEdgeBits &this_plan_forest_edge_bits,
                       USTSampler &ust_sampler,
                       RNGState &rng_state,
               EdgeBitWord const *initial_forest_words = nullptr);

    // We now need to keep track of trees as undirected graphs

    void update_vertex_and_plan_specific_info_from_cut(
        TreeSplitter const &tree_splitter, USTSampler &ust_sampler, EdgeCut const cut_edge,
        const int split_region1_id, const int split_region2_id, bool const add_region) override;

    double get_log_eff_boundary_len(PlanMultigraph &plan_multigraph,
                                    const SplittingSchedule &splitting_schedule, ScheduleContext &schedule_ctx,
                                    USTSampler &ust_sampler, TreeSplitter &tree_splitter,
                                    ScoringFunction const &scoring_function,
                                    const int region1_id, int const region2_id) const override;

    std::vector<std::tuple<RegionID, RegionID, double>>
    get_valid_adj_regions_and_eff_log_boundary_lens(PlanMultigraph &plan_multigraph,
                                                    const SplittingSchedule &splitting_schedule, ScheduleContext &schedule_ctx,
                                                    ScoringFunction const &scoring_function,
                                                    bool const is_final_split,
                                                    USTSampler &ust_sampler,
                                                    TreeSplitter &tree_splitter) const override;
};

#endif