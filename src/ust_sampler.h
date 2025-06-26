#pragma once
#ifndef UST_SAMPLER_H
#define UST_SAMPLER_H

#include <RcppArmadillo.h>
#include "redist_types.h"
#include "tree_op.h"
#include "base_plan_type.h"

class Plan;

// For the future, to avoid needing to create visited and ignore
class USTSampler {

private:


public:

    USTSampler(MapParams const &map_params, SplittingSchedule const &splitting_schedule) : 
        ust(init_tree(map_params.V)), pops_below_vertex(map_params.V, 0), 
        visited(map_params.V), ignore(map_params.V),
        map_params(map_params), splitting_schedule(splitting_schedule){};


    Tree ust;
    std::vector<int> pops_below_vertex;
    std::vector<bool> visited, ignore;
    int root;
    MapParams const &map_params;
    SplittingSchedule const &splitting_schedule;
    

    // Attempts to draw a tree on a region 
    bool attempt_to_draw_tree_on_region(RNGState &rng_state,
        Plan const &plan, const int region_to_draw_tree_on);

    // Attempts to draw a tree on a region formed by merging the two regions
    bool attempt_to_draw_tree_on_merged_region(RNGState &rng_state,
        Plan const &plan, 
        const int region1_to_draw_tree_on, const int region2_to_draw_tree_on);

    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_split(
        RNGState &rng_state, TreeSplitter const &tree_splitter,
        Plan const &plan, int const region_to_split,
        bool const save_selection_prob
    );

    std::pair<bool, EdgeCut> try_to_sample_splittable_tree(
        RNGState &rng_state, TreeSplitter const &tree_splitter,
        int const region_populations, int const region_size,
        bool const save_selection_prob
    );


    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_mergesplit(
        RNGState &rng_state, TreeSplitter const &tree_splitter,
        Plan const &plan, int const merge_region1, int const merge_region2,
        bool const save_selection_prob
    );

};

#endif
