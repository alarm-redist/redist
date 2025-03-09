#pragma once
#ifndef UST_SAMPLER_H
#define UST_SAMPLER_H

#include <RcppArmadillo.h>
#include "gredist_types.h"
#include "tree_op.h"
#include "base_plan_type.h"

// For the future, to avoid needing to create visited and ignore
class USTSampler {

public:

    USTSampler(int V) : 
        ust(init_tree(V)), pops_below_vertex(V, 0), visited(V), ignore(V){};


    Tree ust;
    std::vector<int> pops_below_vertex;
    std::vector<bool> visited, ignore;
    int root;


    // Attempts to draw a tree on a region 
    bool draw_tree_on_region(const MapParams &map_params, RNGState &rng_state,
        Plan const &plan, const int region_to_draw_tree_on);

    std::tuple<bool, EdgeCut, double> attempt_to_find_valid_tree_split(
        const MapParams &map_params, SplittingSchedule const &splitting_schedule,
        RNGState &rng_state, TreeSplitter const &tree_splitter,
        Plan const &plan, int const region_to_split,
        bool const save_selection_prob
    );

};

#endif