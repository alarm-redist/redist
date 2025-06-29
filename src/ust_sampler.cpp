/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/3
* Purpose: Encapsulation of uniform spanning tree sampler functions
********************************************************/

#include "ust_sampler.h"


bool USTSampler::attempt_to_draw_tree_on_region(
    RNGState &rng_state,
    Plan const &plan, const int region_to_draw_tree_on){
    int V = map_params.V;

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_ids[i] != region_to_draw_tree_on;
    }

    // get upper and lower bounds on region pops
    auto min_max_pair = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[
        plan.region_sizes[region_to_draw_tree_on]
    ];

    // clear the tree
    clear_tree(ust);
    // Get a uniform spanning tree drawn on that region
    int result = sample_sub_ust(map_params.g, ust, 
        V, root, visited, ignore, 
        map_params.pop, 
        min_max_pair.first* map_params.lower, min_max_pair.second*map_params.upper, 
        map_params.counties, map_params.cg,
        rng_state);
    // result == 0 means it was successful
    return(result == 0);
}


bool USTSampler::attempt_to_draw_tree_on_merged_region(RNGState &rng_state,
    Plan const &plan, 
    const int region1_to_draw_tree_on, const int region2_to_draw_tree_on){
    int V = map_params.V;

    // Mark it as ignore if its not in either of the two regions
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_ids[i] != region1_to_draw_tree_on && plan.region_ids[i] != region2_to_draw_tree_on;
    }

    int merged_region_size = plan.region_sizes[region1_to_draw_tree_on] + plan.region_sizes[region2_to_draw_tree_on];
    // get upper and lower bounds on region pops
    auto min_max_pair = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[
        merged_region_size
    ];


    // clear the tree
    clear_tree(ust);
    // Get a uniform spanning tree drawn on that region
    int result = sample_sub_ust(map_params.g, ust, 
        V, root, visited, ignore, 
        map_params.pop, 
        min_max_pair.first* map_params.lower, min_max_pair.second*map_params.upper,
        map_params.counties, map_params.cg,
        rng_state);
    // result == 0 means it was successful
    return(result == 0);
}



std::pair<bool, EdgeCut> USTSampler::try_to_sample_splittable_tree(
    RNGState &rng_state, TreeSplitter const &tree_splitter,
    int const region_populations, int const region_size,
    bool const save_selection_prob
){
    // We assume a tree has already been successfully drawn so
    // try to find a valid cut
    auto cut_size_bounds = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // REprintf("Remainder Size: %d Cut sizes:", region_size);
    // for(auto const &v: splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size]){
    //     REprintf("%d, ", v);
    // }
    // REprintf("\n");

    std::pair<bool, EdgeCut> edge_search_result = tree_splitter.attempt_to_find_edge_to_cut(
        map_params, rng_state,
        ust, root, pops_below_vertex, ignore,
        region_populations, region_size,
        min_possible_cut_size, max_possible_cut_size,
        splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size],
        save_selection_prob
    );

    bool search_successful = std::get<0>(edge_search_result);
    // return false if unsuccessful 
    if(!search_successful) return std::make_pair(false, EdgeCut());

    // If successful extract the edge cut info
    EdgeCut cut_edge = std::get<1>(edge_search_result);
    // Now erase the cut edge in the tree
    erase_tree_edge(ust, cut_edge);

    return edge_search_result;
}


std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_split(
    RNGState &rng_state, TreeSplitter const &tree_splitter,
    Plan const &plan, int const region_to_split,
    bool const save_selection_prob
){
    // Try to draw a tree
    bool tree_drawn = attempt_to_draw_tree_on_region(
        rng_state,
        plan, region_to_split
    );
    // return false if unsuccessful 
    if(!tree_drawn) return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size = plan.region_sizes[region_to_split];
    int region_to_split_population = plan.region_pops[region_to_split];

    return try_to_sample_splittable_tree(
        rng_state, tree_splitter,
        region_to_split_population, region_to_split_size,
        save_selection_prob
    );
}



std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_mergesplit(
    RNGState &rng_state, TreeSplitter const &tree_splitter,
    Plan const &plan, int const merge_region1, int const merge_region2,
    bool const save_selection_prob
){
    // Try to draw a tree
    bool tree_drawn = attempt_to_draw_tree_on_merged_region(
        rng_state,
        plan, merge_region1, merge_region2
    );
    // return false if unsuccessful 
    if(!tree_drawn) return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size = plan.region_sizes[merge_region1]+plan.region_sizes[merge_region2];
    int region_to_split_population = plan.region_pops[merge_region1]+plan.region_pops[merge_region2];

    return try_to_sample_splittable_tree(
        rng_state, tree_splitter,
        region_to_split_population, region_to_split_size,
        save_selection_prob
    );
}