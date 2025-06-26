#include "graph_plan_type.h"


constexpr bool DEBUG_GRAPH_PLANS_VERBOSE = false; // Compile-time constant

void GraphPlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id,
    bool const add_region
){
    // Get the root of the tree associated with region 1 and 2
    int split_region1_tree_root, split_region2_tree_root;
    int split_region1_size, split_region2_size;
    int split_region1_pop, split_region2_pop;

    cut_edge.get_split_regions_info(
        split_region1_tree_root, split_region1_size, split_region1_pop,
        split_region2_tree_root, split_region2_size, split_region2_pop
    );

    // update the vertex labels
    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region1_tree_root, split_region1_id);

    assign_region_id_from_tree(ust_sampler.ust, region_ids,
        split_region2_tree_root, split_region2_id);

    return;
}





/*  
 * @title Get Valid Adjacent Regions to Boundary Lengths Map
 * 
 * Returns a unordered_map mapping valid pairs of adjacent regions to 
 * the length of the boundary between them in the graph
 * 
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid 
 * to merge) and returns a hash map mapping the pairs of valid adjacent 
 * regions to the length of the boundary between them in `g`
 * 
 * 
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should 
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find 
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 * @param existing_pair_map A hash map mapping pairs of regions to a double. 
 * This is an optional parameter and if its empty then its ignored. If its not
 * empty then pairs already in the hash won't be counted added to the output.
 * 
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to int) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the length of the boundary
 * between them in `g`.
 */




std::vector<std::tuple<RegionID, RegionID, double>> GraphPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter
) const{

    // build the multigraph 
    plan_multigraph.build_plan_multigraph(*this);
    // remove invalid hierarchical merges if needed
    plan_multigraph.remove_invalid_hierarchical_merge_pairs(*this);
    // remove invalid size pairs 
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    // now make the output vector 
    std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(plan_multigraph.pair_map.num_hashed_pairs);

    for(auto const key_val_pair: plan_multigraph.pair_map.get_all_values()){
        
        double log_boundary_len;
        if(key_val_pair.second.admin_adjacent){
            // if administratively adjacent only count within same county
            log_boundary_len = std::log(key_val_pair.second.within_county_edges);
        }else{
            // else only count across county 
            log_boundary_len = std::log(key_val_pair.second.across_county_edges);
        }

        region_pairs_tuple_vec.emplace_back(
            key_val_pair.first.first, 
            key_val_pair.first.second, 
            log_boundary_len
        );
    }

    return region_pairs_tuple_vec;
    
}


double GraphPlan::get_log_eff_boundary_len(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const{
    // Return the log of the graph theoretic boundary
    auto search_result = plan_multigraph.pair_map.get_value(region1_id, region2_id);

    if(search_result.second.admin_adjacent){
        // if administratively adjacent only count within same county
        return std::log(search_result.second.within_county_edges);
    }else{
        // else only count across county 
        return std::log(search_result.second.across_county_edges);
    }
}