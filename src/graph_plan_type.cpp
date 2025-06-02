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
std::vector<std::tuple<RegionID, RegionID, double>> NEW_get_valid_pairs_and_log_graph_boundary_len_map(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    const Plan &plan,
    CountyComponents &county_components,
    EffBoundaryMap &pair_map
){    
    int const V = map_params.V;

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        auto v_region_num = plan.region_ids[v];
        auto v_region_size = plan.region_sizes[v_region_num];

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!splitting_schedule.check_adj_to_regions.at(v_region_size)) continue;

        // now iterate over its neighbors
        for (auto const v_nbor : map_params.g.at(v)) {
            // find which region neighbor corresponds to
            auto v_nbor_region_num = plan.region_ids[v_nbor];
            // ignore if same region 
            if(v_region_num == v_nbor_region_num) continue;

            auto v_nbor_region_size = plan.region_sizes[v_nbor_region_num];
            // else they are different so regions are adj
            // check if we need to be aware of double counting. IE if both regions
            // are ones where check adjacent is true then its possible to double count edges
            bool const double_counting_possible = splitting_schedule.check_adj_to_regions[v_nbor_region_size];
            // ignore if double counting is possible and v is smaller region 
            if(double_counting_possible && v_region_num < v_nbor_region_num) continue;

            // see if this pair is one we count 
            auto search_result = pair_map.get_value(v_region_num, v_nbor_region_num);
            // ignore if not hashing or if we don't count this boundary 
            if(!search_result.first || 
                !county_components.count_county_boundary(
                    v_region_num, map_params.counties(v),
                    v_nbor_region_num, map_params.counties(v_nbor)
                )
            ) continue;
            // Now we add 1 for this boundary 

            // set the value to be the current boundary length plus 
            pair_map.set_value(
                v_region_num, v_nbor_region_num, 
                search_result.second + 1.0
            );
        }
    }

    // now make the output vector 
    std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(pair_map.num_hashed_values);

    for(auto const &a_pair: pair_map.hashed_pairs){
        auto search_result = pair_map.get_value(a_pair.first, a_pair.second);
        if(search_result.second == 0.0){
            REprintf("Pair (%u,%u) = 0\n", a_pair.first, a_pair.second);
        }
        // add the region ids and the log of the boundary 
        region_pairs_tuple_vec.emplace_back(
            a_pair.first, 
            a_pair.second,
            std::log(search_result.second)
        );
    }

    return region_pairs_tuple_vec;
}



std::vector<std::tuple<RegionID, RegionID, double>> GraphPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, CountyComponents &county_components,
    EffBoundaryMap &pair_map
) const{
    
    return NEW_get_valid_pairs_and_log_graph_boundary_len_map(
        map_params, splitting_schedule, *this,
        county_components, pair_map
    );

}




/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 * 
 * If county constraints are on then it won't count any boundaries in invalid counties
 */
double count_log_graph_boundary(
    PlanVector const &region_ids,
    int const region1_id, int const region2_id,
    CountyComponents &county_components
){
    int const V = county_components.map_params.V;
    double count = 0.0;

    for (int v = 0; v < V; v++) {
        auto v_region = region_ids[v];
        if (v_region != region1_id) continue; // Only count if starting vertex in region 1
        auto v_county = county_components.map_params.counties[v];
        for (int u : county_components.map_params.g[v]) {
            int u_region = region_ids[u];
            // ignore if not the right id
            if (region_ids[u] != region2_id) continue;
            // ignore if we can't count this boundary
            if(!county_components.count_county_boundary(
                region1_id, v_county, 
                region2_id, county_components.map_params.counties[u]
                )) continue;
            count += 1.0;
        }
    }

    return std::log(count);
}

double GraphPlan::get_log_eff_boundary_len(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, CountyComponents &county_components,
    const int region1_id, int const region2_id
) const{
    // Return the log of the graph theoretic boundary 
    return count_log_graph_boundary(region_ids,
        region1_id, region2_id, county_components
    );
}