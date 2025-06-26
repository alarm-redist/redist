#include "forest_plan_type.h"


// constructor for partial plan (more than 1 region)
ForestPlan::ForestPlan(
    int const ndists, int const num_regions,
    const arma::uvec &pop,
    PlanVector &this_plan_region_ids, 
    RegionSizes &this_plan_region_sizes,
    IntPlanAttribute &this_plan_region_pops,
    IntPlanAttribute &this_plan_order_added,
    MapParams const &map_params, 
    Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore,
    RNGState &rng_state, const Rcpp::List &initial_forest_adj_list
):
Plan(num_regions, pop, 
    this_plan_region_ids, this_plan_region_sizes, this_plan_region_pops, this_plan_order_added
){
    // resize the forest graph 
    forest_graph.resize(region_ids.size());

    if(initial_forest_adj_list.size() > 1){
        throw Rcpp::exception("Input Forest list not supported right now\n");
    }else{
        // else just build a forest at random 
        for (size_t region_id = 0; region_id < num_regions; region_id++)
        {
            int root;
            auto result = draw_tree_on_region(
                map_params, region_id, ust, visited, ignore,
                root, rng_state, 1000000
            );

            if(!result.first){
                throw Rcpp::exception("Could not draw a tree on a region after 1000000 attempts");
            } 
        }
    
    } 
}


VertexGraph ForestPlan::get_forest_adj(){
    return forest_graph;
}



void ForestPlan::update_vertex_and_plan_specific_info_from_cut(
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
    // update the vertex labels and the tree
    assign_region_id_and_forest_from_tree(
        ust_sampler.ust, region_ids, forest_graph,
        split_region1_tree_root, split_region1_id);

    assign_region_id_and_forest_from_tree(
        ust_sampler.ust, region_ids, forest_graph,
        split_region2_tree_root, split_region2_id);

    return;
}



/*  
 * @title Get Valid Adjacent Regions to Effective Tree Boundary Lengths Map
 * 
 * Returns a unordered_map mapping valid pairs of adjacent regions to 
 * the length of the effect tree boundary between them in the forest
 * 
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid 
 * to merge) and returns a hash map mapping the pairs of valid adjacent 
 * regions to the length of the effective tree boundary between them in `g`.
 * Recall this is the probability that we chose the actual split edge among
 * the tree made by merging each boundary edge
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
 * 
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to double) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the effective tree boundary length 
 * between them in `g`.
 */
std::vector<std::tuple<RegionID, RegionID, double>> compute_log_tree_eff_boundary_lens(
    const ForestPlan &plan, VertexGraph const &forest_graph,
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    const TreeSplitter &edge_splitter
){
    std::vector<bool> visited(plan_multigraph.map_params.V);
    std::vector<int> pops_below_vertex(plan_multigraph.map_params.V, 0);
    
    int const V = plan_multigraph.map_params.V;

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        auto v_region_num = plan.region_ids[v];
        auto v_region_size = plan.region_sizes[v_region_num];

        // check if its a region we want to find regions adjacent to 
        // and if not keep going
        if(!splitting_schedule.check_adj_to_regions.at(v_region_size)) continue;

        
        // now iterate over its neighbors
        for (auto const v_nbor : plan_multigraph.map_params.g[v]) {
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
            auto search_result = plan_multigraph.pair_map.get_value(v_nbor_region_num, v_region_num);

            // plan_multigraph.Rprint();
            // REprintf("(%d, %d) - (%u, %u)\n", 
            //     v, v_nbor, v_region_num, v_nbor_region_num);

            // REprintf("Its %s and %s\n",
            //     (search_result.first ? "TRUE" : "FALSE"),
            // (std::get<3>(search_result.second) ? "TRUE" : "FALSE") );
            
            // ignore if not a pair in the map or we're ignoring it for other reasons
            if(!search_result.first || !search_result.second.count_pair) continue;


            // ignore if hierarchically adjacent but this edge isn't in the same county
            if(
               (plan_multigraph.map_params.counties[v] != plan_multigraph.map_params.counties[v_nbor]) &&
               search_result.second.admin_adjacent
            ) continue;

            // Now we can get the effective boundary length             
            auto merged_region_size = v_region_size + v_nbor_region_size;
            auto cut_size_bounds = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
            auto min_possible_cut_size = cut_size_bounds.first;
            auto max_possible_cut_size = cut_size_bounds.second;

            double log_edge_selection_prob = edge_splitter.get_log_retroactive_splitting_prob_for_joined_tree(
                    plan_multigraph.map_params, forest_graph,
                    visited, pops_below_vertex,
                    v, v_nbor,
                    plan.region_pops[plan.region_ids[v]], plan.region_pops[plan.region_ids[v_nbor]],
                    plan.region_sizes[plan.region_ids[v]], plan.region_sizes[plan.region_ids[v_nbor]],
                    min_possible_cut_size, max_possible_cut_size,
                    splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]);

            // we add to the boundary 
            plan_multigraph.pair_map.add_to_eff_boundary(
                v_region_num, v_nbor_region_num,
                std::exp(log_edge_selection_prob)
            );

            // Rprintf("Adding (%d,%d) w/ %.4f\n", v, v_nbor, std::exp(log_edge_selection_prob));
        }
    }

    // now make the output vector 
    std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(plan_multigraph.pair_map.num_hashed_pairs);

    for(auto const key_val_pair: plan_multigraph.pair_map.get_all_values()){
        region_pairs_tuple_vec.emplace_back(
            key_val_pair.first.first, 
            key_val_pair.first.second, 
            std::log(key_val_pair.second.eff_boundary_len)
        );
    }


    return region_pairs_tuple_vec;
}



std::vector<std::tuple<RegionID, RegionID, double>> ForestPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter
) const{
    // build the multigraph 
    plan_multigraph.build_plan_multigraph(*this);
    // remove invalid hierarchical merges if needed
    plan_multigraph.remove_invalid_hierarchical_merge_pairs(*this);
    // remove invalid size pairs 
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    // get pairs and log tree effective boundary
    auto region_pairs_tuple_vec = compute_log_tree_eff_boundary_lens(
        *this, forest_graph,
        plan_multigraph, splitting_schedule,
         tree_splitter
    );

    return region_pairs_tuple_vec;

}




double ForestPlan::get_log_eff_boundary_len(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const{

    std::vector<bool> visited(plan_multigraph.map_params.V);
    std::vector<int> pops_below_vertex(plan_multigraph.map_params.V);

    int const V = plan_multigraph.map_params.V;
    int const merged_region_size = region_sizes[region1_id]+ region_sizes[region2_id];

    auto cut_size_bounds = splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // figure out if we are allowed to count across county boundaries 
    bool count_edges_across;
    if(plan_multigraph.counties_on){
        auto search_result = plan_multigraph.pair_map.get_value(region1_id, region2_id);
        // only count across if not admin adjacent 
        count_edges_across = !search_result.second.admin_adjacent;
    }else{
        count_edges_across = true; 
    }
    
    double tree_selection_probs = 0.0;
    for (int v = 0; v < V; v++) {
        if (region_ids[v] != region1_id) continue; // Only count if starting vertex in region 1
        auto v_county = plan_multigraph.map_params.counties[v];
        for (auto nbor : plan_multigraph.map_params.g[v]) {
            // ignore if not region 2
            if (region_ids[nbor] != region2_id) continue;
            // ignore if we can't count this boundary
            if(v_county != plan_multigraph.map_params.counties[nbor] && !count_edges_across) continue;

            double log_edge_selection_prob = tree_splitter.get_log_retroactive_splitting_prob_for_joined_tree(
                plan_multigraph.map_params, forest_graph,
                visited, pops_below_vertex,
                v, nbor,
                region_pops[region_ids[v]], region_pops[region_ids[nbor]],
                region_sizes[region_ids[v]], region_sizes[region_ids[nbor]],
                min_possible_cut_size,max_possible_cut_size,
                splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]);
            
            tree_selection_probs += std::exp(log_edge_selection_prob);
        }
    }

    return std::log(tree_selection_probs);
}