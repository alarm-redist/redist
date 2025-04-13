/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/03
 * Purpose: Implementation of Plan type for Linking Edge space
 ********************************************************/

 #include "linking_edge_plan_type.h"

constexpr bool DEBUG_L_EDGE_PLANS_VERBOSE = false; // Compile-time constant

LinkingEdgePlan::LinkingEdgePlan(arma::subview_col<arma::uword> region_ids_col, 
    arma::subview_col<arma::uword> region_sizes_col, 
    int ndists, int num_regions, const arma::uvec &pop, 
    bool split_district_only,
   const std::vector<std::tuple<int, int, double>> &initial_linking_edge_list
): Plan(region_ids_col, region_sizes_col, ndists, num_regions, pop, split_district_only){

    if(num_regions == 1 || num_regions == ndists){
        linking_edges.reserve(ndists-1);
        forest_graph.resize(region_ids.size());
        // complete hueristic
        for (size_t i = 0; i < forest_graph.size(); i++)
        {
            forest_graph[i].reserve(2);
        }
    }else if(num_regions > 1){
        throw Rcpp::exception("Custom linking edges not ready yet!");
        linking_edges = initial_linking_edge_list;
    }    
}


// NOTE: This really only belongs to Forest and linking edge class, not graph 
double LinkingEdgePlan::get_regions_log_splitting_prob(
    TreeSplitter const &tree_splitter, USTSampler &ust_sampler,
    const int region1_root, const int region2_root
) const{
    int region1_id = region_ids[region1_root]; int region2_id = region_ids[region2_root];
    int region1_size = region_sizes[region1_id]; int region2_size = region_sizes[region2_id];
    int merged_region_size = region1_size + region2_size;
    auto cut_size_bounds = ust_sampler.splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // get the log probability
    return tree_splitter.get_log_retroactive_splitting_prob_for_joined_tree(
        ust_sampler.map_params, forest_graph, 
        ust_sampler.visited, ust_sampler.pops_below_vertex,
        region1_root, region2_root,
        region_pops[region1_id], region_pops[region2_id],
        region1_size, region2_size,
        min_possible_cut_size, max_possible_cut_size,
        ust_sampler.splitting_schedule.all_regions_smaller_cut_sizes_to_try[merged_region_size]
    );
}


void LinkingEdgePlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id,
    bool const add_region
){

    
    // If we're not adding a region we need to erase the old linking edge associated 
    // with these two regions if it exists 
    if(!add_region){
        if(DEBUG_L_EDGE_PLANS_VERBOSE){
        Rprintf("Adding Edge (%d, %d) - Regions (%d, %d)\n",
            cut_edge.cut_vertex, cut_edge.cut_vertex_parent,
            split_region1_id, split_region2_id);
        Rprintf("Current Linking Edges:\n");
        for (auto const &a_edge: linking_edges)
        {
            Rprintf("\tEdge(%d, %d) - Region (%d, %d)- Prob %f \n", 
                std::get<0>(a_edge), std::get<1>(a_edge), 
            region_ids[std::get<0>(a_edge)], region_ids[std::get<1>(a_edge)],
            std::exp(std::get<2>(a_edge)));
        }
        Rprintf("\n");
        }
        int erase_index;
        bool need_to_erase = false;
        std::set<int> merge_region_ids = {split_region1_id, split_region2_id};
        // Go through and find that pair 
        for (int i = 0; i < linking_edges.size(); i++)
        {
            // get the regions associated with the linking edge
            std::set<int> candidate_pairs = {
                static_cast<int>(region_ids[std::get<0>(linking_edges[i])]), 
                static_cast<int>(region_ids[std::get<1>(linking_edges[i])])
            };
            // if they match break
            if(candidate_pairs == merge_region_ids){
                erase_index = i;
                need_to_erase = true;
                break;
            }
        }
        // erase that pair of linking edges if needed
        if(need_to_erase){
            linking_edges.erase(linking_edges.begin()+erase_index);
        }else{
            REprintf("Error could not find old edge!!\n");
            throw Rcpp::exception("NOOOO");
        }
    }


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

    // TODO need to find the edge and update that stuff 
    // Now update the linking edge stuff 

    // First iterate through the old linking edges and update ones touching a split region
    for(auto &a_linking_edge: linking_edges){
        int v = std::get<0>(a_linking_edge); int u = std::get<1>(a_linking_edge);
        int edge_region1 = region_ids[v];
        int edge_region2 = region_ids[u];
        // If either of the vertices in the edge is now in a split region we need 
        // to update the probability that edge was chosen 
        if(edge_region1 == split_region1_id || edge_region2 == split_region1_id ||
            edge_region1 == split_region2_id || edge_region2 == split_region2_id){
            if(DEBUG_L_EDGE_PLANS_VERBOSE){
                Rprintf("(%d,%d) in Region (%d,%d)", v,u, edge_region1, edge_region2);
                Rprintf(" Linking edge prob was %f\n", std::get<2>(a_linking_edge));
            } 
            // update the log probability 
            std::get<2>(a_linking_edge) = get_regions_log_splitting_prob(
                tree_splitter, ust_sampler,
                u, v
            );
            if(DEBUG_L_EDGE_PLANS_VERBOSE) Rprintf("Linking edge prob is now %f\n", std::get<2>(a_linking_edge));
        }
    }

    // Append the new linking edge from the region we just split 
    linking_edges.push_back({cut_edge.cut_vertex, cut_edge.cut_vertex_parent, cut_edge.log_prob});
    if(DEBUG_L_EDGE_PLANS_VERBOSE){
        Rprintf("\nAdding New Edge (%d,%d) in Region (%d,%d), log prob %f\n",
            cut_edge.cut_vertex, cut_edge.cut_vertex_parent,
            region_ids[cut_edge.cut_vertex], region_ids[cut_edge.cut_vertex_parent],
            cut_edge.log_prob);
    } 

    return;
}



double LinkingEdgePlan::get_log_eff_boundary_len(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const {
    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        if(edge_region1 == region1_id && edge_region2 == region2_id){
            return std::get<2>(edge_pair);
        }
    }
    throw Rcpp::exception("Linking Pair not found!\n");
}


std::vector<std::tuple<int, int, double>> LinkingEdgePlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter,
    std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map
) const{
    std::vector<std::tuple<int, int, double>> valid_adj_region_pairs_to_boundary_map;
    valid_adj_region_pairs_to_boundary_map.reserve(linking_edges.size());


    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];
        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            valid_adj_region_pairs_to_boundary_map.push_back(
                {edge_region1, edge_region2, std::get<2>(edge_pair)}
            );
        }
    }
    return(valid_adj_region_pairs_to_boundary_map);
}



int LinkingEdgePlan::count_valid_adj_regions(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule
) const{
    int num_valid_pairs = 0;
    // Go through and count linking edge pairs that can be merged
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];
        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            ++num_valid_pairs;
        }
    }
    return(num_valid_pairs);
}


std::vector<std::pair<int,int>> LinkingEdgePlan::get_valid_adj_regions(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    bool const check_split_constraint
) const{
    std::vector<std::pair<int,int>> valid_adj_region;
    valid_adj_region.reserve(linking_edges.size());


    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        // get the regions associated with the linking edge
        int edge_region1 = std::min(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);
        int edge_region2 = std::max(region_ids[std::get<0>(edge_pair)], region_ids[std::get<1>(edge_pair)]);

        int edge_region1_size = region_sizes[edge_region1]; int edge_region2_size = region_sizes[edge_region2];
        
        // Add if its ok to merge 
        if(splitting_schedule.valid_merge_pair_sizes[edge_region1_size][edge_region2_size]){
            valid_adj_region.push_back(
                {edge_region1, edge_region2}
            );
        }
    }
    return(valid_adj_region);
}