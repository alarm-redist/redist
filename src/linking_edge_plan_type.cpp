/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/03
 * Purpose: Implementation of Plan type for Linking Edge space
 ********************************************************/

 #include "linking_edge_plan_type.h"


LinkingEdgePlan::LinkingEdgePlan(arma::subview_col<arma::uword> region_ids_col, 
    arma::subview_col<arma::uword> region_sizes_col, 
    int ndists, int num_regions, const arma::uvec &pop, 
    bool split_district_only,
   const std::vector<std::tuple<int, int, double>> &initial_linking_edge_list
): Plan(region_ids_col, region_sizes_col, ndists, num_regions, pop, split_district_only){

    if(num_regions == 1 || num_regions == ndists){
        linking_edges.reserve(ndists-1);
        forest_graph.resize(region_ids.n_elem);
    }else if(num_regions > 1){
        throw Rcpp::exception("Custom linking edges not ready yet!");
        linking_edges = initial_linking_edge_list;
    }    
}


void LinkingEdgePlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter,
    USTSampler &ust_sampler, EdgeCut const cut_edge, 
    const int split_region1_id, const int split_region2_id
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

    // TODO need to find the edge and update that stuff 
    // Now update the linking edge stuff 

    // First iterate through the old linking edges and update ones touching a split region
    for(auto &a_linking_edge: linking_edges){
        int edge_region1 = region_ids(std::get<0>(a_linking_edge));
        int edge_region2 = region_ids(std::get<0>(a_linking_edge));
        // If either of the vertices in the edge is now in a split region we need 
        // to update the probability 
        if(edge_region1 == split_region1_id || edge_region2 == split_region1_id){
            // std::get<2>(a_linking_edge) = tree_splitter.get_log_retroactive_splitting_prob_for_joined_tree(
                
            // )
        }
    }


    

    // Append the new linking edge
    linking_edges.push_back({cut_edge.cut_vertex, cut_edge.cut_vertex_parent, cut_edge.log_prob});
    // 

    return;
}



double LinkingEdgePlan::get_log_eff_boundary_len(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    TreeSplitter const &tree_splitter, 
    const int region1_id, int const region2_id
) const {
    // Go through and find that pair 
    for (auto const &edge_pair: linking_edges){
        if(std::get<0>(edge_pair) == region1_id && std::get<1>(edge_pair) == region2_id){
            return std::get<2>(edge_pair);
        }
    }
    throw Rcpp::exception("Linking Pair not found!\n");
}