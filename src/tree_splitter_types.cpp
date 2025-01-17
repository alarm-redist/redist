/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Implementation of classes for tree splitting 
functions
********************************************************/

#include "tree_splitter_types.h"


void NaiveTopKSplitter::update_single_int_param(int int_param){
    if(int_param <= 0) throw Rcpp::exception("Splitting k must be at least 1!\n");
    k_param = int_param;
}

std::pair<bool,EdgeCut> NaiveTopKSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){

    // Get population below and parent vertices
    int V = static_cast<int>(plan.region_ids.n_elem);

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V);
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);


    return get_naive_top_k_edge(root, pop_below, tree_vertex_parents,
                     k_param, min_potential_cut_size, max_potential_cut_size,
                     plan.region_ids, region_id_to_split, 
                     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
                     //total_region_pop, total_region_size,
                     map_params.lower, map_params.upper, map_params.target);

}



std::pair<bool,EdgeCut> UniformValidSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){

    // Get population below and parent vertices
    int V = map_params.V;

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V);
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);

    // get edge unif at random
    return get_unif_valid_edge(root, pop_below, tree_vertex_parents,
                    min_potential_cut_size, max_potential_cut_size,
                     plan.region_ids, region_id_to_split, 
                     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
                     map_params.lower, map_params.upper, map_params.target);                   

}


std::pair<bool,EdgeCut> ExpoWeightedSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){

    // Get population below and parent vertices
    int V = map_params.V;

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V);
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);


    // get edge with exponential probability
    return get_valid_edge_w_expo_prob(root, pop_below, tree_vertex_parents,
                    alpha,
                    min_potential_cut_size, max_potential_cut_size,
                     plan.region_ids, region_id_to_split, 
                     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
                     map_params.lower, map_params.upper, map_params.target);

                     

}