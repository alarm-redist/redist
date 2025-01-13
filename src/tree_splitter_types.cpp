/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Implementation of classes for tree splitting 
functions
********************************************************/

#include "tree_splitter_types.h"

std::pair<bool,EdgeCut> NaiveTopKSplitter::select_edge_to_cut(
        const MapParams &map_params, 
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const arma::subview_col<arma::uword> &region_ids, 
        const int region_id_to_split, const int total_region_pop, const int total_region_size){

    // Get population below and parent vertices
    int V = static_cast<int>(region_ids.n_elem);

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V);
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);


    return get_naive_top_k_edge(root, pop_below, tree_vertex_parents,
                     k_param, min_potential_cut_size, max_potential_cut_size,
                     region_ids, region_id_to_split, 
                     total_region_pop, total_region_size,
                     map_params.lower, map_params.upper, map_params.target);

}