/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Implementation of classes for tree splitting 
functions
********************************************************/

#include "tree_splitter_types.h"

std::vector<EdgeCut> TreeSplitter::get_all_valid_pop_edge_cuts_in_directed_tree(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const int region_id_to_split
){
    // reset pops_below_vertex
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    // get population below
    get_tree_pops_below(ust, root, map_params.pop, pops_below_vertex);

    std::vector<EdgeCut> valid_edges = NEW2_get_all_valid_edges_in_directed_tree(ust, root, 
        pops_below_vertex, 
        min_potential_cut_size, max_potential_cut_size,
        plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
        map_params.lower, map_params.upper, map_params.target);


    // OLD JUST UNTIL I WRITE TESTS then everything below can be deleted
    // Get population below and parent vertices

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(map_params.V, 0);
    std::vector<int> tree_vertex_parents(map_params.V);
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);

    std::vector<EdgeCut> valid_edges2 = get_all_valid_edges_in_directed_tree(root, 
                        pop_below, tree_vertex_parents,
                        min_potential_cut_size, max_potential_cut_size,
                        plan.region_ids,
                        region_id_to_split, 
                        plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
                        map_params.lower, map_params.upper, map_params.target);

    std::set<EdgeCut> valid_edges_set1(valid_edges.begin(), valid_edges.end());
    std::set<EdgeCut> valid_edges_set2(valid_edges2.begin(), valid_edges2.end());

    if(valid_edges_set1 != valid_edges_set2){
        throw Rcpp::exception("Not equal!!");
    }

    return valid_edges;
}


void NaiveTopKSplitter::update_single_int_param(int int_param){
    if(int_param <= 0) throw Rcpp::exception("Splitting k must be at least 1!\n");
    k_param = int_param;
}

std::pair<bool,EdgeCut> NaiveTopKSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){

    // create vector that points to parents & get population below each vtx
    // reset pops_below_vertex
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    // don't need to reset parent vector 
    vertex_parents.at(root) = -1;
    tree_pop(ust, root, map_params.pop, pops_below_vertex, vertex_parents);


    return get_naive_top_k_edge(root, pops_below_vertex, vertex_parents,
                     k_param, min_potential_cut_size, max_potential_cut_size,
                     plan.region_ids, region_id_to_split, 
                     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
                     //total_region_pop, total_region_size,
                     map_params.lower, map_params.upper, map_params.target);

    // // create list that points to parents & computes population below each vtx
    // std::vector<int> pop_below(map_params.V, 0);
    // std::vector<int> tree_vertex_parents(map_params.V);
    // tree_vertex_parents.at(root) = -1;
    // tree_pop(ust, root, map_params.pop, pop_below, tree_vertex_parents);

    // std::vector<EdgeCut> valid_edges = NEW2_get_all_valid_edges_in_directed_tree(ust, root, 
    //     pop_below, 
    //     min_potential_cut_size, max_potential_cut_size,
    //     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
    //     map_params.lower, map_params.upper, map_params.target);

    // std::vector<EdgeCut> valid_edges2 = get_all_valid_edges_in_directed_tree(root, 
    //                     pop_below, tree_vertex_parents,
    //                     min_potential_cut_size, max_potential_cut_size,
    //                     plan.region_ids,
    //                     region_id_to_split, 
    //                     plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
    //                     map_params.lower, map_params.upper, map_params.target);

    
    // std::set<EdgeCut> valid_edges_set1(valid_edges.begin(), valid_edges.end());
    // std::set<EdgeCut> valid_edges_set2(valid_edges2.begin(), valid_edges2.end());

    // if(valid_edges_set1 != valid_edges_set2){
    //     throw Rcpp::exception("Not equal!!");
    // }

}



std::pair<bool,EdgeCut> UniformValidSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){

    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        region_id_to_split
    );

    int num_valid_edges = static_cast<int>(valid_edges.size());

    // if no valid edges reject immediately 
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }else if(num_valid_edges == 1){
        // if only 1 just return that
        return std::make_pair(true, valid_edges.at(0));
    }else{
        int random_idx = r_int(num_valid_edges);
        return std::make_pair(true, valid_edges.at(random_idx));
    }


    // get edge unif at random
    // return get_unif_valid_edge(root, pop_below, tree_vertex_parents,
    //                 min_potential_cut_size, max_potential_cut_size,
    //                  plan.region_ids, region_id_to_split, 
    //                  plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
    //                  map_params.lower, map_params.upper, map_params.target);                   

}


std::pair<bool,EdgeCut> ExpoWeightedSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split){
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        region_id_to_split
    );

    int num_valid_edges = static_cast<int>(valid_edges.size());

    // if no valid edges reject immediately 
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }else if(num_valid_edges == 1){
        // if only 1 just return that
        return std::make_pair(true, valid_edges.at(0));
    }

    // get the weights 
    arma::vec unnormalized_wgts = compute_expo_prob_weights_on_edges(
        valid_edges, alpha, map_params.target);
    
    // select with prob proportional to the weights
    int idx = r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // EdgeCut selected_edge_cut = valid_edges.at(unnormalized_wgts.index_max());

    return std::make_pair(true, selected_edge_cut);

               

}