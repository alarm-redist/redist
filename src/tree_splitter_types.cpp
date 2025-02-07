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
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int region_id_to_split
){
    // reset pops_below_vertex
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    // get population below
    get_tree_pops_below(ust, root, map_params.pop, pops_below_vertex);

    std::vector<EdgeCut> valid_edges = get_all_valid_edges_in_directed_tree(ust, root, 
        pops_below_vertex, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        plan.region_pops.at(region_id_to_split), plan.region_sizes(region_id_to_split),
        map_params.lower, map_params.upper, map_params.target);


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
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int region_id_to_split){


    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        region_id_to_split
    );

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // iif no valid edges immediately return false
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }// else if(num_valid_edges > k_param){
    //     REprintf("k was %d but found %d valid edges\n", k_param, num_valid_edges);
    //     // throw Rcpp::exception("K not big enough!\n");
    // }

    int idx = r_int(k_param);

    // if we selected k greater than number of edges failure
    if(idx >= num_valid_edges){
        return std::make_pair(false, EdgeCut()); 
    }else{
        return std::make_pair(true, valid_edges[idx]);
    }

}



std::pair<bool,EdgeCut> UniformValidSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try, 
        const int region_id_to_split){

    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
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

}


std::pair<bool,EdgeCut> ExpoWeightedSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int region_id_to_split){
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
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


double ExpoWeightedSplitter::get_log_selection_prob(
    const MapParams &map_params,
    const std::vector<EdgeCut> &valid_edges,
    int idx
    ) const{
    // get the weights 
    arma::vec unnormalized_wgts = compute_expo_prob_weights_on_edges(
        valid_edges, alpha, map_params.target);
    
    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));

}


std::pair<bool,EdgeCut> ExpoWeightedSmallerDevSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int region_id_to_split){
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
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
    arma::vec unnormalized_wgts = compute_expo_prob_weights_on_smaller_dev_edges(
        valid_edges, alpha, map_params.target);
    
    // select with prob proportional to the weights
    int idx = r_int_unnormalized_wgt(unnormalized_wgts);
    // idx = static_cast<int>(unnormalized_wgts.index_max());
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // EdgeCut selected_edge_cut = valid_edges.at(unnormalized_wgts.index_max());
    return std::make_pair(true, selected_edge_cut);

}


double ExpoWeightedSmallerDevSplitter::get_log_selection_prob(
    const MapParams &map_params,
    const std::vector<EdgeCut> &valid_edges,
    int idx
    ) const{
    // get the weights 
    arma::vec unnormalized_wgts = compute_expo_prob_weights_on_smaller_dev_edges(
        valid_edges, alpha, map_params.target);
    
    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));

}


std::pair<bool,EdgeCut> ExperimentalSplitter::select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int region_id_to_split){
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, plan,
        ust, root, 
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
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
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, map_params.target);
    
    // select with prob proportional to the weights
    int idx = r_int_unnormalized_wgt(unnormalized_wgts);
    // idx = static_cast<int>(unnormalized_wgts.index_max());
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // EdgeCut selected_edge_cut = valid_edges.at(unnormalized_wgts.index_max());
    return std::make_pair(true, selected_edge_cut);

}


double ExperimentalSplitter::get_log_selection_prob(
    const MapParams &map_params,
    const std::vector<EdgeCut> &valid_edges,
    int idx
    ) const{
    // get the weights 
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, map_params.target);
    
    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));

}