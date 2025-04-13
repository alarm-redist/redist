/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Implementation of classes for tree splitting 
functions
********************************************************/

#include "tree_splitter_types.h"




std::vector<EdgeCut> TreeSplitter::get_all_valid_pop_edge_cuts_in_directed_tree(
    const MapParams &map_params, 
    Tree const &ust, const int root, 
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    int const region_population, int const region_size,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
) const {

    // reset pops_below_vertex and valid edges thing
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    std::fill(no_valid_edges_vertices.begin(), no_valid_edges_vertices.end(), false);
    std::vector<EdgeCut> valid_edges = get_all_valid_edges_in_directed_tree(
        ust, root, map_params.pop, 
        pops_below_vertex, no_valid_edges_vertices,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        region_population, region_size,
        map_params.lower, map_params.upper, map_params.target);


    return valid_edges;
}


std::pair<bool, EdgeCut> TreeSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, RNGState &rng_state,
    Tree const &ust, const int root, 
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    int const region_population, int const region_size,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    bool save_selection_prob
) const {
    // get all the valid edges 
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, 
        ust, root, 
        pops_below_vertex, no_valid_edges_vertices,
        region_population, region_size,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try
    );

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if(num_valid_edges == 0){
        return std::make_pair(false, EdgeCut());
    }else{ // else have derived class choose according to its rule
        return select_edge_to_cut(rng_state, valid_edges, save_selection_prob);
    }
}

// returns edge cut and log probability it was chosen
std::pair<bool, EdgeCut> TreeSplitter::select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob
    ) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately 
    if(num_valid_edges == 1){
        // if only 1 just return that
        // selection prob is just 1 so don't touch
        // if(save_selection_prob){
        //     Rprintf("Save true: %d valid, only 1 edge, log prob is %f \n", 
        //         num_valid_edges, valid_edges[0].log_prob);
        // }
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights 
    arma::vec unnormalized_wgts(num_valid_edges);

    for (size_t i = 0; i < num_valid_edges; i++)
    {
        unnormalized_wgts(i) = compute_unnormalized_edge_cut_weight(
            valid_edges[i]
        );
    }
    
    
    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    double log_selection_prob = 0.0;
    if(save_selection_prob){
        selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
        // Rprintf("Save, %d valid, log prob is %f and %f\n", num_valid_edges, selected_edge_cut.log_prob, 
        //     std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts)));
    }

    return std::make_pair(true, selected_edge_cut);
}


// Takes a vector of valid edge cuts and returns the log probability 
    // the one an index idx would have been chosen 
double TreeSplitter::get_log_selection_prob(
        const std::vector<EdgeCut> &valid_edges,
        int idx
) const{
    auto num_valid_edges = valid_edges.size();
    // get the weights 
    double weight_sum = 0.0;
    // get idx weight
    double idx_weight = compute_unnormalized_edge_cut_weight(valid_edges[idx]);

    // get sum of weights 
    for (size_t i = 0; i < num_valid_edges; i++)
    {
        weight_sum += compute_unnormalized_edge_cut_weight(
            valid_edges[i]
        );
    }

    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return std::log(idx_weight) - std::log(weight_sum);
}


double TreeSplitter::get_log_retroactive_splitting_prob_for_joined_tree(
    MapParams const &map_params,
    VertexGraph const &forest_graph,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_root, const int region2_root,
    const int region1_population, const int region2_population,
    const int region1_size, const int region2_size,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try
) const{
    int total_merged_region_size = region1_size+region2_size;

    // Get all the valid edges in the joined tree 
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_tree(
        map_params, forest_graph, 
        pops_below_vertex, visited,
        region1_root, region1_population,
        region2_root, region2_population,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_size
    );


    // find the index of the actual edge we cut 
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(
        region1_root, region2_root, region1_root, 
        region2_size, region2_population,
        region1_size, region1_population
    );

    // if(TREE_SPLITTING_DEBUG_VERBOSE)
    // Rprintf("%d valid edges!\n", valid_edges.size());

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    if(TREE_SPLITTING_DEBUG_VERBOSE){
    REprintf("Actual Cut Edge at Index %d and so prob is %f \n", 
        actual_cut_edge_index,
        get_log_selection_prob(valid_edges, actual_cut_edge_index));
    }

    return get_log_selection_prob(valid_edges, actual_cut_edge_index);
}


void NaiveTopKSplitter::update_single_int_param(int int_param){
    if(int_param <= 0) throw Rcpp::exception("Splitting k must be at least 1!\n");
    k_param = int_param;
}

std::pair<bool, EdgeCut> NaiveTopKSplitter::select_edge_to_cut(
    RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
    bool save_selection_prob
) const{

    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // if(num_valid_edges > k_param){
    //     REprintf("k was %d but found %d valid edges\n", k_param, num_valid_edges);
    //     // throw Rcpp::exception("K not big enough!\n");
    // }

    int idx = rng_state.r_int(k_param);
    // if we selected k greater than number of edges failure
    if(idx >= num_valid_edges){
        return std::make_pair(false, EdgeCut()); 
    }else{
        // we always store selection probability since its so cheap to compute
        EdgeCut selected_edge_cut = valid_edges[idx];
        selected_edge_cut.log_prob = - std::log(k_param);
        return std::make_pair(true, selected_edge_cut);
    }

}



std::pair<bool, EdgeCut> UniformValidSplitter::select_edge_to_cut(
    RNGState &rng_state,std::vector<EdgeCut> const &valid_edges,
    bool save_selection_prob
) const{
    int num_valid_edges  = static_cast<int>(valid_edges.size());
    // pick one unif at random 
    int idx = rng_state.r_int(num_valid_edges);
    // we always store selection probability since its so cheap to compute
    EdgeCut selected_edge_cut = valid_edges[idx];
    selected_edge_cut.log_prob = - std::log(num_valid_edges);
    return std::make_pair(true, selected_edge_cut);
}



double ExpoWeightedSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut
) const{
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double bigger_dev = std::max(devs.at(0), devs.at(1));
    return std::exp(-alpha*bigger_dev);
}


double ExpoWeightedSmallerDevSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut
) const{
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double smaller_dev = std::min(devs.at(0), devs.at(1));
    return std::exp(-alpha*smaller_dev);
}



std::pair<bool, EdgeCut> ExperimentalSplitter::select_edge_to_cut(
        RNGState &rng_state,std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob
    ) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately 
    if(num_valid_edges == 1){
        // if only 1 just return that
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights 
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, target);
    
    
    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    double log_selection_prob = 0.0;
    if(save_selection_prob){
        selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
    }

    return std::make_pair(true, selected_edge_cut);
}


double ExperimentalSplitter::get_log_selection_prob(
    const std::vector<EdgeCut> &valid_edges,
    int idx
    ) const{
    // get the weights 
    arma::vec unnormalized_wgts = compute_almost_best_weights_on_smaller_dev_edges(
        valid_edges, epsilon, target);
    
    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));
}