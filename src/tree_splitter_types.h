#pragma once
#ifndef TREE_SPLITTER_TYPES_H
#define TREE_SPLITTER_TYPES_H

#include <RcppArmadillo.h>
#include "redist_types.h"
#include "tree_op.h"
#include "tree_splitting.h"


// [[Rcpp::depends(RcppArmadillo)]]



// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut 
class TreeSplitter {

public:
    // Default Constructor 
    TreeSplitter(int V) {};
    virtual ~TreeSplitter() = default; 

    // Returns a vector of all the valid edges in the tree 
    std::vector<EdgeCut> get_all_valid_pop_edge_cuts_in_directed_tree(
        MapParams const &map_params, 
        Tree const &ust, const int root, 
        std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
        int const region_population, int const region_size,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try
    ) const;
    
    // Takes a spanning tree and returns the edge to cut if successful
    std::pair<bool, EdgeCut> attempt_to_find_edge_to_cut(
        const MapParams &map_params, RNGState &rng_state,
        Tree const &ust, const int root, 
        std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
        int const region_population, int const region_size,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try,
        bool save_selection_prob = false
    ) const;

    // Get probability a specific edge was cut in the tree made by joining
    // the trees in the two regions 
    double get_log_retroactive_splitting_prob_for_joined_tree(
        MapParams const &map_params,
        VertexGraph const &forest_graph, 
        std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
        const int region1_root, const int region2_root,
        const int region1_population, const int region2_population,
        const int region1_size, const int region2_size,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try
    ) const;



    // Takes a vector of valid edge cuts and returns the log probability 
    // the one an index idx would have been chosen 
    virtual double get_log_selection_prob(
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const;

    // used to update the k parameter for top k splitter
    virtual void update_single_int_param(int int_param){
        throw Rcpp::exception("Update single int param not implemented!\n");
    };

    // used to get the k parameter for top k splitter
    virtual int get_single_int_param() const {
        throw Rcpp::exception("Update single int param not implemented!\n");
        return -1;
    };

    // returns edge cut and log probability it was chosen
    virtual std::pair<bool, EdgeCut> select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob
    ) const;

    virtual double compute_unnormalized_edge_cut_weight(
        EdgeCut const &edge_cut
    ) const {throw Rcpp::exception("Not implemented for this class!");};

};


// Splitting method that just tries to pick one of the top k edges unif 
class NaiveTopKSplitter : public TreeSplitter{

public:
    // Constructor for NaiveTopKSplitter
    NaiveTopKSplitter(int V, int k_param)
        : TreeSplitter(V), k_param(k_param) {
    }

    // Attributes specific to NaiveTopKSplitter
    int k_param;        // Top k valuex

    // how to update the k param
    void update_single_int_param(int int_param) override;


    int get_single_int_param()const override{
        return k_param;
    };

    std::pair<bool, EdgeCut> select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob
    ) const override;

    double get_log_selection_prob(
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const override {throw Rcpp::exception("No log selection prob implemented for naive k!\n");}


};



// Splitting method that just tries to pick one of the top k edges unif 
class UniformValidSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    UniformValidSplitter(int V): TreeSplitter(V){};


    std::pair<bool, EdgeCut> select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob 
    ) const override;

    // since uniform log prob is just -log(# of candidates)
    double get_log_selection_prob(
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const override{return -std::log(valid_edges.size());}
};


// Splitting method that picks edge w/ prob ∝ exp(-alpha*bigger dev)
class ExpoWeightedSplitter : public TreeSplitter{

public:

    ExpoWeightedSplitter(int V, double alpha, double target)
        : TreeSplitter(V), alpha(alpha), target(target) {
        if(alpha < 0.0) throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;
    double target;

    virtual double compute_unnormalized_edge_cut_weight(
        EdgeCut const &edge_cut
    ) const override;
};



class ExpoWeightedSmallerDevSplitter : public TreeSplitter{

public:

    ExpoWeightedSmallerDevSplitter(int V, double alpha, double target)
        : TreeSplitter(V), alpha(alpha), target(target) {
        if(alpha < 0.0) throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;
    double target;


    virtual double compute_unnormalized_edge_cut_weight(
        EdgeCut const &edge_cut
    ) const override;
};


class ExperimentalSplitter : public TreeSplitter{

public:

    ExperimentalSplitter(int V, double epsilon, double target)
        : TreeSplitter(V), epsilon(epsilon), target(target) {
        if(epsilon < 0.0) throw Rcpp::exception("Epsilon must be greater than zero!");
    }

    double epsilon;
    double target;

    std::pair<bool, EdgeCut> select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> const &valid_edges,
        bool save_selection_prob
    ) const override;

    double get_log_selection_prob(
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const override;
};


#endif
