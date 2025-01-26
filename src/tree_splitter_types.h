#pragma once
#ifndef TREE_SPLITTER_TYPES_H
#define TREE_SPLITTER_TYPES_H

#include <RcppArmadillo.h>
#include "gredist_types.h"
#include "tree_op.h"
#include "tree_splitting.h"
#include "base_plan_type.h"


// [[Rcpp::depends(RcppArmadillo)]]

class Plan;

// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut 
class TreeSplitter {

public:
    // Default Constructor 
    TreeSplitter(int V) : pops_below_vertex(V, 0) {};
    virtual ~TreeSplitter() = default; 

    // for storing population below each node
    std::vector<int> pops_below_vertex;

    // Returns a vector of all the valid edges in the tree 
    std::vector<EdgeCut> get_all_valid_pop_edge_cuts_in_directed_tree(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const int region_id_to_split
    );
    
    // Takes a spanning tree and returns the edge to cut if successful
    virtual std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const int region_id_to_split) = 0;

    // Takes a vector of valid edge cuts and returns the log probability 
    // the one an index idx would have been chosen 
    virtual double get_log_selection_prob(
        const MapParams &map_params,
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const = 0;

    // used to update the k parameter for top k splitter
    virtual void update_single_int_param(int int_param){
        throw Rcpp::exception("Update single int param not implemented!\n");
    };

};


// Splitting method that just tries to pick one of the top k edges unif 
class NaiveTopKSplitter : public TreeSplitter{

public:
    // Constructor for NaiveTopKSplitter
    NaiveTopKSplitter(int V, int k_param)
        : TreeSplitter(V), k_param(k_param), vertex_parents(V,0) {
        std::random_device rd;  // A seed source for the random number engine
        gen = std::mt19937(rd()); // Mersenne Twister engine seeded with rd()
    }

    // Attributes specific to NaiveTopKSplitter
    int k_param;        // Top k value
    std::vector<int> vertex_parents;
    std::mt19937 gen;   // Random number generator

    // how to update the k param
    void update_single_int_param(int int_param);


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);

    double get_log_selection_prob(
        const MapParams &map_params,
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const {throw Rcpp::exception("No log selection prob implemented for naive k!\n");}
};



// Splitting method that just tries to pick one of the top k edges unif 
class UniformValidSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    UniformValidSplitter(int V): TreeSplitter(V){};


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);

    // since uniform log prob is just -log(# of candidates)
    double get_log_selection_prob(
        const MapParams &map_params,
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const {return -std::log(valid_edges.size());}
};


// Splitting method that picks edge w/ prob ∝ exp(-alpha*bigger dev)
class ExpoWeightedSplitter : public TreeSplitter{

public:

    ExpoWeightedSplitter(int V, double alpha)
        : TreeSplitter(V), alpha(alpha) {
        if(alpha < 0.0) throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);

    double get_log_selection_prob(
        const MapParams &map_params,
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const;
};



class ExpoWeightedSmallerDevSplitter : public TreeSplitter{

public:

    ExpoWeightedSmallerDevSplitter(int V, double alpha)
        : TreeSplitter(V), alpha(alpha) {
        if(alpha < 0.0) throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);

    double get_log_selection_prob(
        const MapParams &map_params,
        const std::vector<EdgeCut> &valid_edges,
        int idx
    ) const;
};


#endif