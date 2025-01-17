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
    virtual ~TreeSplitter() = default; 
    
    // Takes a spanning tree and returns the edge to cut if successful
    virtual std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const int region_id_to_split) = 0;

    // used to update the k parameter for top k splitter
    virtual void update_single_int_param(int int_param){
        throw Rcpp::exception("Update single int param not implemented!\n");
    };

};


// Splitting method that just tries to pick one of the top k edges unif 
class NaiveTopKSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    NaiveTopKSplitter(int k_param): k_param(k_param) {}
    

    int k_param; // top k value

    // how to update the k param
    void update_single_int_param(int int_param);


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);
};



// Splitting method that just tries to pick one of the top k edges unif 
class UniformValidSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    UniformValidSplitter() = default;


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);
};


// Splitting method that picks edge w/ prob ∝ exp(-alpha*bigger dev)
class ExpoWeightedSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    ExpoWeightedSplitter(double alpha): alpha(alpha){
        if(alpha < 0.0) throw Rcpp::exception("Alpha must be less than zero!");
    }

    double alpha;


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params, Plan &plan,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size, 
        const int region_id_to_split);
};


#endif