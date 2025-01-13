#pragma once
#ifndef TREE_SPLITTER_TYPES_H
#define TREE_SPLITTER_TYPES_H

#include <RcppArmadillo.h>
#include "gredist_types.h"
#include "tree_op.h"
#include "tree_splitting.h"

// [[Rcpp::depends(RcppArmadillo)]]



// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut 
class TreeSplitter {

public:
    // Default Constructor 
    virtual ~TreeSplitter() = default; 
    
    // Takes a spanning tree and returns the edge to cut if successful
    virtual std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const arma::subview_col<arma::uword> &region_ids, 
        const int region_id_to_split, const int total_region_pop, const int total_region_size) = 0;

};


// Splitting method that just tries to pick one of the top k edges unif 
class NaiveTopKSplitter : public TreeSplitter{

public:
    // implementation of the pure virtual function
    NaiveTopKSplitter(int k_param): k_param(k_param) {}
    

    int k_param; // top k value


    std::pair<bool,EdgeCut> select_edge_to_cut(
        const MapParams &map_params,
        Tree &ust, const int root, 
        const int min_potential_cut_size, const int max_potential_cut_size,
        const arma::subview_col<arma::uword> &region_ids, 
        const int region_id_to_split, const int total_region_pop, const int total_region_size);
};


#endif