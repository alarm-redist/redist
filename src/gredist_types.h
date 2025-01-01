#pragma once
#ifndef TYPES_H
#define TYPES_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#define PRINT_LN Rcout << __func__ << "(), " << __FILE__ << ":" << __LINE__ << "\n";

#include <vector>
#include <iostream>
#include <RcppThread.h>
#include <RcppArmadillo.h>
#include <string>
#include <cassert>
#include <map>

// [[Rcpp::depends(RcppArmadillo)]]

typedef std::vector<std::vector<int>> Tree;
typedef std::vector<std::vector<int>> Graph;
typedef std::vector<std::vector<std::vector<int>>> Multigraph;



class Plan
{
public:
    // constructor 
    Plan(arma::subview_col<arma::uword> region_ids_col, 
             arma::subview_col<arma::uword> region_dvals_col, 
             int N, int total_map_pop, bool split_district_only=false,
             int num_regions=1, int num_districts=0,
             const arma::uvec &pop = {}); // constructor for 1 region plan
    Plan(const Plan& other); // Copy constructor


    // attributes
    int N; // Number all d_nk must sum to
    int V; // Number of nodes in graph
    int num_regions; // Number of regions in the plan
    int num_districts; // Number of districts in the plan
    int num_multidistricts; // Number of multidistricts, always `num_regions` - `num_districts`
    int map_pop; // The population of the entire map
    int remainder_region; // ID of the remainder region of the plan (if it has one)

    // basically vectors with length V
    arma::subview_col<arma::uword> region_ids; // Representation of regions in integer form (not as easy to trace
    // lineage). This is a length V vector where ith value maps it the integer id of its region

    // Think of this like the pointer to some column in a matrix 


    // vectors with length num_regions
    arma::subview_col<arma::uword> region_dvals; //Vector of length num_regions mapping region ids to their d values
    // Regions have R prefix whereas districts end in an integer (in string form).
    std::vector<int> region_pops; // Vector of length num_regions mapping region ids
    // to the regions population

    std::vector<int> region_added_order; // Vector of length num_regions that
    // tracks the relative order that regions were added. As in if we have
    // region_added_order[i] > region_added_order[j] then that means region i
    // was added after region j
    
    int region_order_max; // value of current largest region order number

    
    // shallow copy another plan
    void shallow_copy(const Plan& other);

    // methods
    void Rprint() const;
    void reorder_plan_by_oldest_split(Plan &dummy_plan);

};



// Custom hash function for hashing pairs of regions 
// N should be the number of regions minus 1, ie the biggest
// region id
struct bounded_hash {
    int N;
    bounded_hash(int max_value) : N(max_value) {}
    std::size_t operator()(const std::pair<int, int>& p) const {
        return p.first * (N + 1) + p.second;
    }
};

#endif
