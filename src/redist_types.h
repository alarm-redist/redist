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
    Plan(int V, int N, double total_map_pop); // constructor for 1 region plan
    // attributes
    int N; // Number all d_nk must sum to
    int V; // Number of nodes in graph
    int num_regions; // Number of regions in the plan
    int num_districts; // Number of districts in the plan
    int num_multidistricts; // Number of multidistricts, always `num_regions` - `num_districts`
    double map_pop; // The population of the entire map

    // vectors with length V
    std::vector<int> region_num_ids; // Representation of regions in integer form (not as easy to trace
    // lineage). This is a length V vector where ith value maps it the integer id of its region


    // vectors with length num_regions
    std::vector<int> region_dvals; //Vector of length num_regions mapping region ids to their d values
    std::vector<std::string> region_str_labels; // Vector of length num_regions mapping region id to string label (if tracked)
    // Regions have R prefix whereas districts end in an integer (in string form).
    std::vector<double> region_pops; // Vector of length num_regions mapping region ids
    // to the regions population

    // methods
    void Rprint() const;

};

#endif
