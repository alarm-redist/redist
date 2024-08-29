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

    std::vector<std::string> region_labels; // Representation of regions in string form (easy to trace lineage)
    // Regions have R prefix whereas districts are just an integer (in string form). This is a vector length V
    // which should have num_regions district values (I DO NOT CHECK THIS RIGHT NOW)
    std::vector<int> region_num_ids; // Representation of regions in integer form (not as easy to trace
    // lineage). This is a length V vector where ith value maps it the integer id of its region
    std::vector<int> region_dval; //Vector of length V mapping vertex to the d value of its region
    std::vector<double> region_pop; //Vector of length V mapping vertex to the population of ITS REGION

    std::map<std::string, int> r_to_d_map; // Map of region label values to the d_n,k value (number of districts
    std::map<std::string, double> r_to_pop_map; // Map of region label values to The population of that region

    std::map<std::string, int> str_label_to_num_id_map; //Maps region label values to integer id number
    std::map<int, std::string> num_id_to_str_label_map; //Maps region label values to integer id number


    // methods
    void Rprint() const;

};

#endif
