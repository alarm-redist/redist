#pragma once
#ifndef BASE_PLAN_TYPE_H
#define BASE_PLAN_TYPE_H

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#include <vector>
#include <iostream>
#include <RcppThread.h>
#include <RcppArmadillo.h>
#include <string>
#include <cassert>
#include <map>

#include "gredist_types.h"
#include "splitting_schedule_types.h"
#include "tree_splitter_types.h"
#include "tree_op.h"
#include "wilson.h"

// [[Rcpp::depends(RcppArmadillo)]]

// forward declaration for compilation
class TreeSplitter;

class Plan {

private:
    // checks inputted plan is valid
    void check_inputted_region_ids(); 
    void check_inputted_region_sizes(bool split_district_only);

public:
    // constructor 
    Plan(arma::subview_col<arma::uword> region_ids_col, 
        arma::subview_col<arma::uword> region_sizes_col, 
        int ndists, int num_regions, const arma::uvec &pop, bool split_district_only
    ); // constructor for plan

    void shallow_copy(const Plan& plan_to_copy); // Shallow copies everything and doesn't change the underlying 
    // umat columns that region_ids and region_sizes

    // Virtual deep_clone method to create a copy of the object
    virtual std::unique_ptr<Plan> deep_clone() const = 0;


    // attributes
    int ndists; // Number all d_nk must sum to
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
    arma::subview_col<arma::uword> region_sizes; //Vector of length num_regions mapping region ids to their d values
    // Regions have R prefix whereas districts end in an integer (in string form).
    std::vector<int> region_pops; // Vector of length num_regions mapping region ids
    // to the regions population

    std::vector<int> region_added_order; // Vector of length num_regions that
    // tracks the relative order that regions were added. As in if we have
    // region_added_order[i] > region_added_order[j] then that means region i
    // was added after region j
    
    int region_order_max; // value of current largest region order number


    // not actually used by all plan types
    Graph forest_graph;    


    virtual ~Plan() = default; 

    // methods
    void Rprint() const;
    void reorder_plan_by_oldest_split(Plan &dummy_plan);
    virtual Graph get_forest_adj(){throw Rcpp::exception("Get Forest Adj not Supported for this!\n");};

    // redist_smc related methods 
    double choose_multidistrict_to_split(int &region_id_to_split, std::vector<bool> const &valid_region_sizes_to_split);
    bool draw_tree_on_region(const MapParams &map_params, const int region_to_draw_tree_on,
        Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore, int &root);


    void update_region_info_from_cut(
        EdgeCut cut_edge, bool split_district_only,
        const int split_region1_id, const int split_region2_id
    );

    virtual void update_vertex_info_from_cut(
        Tree &ust, EdgeCut cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool split_district_only
    ) = 0;

    

    // virtual redist_smc methods

    bool attempt_split(const MapParams &map_params, const SplittingSchedule &splitting_schedule,
                Tree &ust, TreeSplitter &tree_splitter,
                 std::vector<bool> &visited, std::vector<bool> &ignore, 
                 int const min_region_cut_size, int const max_region_cut_size, 
                 std::vector<int> const &smaller_cut_sizes_to_try,
                 const bool split_district_only, 
                 const int region_id_to_split, const int new_region_id);

};



// Custom hash function for hashing pairs of regions 
// ndists should be the number of regions minus 1, ie the biggest
// region id
struct bounded_hash {
    int ndists;
    bounded_hash(int max_value) : ndists(max_value) {}
    std::size_t operator()(const std::pair<int, int>& p) const {
        return p.first * (ndists + 1) + p.second;
    }
};

#endif