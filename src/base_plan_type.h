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
#include "map_calc.h"
#include "ust_sampler.h"
#include "graph_ops.h"

// [[Rcpp::depends(RcppArmadillo)]]




/*
 * Custom Hash for hashing a pair of region ids 
 *
 * This class implements a hash function for a pair a region ids. 
 * Using the fact that all pair ids will be between 0 and 
 * `max_num_regions-1` (inclusive) it uses a simple grid id hash. 
 */
struct bounded_hash {
    int max_num_regions;
    bounded_hash(int max_num_regions) : max_num_regions(max_num_regions) {}
    std::size_t operator()(const std::pair<int, int>& p) const {
        return p.first * (max_num_regions + 1) + p.second;
    }
    // DANGEROUS DO NOT ACTUALLY USE!!
    // Should probably remove in the future, just needed it to be able to pass default 
    // version to boundary length function
    bounded_hash() : max_num_regions(0) {}
};

class USTSampler;

/*
 * Abstract Class implementation of plan 
 * 
 * The `Plan` object is the abstract class implementation of a plan. The 
 * abstract class encapsulates all the attributes and methods a plan needs
 * to be able to do regardless of sampling space. 
 * 
 * Attributes:
 *  - num_regions The number of regions the plan currently has 
 *  - region_ids A vector of length V mapping each vertex of the graph to the id
 *      of the region it is associated with. 
 *  - region_sizes A vector of length `ndists` storing the size of each region 
 *      (indexed by region id). The sum of all entries should always equal `ndists`.
 *  - region_pops A vector of length `ndists` storing the population of each region 
 *    (indexed by region id). 
 *  - region_added_order A vector of length `ndists` used to deduce the relative order
 *    of which regions were most recently split in the plan. The relative ordering of
 *    the values in this vector indicate which region was split more recently. e.g. 
 *    if `region_added_order[i] < region_added_order[j]` then that means region j 
 *    was split more recently than region i. 
 * 
 * Methods:
 * 
 * 
 * Abstract Methods: 
 */
class Plan {

private:
    // checks inputted plan is valid
    void check_inputted_region_ids(int ndists) const;
    void check_inputted_region_sizes(int ndists, bool split_district_only) const;

protected:
    VertexGraph forest_graph; 
    std::vector<std::tuple<int, int, double>> linking_edges;

public:
    // constructor for a blank plan 
    Plan(int const V, int const ndists, int const nsim_number,
        int const total_pop,
        AllPlansVector &all_plans_vec, 
        AllRegionSizesVector &all_region_sizes_vec,
        std::vector<int> &all_region_pops_vec,
        std::vector<int> &all_region_order_added_vec
   );
   // constructor for partial plan (more than 1 region)
    Plan(int const V, int const ndists, int const num_regions,
        int const nsim_number, const arma::uvec &pop,
        AllPlansVector &all_plans_vec, 
        AllRegionSizesVector &all_region_sizes_vec,
        std::vector<int> &all_region_pops_vec,
        std::vector<int> &all_region_order_added_vec
    ); 

    // shallow copy methods 
    void shallow_copy(Plan const &plan_to_copy);

    // attributes
    PlanVector region_ids;
    RegionSizes region_sizes;
    IntPlanAttribute region_pops;
    IntPlanAttribute region_added_order;
    int num_regions; // Number of regions in the plan
    int region_order_max; 


    
    virtual ~Plan() = default; 

    // methods
    void Rprint() const;
    void reorder_plan_by_oldest_split(Plan &dummy_plan);
    std::pair<int, int> get_most_recently_split_regions() const;
    std::pair<int, int> get_num_district_and_multidistricts() const;
 
    virtual VertexGraph get_forest_adj(){throw Rcpp::exception("Get Forest Adj not Supported for this!\n");};

    virtual std::vector<std::array<double, 3>> get_linking_edges(){
        throw Rcpp::exception("Get Linking edges not Supported for this concrete Plan class!\n");
    };

    // Compute the log number of spanning trees on a region 
    double compute_log_region_spanning_trees(MapParams const &map_params,
        int const region_id) const;

    // Compute the log number of spanning trees on a merged region 
    double compute_log_merged_region_spanning_trees(MapParams const &map_params,
        int const region1_id, int const region2_id) const;

    // Compute the log number of spanning trees on the entire plan
    double compute_log_plan_spanning_trees(MapParams const &map_params) const;

    double compute_log_linking_edge_count(MapParams const &map_params) const;

    // count global number of county splits the map creates
    int count_county_splits(MapParams const &map_params, std::vector<bool> &visited) const;

    int count_merged_county_splits(MapParams const &map_params, std::vector<bool> &visited,
        int const region1_id, int const region2_id) const;

    std::pair<int, std::vector<std::pair<int,int>>> get_or_count_valid_adj_regions_ignore_counties(
        MapParams const &map_params, SplittingSchedule const &splitting_schedule,
        bool const count_only = false
    ) const;

    // Count the number of valid adj regions in a map
    virtual int count_valid_adj_regions(
        MapParams const &map_params, SplittingSchedule const &splitting_schedule
    ) const;

    // Get a vector of all valid adj region pairs
    virtual std::vector<std::pair<int,int>> get_valid_adj_regions(
        MapParams const &map_params, SplittingSchedule const &splitting_schedule,
        bool const check_split_constraint = true
    ) const;


    // redist_smc related methods 
    int choose_multidistrict_to_split(std::vector<bool> const &valid_region_sizes_to_split,
        RNGState &rng_state,
        double const selection_alpha = SELECTION_ALPHA) const;

    bool draw_tree_on_region(const MapParams &map_params, const int region_to_draw_tree_on,
        Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore, int &root,
        RNGState &rng_state);


    void update_region_info_from_cut(
        EdgeCut cut_edge,
        const int split_region1_id, const int split_region2_id,
        bool const add_region
    );


    void update_from_successful_split(
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler, EdgeCut const &cut_edge,
        int const new_region1_id, int const new_region2_id,
        bool const add_region
    );

    // virtual redist_smc methods
    virtual void update_vertex_and_plan_specific_info_from_cut(
        TreeSplitter const &tree_splitter,
        USTSampler &ust_sampler, EdgeCut const cut_edge, 
        const int split_region1_id, const int split_region2_id,
        bool const add_region
    ) = 0;

    // Computes the log effective boundary length between two regions
    // The specifics depend on the sampling space
    virtual double get_log_eff_boundary_len(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const = 0;

    // For a given plan and splitting schedule this finds all pairs of adjacent regions
    // and the log eff boundary length
    // - for graph sampling its just the log of the graph theoretic boundary legnth
    // - for forest sampling its the effective tree boundary length
    virtual std::vector<std::tuple<int, int, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter,
        std::unordered_map<std::pair<int, int>, double, bounded_hash> const &existing_pair_map = {}
    ) const = 0;  
};



#endif