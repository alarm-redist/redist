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

#include "redist_types.h"
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
class PlanMultigraph;



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
    Plan(int const total_seats,
        int const total_pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added
   );
   // constructor for partial plan (more than 1 region)
    Plan(int const num_regions,
        const arma::uvec &pop,
        PlanVector &this_plan_region_ids, 
        RegionSizes &this_plan_region_sizes,
        IntPlanAttribute &this_plan_region_pops,
        IntPlanAttribute &this_plan_order_added
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
    virtual void Rprint(bool verbose = false) const;
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

    double compute_log_linking_edge_count(PlanMultigraph &plan_multigraph) const;

    // count global number of county splits the map creates
    int count_county_splits(MapParams const &map_params, std::vector<bool> &visited) const;

    int count_merged_county_splits(MapParams const &map_params, std::vector<bool> &visited,
        int const region1_id, int const region2_id) const;

    
    // attempts to build a plan multigraph and return valid merge split pairs 
    virtual std::pair<bool, std::vector<std::pair<RegionID,RegionID>>> attempt_to_get_valid_mergesplit_pairs(
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
    ) const;

    virtual std::vector<std::pair<RegionID,RegionID>> get_valid_smc_merge_regions(
        PlanMultigraph &plan_multigraph, SplittingSchedule const &splitting_schedule
    ) const;


    // redist_smc related methods 
    int choose_multidistrict_to_split(std::vector<bool> const &valid_region_sizes_to_split,
        RNGState &rng_state,
        double const selection_alpha = SELECTION_ALPHA) const;

    std::pair<bool, int> draw_tree_on_region(
        const MapParams &map_params, const int region_to_draw_tree_on,
        Tree &ust, std::vector<bool> &visited, std::vector<bool> &ignore, int &root,
        RNGState &rng_state, int const attempts_to_make
    );


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
    // assumes the muligraph was already generated 
    virtual double get_log_eff_boundary_len(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter, 
        const int region1_id, int const region2_id
    ) const = 0;

    // For a given plan and splitting schedule this finds all pairs of adjacent regions
    // and the log eff boundary length
    // - for graph sampling its just the log of the graph theoretic boundary legnth
    // - for forest sampling its the effective tree boundary length
    virtual std::vector<std::tuple<RegionID, RegionID, double>> get_valid_adj_regions_and_eff_log_boundary_lens(
        PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
        TreeSplitter const &tree_splitter
    ) const = 0;  
};



// simple struct 
// storing data on a pair of adjacent regions 
struct PairHashData{

    // default constructor 
    PairHashData():
        admin_adjacent(false),
        shared_county(-1),
        within_county_edges(0),
        across_county_edges(0),
        merge_is_hier_valid(true),
        count_pair(true),
        eff_boundary_len(0.0)
    {};

    // specific constructor
    PairHashData(
        bool const admin_adjacent, int const shared_county,
        int const within_county_edges, int const across_county_edges,
        bool const merge_is_hier_valid,
        bool const count_pair, double const eff_boundary_len
    ): 
        admin_adjacent(admin_adjacent),
        shared_county(shared_county),
        within_county_edges(within_county_edges),
        across_county_edges(across_county_edges),
        merge_is_hier_valid(merge_is_hier_valid),
        count_pair(count_pair),
        eff_boundary_len(eff_boundary_len)
    {};


    bool admin_adjacent; // whether or not adjacent within a county
    int shared_county; // -1 if the two regions do not share a county, 0 indexed counties otherwise
    int within_county_edges;
    int across_county_edges;
    bool merge_is_hier_valid; // If merging the two regions 
    bool count_pair; // only used for forest space, whether or not to compute probabilties for pairs
    double eff_boundary_len; // NOT THE LOG, THE ACTUAL BOUNDARY
    
};

// Hashes pairs (x,y) where 0 <= x != y < num_elements to the
// values of type U
// under the hood just implements as vector of size num_elements choose 2
class RegionPairHash{

    public:
        // actual constructor
        RegionPairHash(int const ndists):
        ndists(ndists),
        num_hashed_pairs(0),
        num_hier_smc_merge_valid_pairs(0),
        hash_table_size((ndists*(ndists-1))/2),
        values(hash_table_size),
        hashed(hash_table_size, false){
            if(ndists < 2) throw Rcpp::exception("Region Pair hash must have at least 2 elements!\n");
        }

        
        int const ndists;
        int num_hashed_pairs;
        int num_hier_smc_merge_valid_pairs;
        int const hash_table_size;
        // Table where each tuple is 
        // - bool: Whether or not regions are administratively adjacent 
        // - int: Number of edges across county boundary
        // - int: Number of edges within county boundary
        // - bool: Whether or not to count this pair (only relevant for forest space plans)
        // - double: effective boundary (only relevant for forest space plans)
        std::vector<PairHashData> values;


        std::vector<bool> hashed;
        std::vector<std::pair<RegionID, RegionID>> hashed_pairs;


        // Hash function that returns the index of the pair (x, y)
        size_t pair_hash(RegionID region1_id, RegionID region2_id) const{
            // region1_id must always be smaller 
            return ((region1_id * (2 * ndists - region1_id - 1)) / 2 ) + (region2_id - region1_id - 1);
        }

        // resets the hash map to blank state
        void reset(){
            std::fill(
                values.begin(),
                values.end(),
                PairHashData()
            );
            std::fill(
                hashed.begin(),
                hashed.end(),
                false
            );
            hashed_pairs.clear();
            num_hashed_pairs = 0;
            num_hier_smc_merge_valid_pairs = 0;
        }
        

        // gets the value a pair is hashed to if its been hashed
        std::pair<bool, PairHashData> get_value(
            RegionID const x, RegionID const y
        ) const{
            auto hash_index = pair_hash(x, y);
            // check if we've hashed this pair or not
            if(hashed[hash_index]){
                return std::make_pair(true, values[hash_index]);
            }else{
                return std::make_pair(false, PairHashData());
            }
        }

        // void remove_value(
        //     RegionID const region1_id, RegionID const region2_id
        // ){
        //     auto hash_index = pair_hash(region1_id, region2_id);
        //     // do nothing if not in the table
        //     if(!hashed[hash_index]) return;
        //     // else reset the value and decrease the count 
        //     --num_hashed_values;
        //     hashed[hash_index] = false;
        //     values[hash_index] = PairHashData();
        // }

        std::vector<std::pair<
            std::pair<RegionID, RegionID>, PairHashData
        >> get_all_values() const{
            std::vector<std::pair<
            std::pair<RegionID, RegionID>, 
            PairHashData
            >> all_data;
            all_data.reserve(num_hashed_pairs);

            for(auto const a_pair: hashed_pairs){
                auto val = get_value(a_pair.first, a_pair.second);
                all_data.push_back({a_pair, val.second});
            }

            return all_data;
        }



        // increase the graph edge count between a pair of regions 
        void count_graph_edge(
            RegionID region1_id, RegionID region2_id,
            bool const same_county, int const county
        ){
            // swap if region2_id < region1_id
            if(region2_id < region1_id) std::swap(region1_id, region2_id);
            // get the index
            auto hash_index = pair_hash(region1_id, region2_id);
            // add to the pair list if not already there
            if(!hashed[hash_index]){
                hashed[hash_index] = true;
                hashed_pairs.push_back({region1_id, region2_id});
                ++num_hashed_pairs;
            }
            
            if(same_county){
                // if edges in the same county increase within county count
                // and mark as administratively adjacent 
                values[hash_index].admin_adjacent = true;
                values[hash_index].shared_county = county;
                values[hash_index].within_county_edges++;
            }else{
                // else increase count of across county stuff 
                values[hash_index].across_county_edges++;
            }
            return;
        }

        // Add this amount to the effective boundary length
        // only used for forest plans 
        void add_to_eff_boundary(
            RegionID region1_id, RegionID region2_id,
            double const value_to_add
        ){
            // swap if region2_id < region1_id
            if(region2_id < region1_id) std::swap(region1_id, region2_id);
            // get the index
            auto hash_index = pair_hash(region1_id, region2_id);
            // add to the pair list if not already there
            if(!hashed[hash_index]){
                hashed[hash_index] = true;
                hashed_pairs.push_back({region1_id, region2_id});
                ++num_hashed_pairs;
            }
            // increase the count 
            values[hash_index].eff_boundary_len += value_to_add;
            return;
        }

        void Rprint() const;


        RegionMultigraphCount get_multigraph_counts(int const num_regions) const;
        RegionMultigraphCount get_merged_multigraph_counts(
            int const num_regions, std::vector<RegionID> &merge_index_reshuffle,
            RegionID const region1_id, RegionID const region2_id 
        ) const;



};


void swap_pair_maps(RegionPairHash &a, RegionPairHash &b);

// Class for managing the plan region multigraph 
// primarily used for calculating weights and finding adjacent pairs 
class PlanMultigraph{
    public:

        PlanMultigraph(MapParams const &map_params);

        MapParams const &map_params;
        bool const counties_on;
        std::vector<bool> vertices_visited;

        
        std::vector<int> county_component; // maps regions to their county adj component
        std::vector<int> component_split_counts; // number of splits in each county adj component
        std::vector<int> component_region_counts; // number of regions in each county adj component 
        std::vector<std::unordered_set<CountyID>> region_overlap_counties; // Stores which counties a region overlaps in
        
        int num_county_region_components; // number of connected components 

        RegionPairHash pair_map; // used for tracking adj regions 

        // Prints relevant info - for debugging
        void Rprint() const;
        void Rprint_detailed(Plan const &plan);


        // Checks if the a plan is hierarchically connected, ie
        // Each county intersect region has at most 1 component 
        // and if it is then it also counts the number of county region
        // components in the plan 
        std::pair<bool, int> is_hierarchically_connected(
            Plan const &plan, std::vector<bool> component_lookup
        );


        bool build_plan_hierarchical_multigraph(
            Plan const &plan
        );

        // builds the multigraph ignoring counties 
        void build_plan_non_hierarchical_multigraph(
            Plan const &plan
        );


        bool build_plan_multigraph(
            Plan const &plan
        );

        // Removes the pairs with invalid merge sizes 
        // from the hash map 
        void remove_invalid_size_pairs(
            Plan const &plan, const SplittingSchedule &splitting_schedule
        );

        // Removes pairs where merging them would create a hierarchicall
        // invalid plan
        void remove_invalid_hierarchical_merge_pairs(
            Plan const &plan
        );

        // removes pairs where the mergesplit kernel probability is zero
        // Currently this is only 
        //  - regions where merging them results in a non-hierarchically connected 
        void remove_invalid_mergesplit_pairs(
            Plan const &plan
        );

        // Takes an arbitrary plan and checks if its hierarchically valid meaning
        //  - hierarchically connected
        //  - The administratively adjacent plan quotient graph has no cycles
        bool is_hierarchically_valid(
            Plan const &plan, std::vector<bool> component_lookup
        );
};

// swap function 
void swap_plan_multigraphs(PlanMultigraph &a, PlanMultigraph &b);

#endif
