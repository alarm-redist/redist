#pragma once
#ifndef COUNTY_COMPONENT_H
#define COUNTY_COMPONENT_H


#include <vector>
#include <RcppArmadillo.h>
#include "gredist_types.h"
#include "base_plan_type.h"

// [[Rcpp::depends(RcppArmadillo)]]

class Plan;

class CountyComponents{

public:

    CountyComponents(
        MapParams const &map_params, int const num_regions
    );

    MapParams const &map_params;
    bool const counties_on;
    int const num_regions;
    int const max_possible_num_componets; // This is just num_regions-1 + num_counties
    std::vector<CountyRegion> county_district_lookup_table; // lookup table 
    // indexed by (region_num, county_num-1) 
    std::vector<bool> vertices_visited;
    std::vector<bool> components_visited;
    std::vector<bool> merged_components_visited;
    std::vector<bool> component_pairs_visited;
    std::vector<bool> counties_component_adj; // used for checking if there's more than one edge between counties on the components graph
    std::vector<bool> neighbor_regions_checked; // used for checking if more than one county overlap with neighbors when checking merge
    
    CountyComponentGraph county_component_graph;
    std::vector<CountyComponentVertex> region_vertices; // gives the vertex of some component of each region
    std::vector<CountyRegion> region_component_counts;
    std::vector<CountyComponentVertex> component_graph_vertices; // maps component graph index to the actual vertex

    DistinctPairHash<RegionID, CountyID> shared_county_hash; // hash for which county (if any) regions overlap in
    
    int num_components; // number of connected components 

    // resets everything 
    void reset();


    // This tries building the county component graph and tree for a hierarchical plan
    // If any region intersect county has more than one connected component or there
    // are more than num_region-1 components total then this will 
    // return false and stop building the graph 
    bool build_component_graph(PlanVector const &region_ids);

    // counts splits in a plan and checks if any county intersect district has more than 1 
    // connected component 
    std::pair<bool, int> count_county_splits(Plan const &plan);

    // given the number of splits and components is ok 
    // checks the component graph reduced to a county multigraph is ok
    bool check_is_county_component_multigraph_valid(Graph &county_graph);

    // checks if a plan is possible to generate hiearchically with 
    // respect to the given counties 
    bool check_valid_hiearchical_plan(Plan const &plan, 
        Graph &county_graph);

    // checks if two adjacent regions can be merged
    bool check_merging_regions_is_ok(
        int const region1_id, int const region2_id
    );
    
    bool count_county_boundary(
        RegionID region1_id, CountyID county_A,
        RegionID region2_id, CountyID county_B
    ) const;
};

#endif
