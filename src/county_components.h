#pragma once
#ifndef COUNTY_COMPONENT_H
#define COUNTY_COMPONENT_H


#include <vector>
#include <RcppArmadillo.h>
#include "gredist_types.h"

// [[Rcpp::depends(RcppArmadillo)]]



class CountyComponents{

public:

    CountyComponents(
        MapParams const &map_params, int const num_regions
    );

    MapParams const &map_params;
    int num_regions;
    int max_possible_num_componets; // This is just num_regions-1 + num_counties
    std::vector<CountyRegion> county_district_lookup_table; // lookup table 
    // indexed by (region_num, county_num-1) 
    std::vector<bool> vertices_visited;
    std::vector<bool> component_pairs_visited;
    
    CountyComponentGraph county_component_graph;
    CountyComponentGraph county_component_forest;
    std::vector<CountyRegion> region_component_counts;
    int num_components; // number of connected components 


    // This builds the county component graph and tree 
    void build_component_graph_and_tree(PlanVector const &region_ids);
    
};

#endif
