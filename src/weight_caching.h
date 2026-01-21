#pragma once
#ifndef WEIGHTS_CACHING_H
#define WEIGHTS_CACHING_H

#include "redist_types.h"
#include "base_plan_type.h"

// Class used to save expensive computations like compactness for specific plans 
class WeightCache {

    public:
        WeightCache(
            int const ndists,
            IntPlanAttribute &this_plan_order_added,
            DoublePlanAttribute &region_cache_values
        ):
            this_plan_order_added(this_plan_order_added),
            region_cache_values(region_cache_values)
            // pair_map(ndists*3 - 6)
            {};

    IntPlanAttribute this_plan_order_added; // This stores the region order counter when region value was cached
    DoublePlanAttribute region_cache_values; // This stores compactness + constraint score for each region 
    // std::unordered_set<std::pair<RegionID, RegionID>, std::pair<int, double>> pair_map;
    // This stores compactness + constraint score for region pairs 

    // This stores compactness + constraint score for each region 
    
    
    // Copy cached data from another WeightCache
    void copy_from(const WeightCache &other) {
        this_plan_order_added.copy(other.this_plan_order_added);
        region_cache_values.copy(other.region_cache_values);
    }

    // For a specific region it fetches the log compactness + score if its still valid
    // or if stale it returns false 
    std::pair<bool, double> attempt_to_get_region_value(
        RegionID const region_id, 
        int const current_region_order_added_num
    );

};


// Wrapper for all the WeightCache's
class WeightCacheEnsemble {

    public:

        WeightCacheEnsemble(
            bool const using_caching,
            MapParams const &map_params,
            int const nsims, double const rho,
            SamplingSpace const sampling_space
        );

    bool const using_caching;
    int const ndists;
    double const rho; // compactness 
    int const nsims;
    std::vector<int> flattened_all_region_order_added;
    std::vector<double> flattened_all_cached_weight_values;
    std::vector<std::unique_ptr<WeightCache>> weight_cache_ptr_vec;

    // This stores compactness + constraint score for each region 
    // This stores compactness + constraint score for region pairs 
    // This is the region update counter 

};


// For a region 
// If using caching it either retrieves old computed value if still fresh
// or computes new one and stores it if stale
// If no caching then just computes the value 
double compute_or_fetch_log_region_compactness(
    Plan const &plan, RegionID const region_id,
    MapParams const &map_params,
    double const rho, 
    bool const using_caching, WeightCache *weight_cache
);



#endif
