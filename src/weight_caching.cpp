/********************************************************
* Author: Philip O'Sullivan
* Institution: Harvard University
* Date Created: 2025/12
* Purpose: Class for saving expsensive computations 
********************************************************/


#include "weight_caching.h"



WeightCacheEnsemble::WeightCacheEnsemble(
    bool const using_caching,
    MapParams const &map_params,
    int const nsims, double const rho,
    SamplingSpace const sampling_space
):
    using_caching(using_caching),
    ndists(map_params.ndists),
    rho(rho),
    nsims(nsims),
    flattened_all_region_order_added(using_caching ? map_params.ndists*nsims : 0, -2),
    flattened_all_cached_weight_values(using_caching ? map_params.ndists*nsims : 0, 0),
    weight_cache_ptr_vec(using_caching ? nsims : 0)
{
    if(!using_caching){
        return;
    }

    weight_cache_ptr_vec.reserve(nsims);
    for (size_t i = 0; i < nsims; i++)
    {
        DoublePlanAttribute weight_cache_values(flattened_all_cached_weight_values, ndists * i, ndists * (i+1));
        IntPlanAttribute cache_region_order_added(flattened_all_region_order_added, ndists * i, ndists * (i+1));
        
        weight_cache_ptr_vec[i] = std::make_unique<WeightCache>(
            ndists, 
            cache_region_order_added,
            weight_cache_values
        );
    }
    
};

// this computes rho-1 compactness
double compute_or_fetch_log_region_compactness(
    Plan const &plan, RegionID const region_id,
    MapParams const &map_params,
    double const rho, 
    bool const using_caching, WeightCache *weight_cache
){
    // If there are 
    // if caching and 
    // region order added number is the same for both then value is still fresh so return that
    if(using_caching && 
        plan.region_added_order[region_id] == weight_cache->this_plan_order_added[region_id]
    ){
        return weight_cache->region_cache_values[region_id];
    }
    // else we need to compute the values 
    double compactness_term = (rho-1) * plan.compute_log_region_spanning_trees(
        map_params, region_id
    );

    if(using_caching){
        // store the new value and update the counter 
        weight_cache->this_plan_order_added[region_id] = plan.region_added_order[region_id];
        weight_cache->region_cache_values[region_id] = compactness_term;
    }

    return compactness_term;

}