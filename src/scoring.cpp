/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2025/02
 * Purpose: Implement constraint and scoring function
 ********************************************************/

 #include "scoring.h"




 double compute_log_pop_temper(
    double const target, double const pop_temper, int const ndists,
    int const region_pop, int const region_size
){
    // get population deviation 
    double const pop_dev = std::fabs(
        region_pop/static_cast<double>(region_size) - target
    )/target;

    double const pop_pen = std::sqrt(static_cast<double>(ndists) - 2) * std::log(1e-12 + pop_dev);

    // now return the values for the old region minus the two new ones
    return pop_pen * pop_temper;
}


double PopTemperConstraint::compute_log_region_constraint(const Plan &plan, int const region_id) const{
    int const region_pop = plan.region_pops[region_id];
    int const region_size = plan.region_sizes(region_id);
    return compute_log_pop_temper(target, pop_temper, ndists, region_pop, region_size);
}


double PopTemperConstraint::compute_log_merged_region_constraint(const Plan &plan, int const region1_id, int const region2_id) const{
    int const merged_region_pop = plan.region_pops[region1_id] + plan.region_pops[region2_id];
    int const merged_region_size = plan.region_sizes(region1_id) + plan.region_sizes(region2_id);
    return compute_log_pop_temper(target, pop_temper, ndists, merged_region_pop, merged_region_size);
}







// Scoring function
ScoringFunction::ScoringFunction(
    MapParams const &map_params,
    Rcpp::List const &constraints, bool const do_pop_temper, double const pop_temper
){
    // add pop temper if doing that 
    if(do_pop_temper){
        non_final_plan_constraint_ptrs.emplace_back(
            std::make_unique<PopTemperConstraint>(
                map_params.target, map_params.ndists, 
                map_params, pop_temper
            ));
    }
}



double ScoringFunction::compute_log_region_score(const Plan &plan, int const region_id, bool const is_final) const{
    // start out with log score of zero
    double log_score = 0.0;

    // get the log constraint 
    for(auto const &constraint_ptr: constraint_ptrs){
        log_score += constraint_ptr->compute_log_region_constraint(plan, region_id);
    }
    // if final then return 
    if(is_final) return log_score;
    
    // if not final then add those constraints 
    for(auto const &constraint_ptr: non_final_plan_constraint_ptrs){
        log_score += constraint_ptr->compute_log_region_constraint(plan, region_id);
    }

    return log_score;
}