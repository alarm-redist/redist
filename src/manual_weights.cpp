/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes weight calculation functions to R code
********************************************************/



#include "manual_weights.h"

arma::vec compute_log_unnormalized_plan_target_density(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double pop_temper,  double rho,
    int ndists, int num_regions,
    double lower, double target, double upper,
    arma::umat region_ids, arma::umat region_sizes,
    int num_threads
){
    // Create the plan objects
    int num_plans = region_ids.n_cols;
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;
    std::vector<GraphPlan> plans_vec; plans_vec.reserve(num_plans);

    for (size_t i = 0; i < num_plans; i++)
    {
        plans_vec.emplace_back(region_ids.col(i), region_sizes.col(i), ndists, num_regions, pop, false);
    }

    // create the map param object
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    arma::vec log_unnormalized_density(num_plans, arma::fill::zeros);
    
    // create thread pool
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 1) num_threads = 0;

    RcppThread::parallelFor(0, num_plans, [&] (int i) {
        // compute the tau for the plan
        log_unnormalized_density(i) += rho * plans_vec[i].compute_log_plan_spanning_trees(map_params);
        
        // subtract score from each region
        for (size_t region_id = 0; region_id < num_regions; region_id++)
        {
            log_unnormalized_density(i) -= scoring_function.compute_region_score(
                plans_vec[i], region_id, is_final
            );
        }
    }, num_threads > 0 ? num_threads : 0);


    return log_unnormalized_density;
}


double compute_a_log_optimal_weight(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List control, 
    int ndists, int num_regions, 
    double lower, double target, double upper,
    arma::umat region_ids, arma::umat region_sizes
){

    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);

    SplittingSizeScheduleType splitting_size_regime = get_splitting_size_regime(
        static_cast<std::string>(control["splitting_size_regime"])
    );

    bool split_district_only = splitting_size_regime == SplittingSizeScheduleType::DistrictOnly;

    auto splitting_schedule_ptr = get_splitting_schedule(
        1, ndists, num_regions, splitting_size_regime, control
    );

    // set the splitting schedule 
    splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
        0, num_regions);

    
    // Create a plan object
    Plan * plan = new GraphPlan(
        region_ids.col(0), region_sizes.col(0), 
        ndists, num_regions, pop, split_district_only);

    
    



    // double compute_log_optimal_weights(
    //     const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    //     ScoringFunction const &scoring_function, Plan const &plan, 
    //     const TreeSplitter &edge_splitter, bool compute_log_splitting_prob,
    //     bool is_final_plan
    // )

    return 3.14;

}