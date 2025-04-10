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

    RcppThread::ThreadPool pool(num_threads);

    pool.parallelFor(0, num_plans, [&] (int i) {
        static thread_local std::vector<bool> visited(map_params.V);
        // check number of counties is valid
        // if too many then log(0) = -Inf
        if(plans_vec[i].count_county_splits(map_params, visited) > num_regions - 1){
            log_unnormalized_density(i) = -arma::math::inf();
        }else{
            // compute the tau for the plan
            log_unnormalized_density(i) += rho * plans_vec[i].compute_log_plan_spanning_trees(map_params);
            
            // subtract score from each region
            for (size_t region_id = 0; region_id < num_regions; region_id++)
            {
                log_unnormalized_density(i) -= scoring_function.compute_region_score(
                    plans_vec[i], region_id, is_final
                );
            }
        }
    });

    pool.wait();

    return log_unnormalized_density;
}




arma::vec compute_plans_log_optimal_weights(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double pop_temper,  double rho,
    std::string const &splitting_schedule_str,
    int ndists, int num_regions,
    double lower, double target, double upper,
    arma::umat region_ids, arma::umat region_sizes,
    int num_threads
){
    // Create the plan objects
    int num_plans = region_ids.n_cols;

    // sampling space
    SamplingSpace sampling_space = SamplingSpace::GraphSpace;
    
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
    // get splitting schedule 
    Rcpp::List control;
    SplittingSizeScheduleType splitting_schedule_type = get_splitting_size_regime(splitting_schedule_str);
    auto splitting_schedule_ptr = get_splitting_schedule(
        1, ndists, splitting_schedule_type, control
    );
    splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
        0, num_regions-1
    );

    

    // create the splitter
    NaiveTopKSplitter tree_splitter(map_params.V, 1);

    bool compute_log_splitting_prob = splitting_schedule_ptr->schedule_type != SplittingSizeScheduleType::DistrictOnly &&
    plans_vec[0].num_regions != ndists;
    
    // create thread pool
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 1) num_threads = 0;

    RcppThread::ThreadPool pool(num_threads);

    arma::vec log_weights(num_plans, arma::fill::none);

    const int nsims = static_cast<int>(plans_vec.size());
    const int check_int = 50; // check for interrupts every _ iterations

    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        // REprintf("I=%d\n", i);
        log_weights(i) = compute_log_optimal_weights(
            map_params, *splitting_schedule_ptr, sampling_space,
            scoring_function, rho,
            plans_vec.at(i), 
            tree_splitter,
            compute_log_splitting_prob,
            is_final
        );

        ++bar;

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();
    return log_weights;
}
