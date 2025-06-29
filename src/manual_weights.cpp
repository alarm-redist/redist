/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes weight calculation functions to R code
********************************************************/

bool DEBUG_MANUAL_WEIGHTS_VERBOSE = false;

#include "manual_weights.h"

// TODO need to add checks for 
// - population tolerance is ok
// - plans are connected
Rcpp::NumericVector compute_log_unnormalized_plan_target_density(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double const pop_temper,  double const rho,
    int const ndists, int const total_seats, int const num_regions,
    Rcpp::IntegerVector const &district_seat_sizes,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int num_threads
){
    // create the map param object
    MapParams map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper);
    // Add hard constraints to scoring function 
    constraints["plan_valid_district_sizes"] = true;
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);
    // Create the plan objects
    int num_plans = region_ids.ncol();

    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);

    PlanEnsemble plan_ensemble(
        map_params, num_regions,
        num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;
    Rcpp::NumericVector log_unnormalized_density(num_plans);

    const int check_int = 50; // check for interrupts every _ iterations
    Rcpp::Rcout << "Computing Log Target Density!" << std::endl;
    RcppThread::ProgressBar bar(num_plans, 1);
    pool.parallelFor(0, num_plans, [&] (int i) {
        static thread_local PlanMultigraph plan_multigraph(map_params);
        static thread_local std::vector<bool> county_component_lookup(
            num_regions * map_params.num_counties, false
        );

        if(DEBUG_MANUAL_WEIGHTS_VERBOSE){
        auto test_result = plan_multigraph.is_hierarchically_connected(
            *plan_ensemble.plan_ptr_vec[i], county_component_lookup
        );
        Rcpp::Rcout << "Result is " << (test_result.first ? "TRUE" : "FALSE") 
                    << " and " << test_result.second << " components!" << std::endl; 
        }

        bool hierarchically_valid = plan_multigraph.is_hierarchically_valid(
            *plan_ensemble.plan_ptr_vec[i], county_component_lookup
        );

        if(DEBUG_MANUAL_WEIGHTS_VERBOSE){
        Rcpp::Rcout << "Result is " << (hierarchically_valid ? "TRUE" : "FALSE") 
                    << std::endl; 
        }


        // If not hierarchically valid then set log(target) = -Inf
        if(!hierarchically_valid){
            log_unnormalized_density[i] = -arma::math::inf();
            ++bar;
            RcppThread::checkUserInterrupt(i % check_int == 0);
            return; // return to break the lambda 
        }
        // Now check hard constraints 
        auto hard_score = scoring_function.compute_hard_plan_constraints_score(*plan_ensemble.plan_ptr_vec[i]);
        if(!hard_score.first){
            log_unnormalized_density[i] = -arma::math::inf();
            ++bar;
            RcppThread::checkUserInterrupt(i % check_int == 0);
            return; // return to break the lambda 
        }

        
        log_unnormalized_density[i] = - hard_score.second;
        // compute the tau for the plan if we care about it
        if(rho != 0){
            log_unnormalized_density[i] += rho * plan_ensemble.plan_ptr_vec[i]->compute_log_plan_spanning_trees(map_params);
        }
        // subtract score from each region
        for (size_t region_id = 0; region_id < num_regions; region_id++)
        {
            log_unnormalized_density[i] -= scoring_function.compute_region_score(
                *plan_ensemble.plan_ptr_vec[i], region_id, is_final
            );
        }
        
        ++bar;
        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    pool.wait();

    return log_unnormalized_density;
}


// returns target density term for each region 
Rcpp::NumericMatrix compute_log_unnormalized_region_target_density(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double const pop_temper,  double const rho,
    int const ndists, int const total_seats, int const num_regions,
    Rcpp::IntegerVector const &district_seat_sizes,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int const num_threads
){
    // create the map param object
    MapParams map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper);
    
    // Add hard constraints to scoring function 
    constraints["plan_valid_district_sizes"] = true;
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);

    // Create the plan objects
    int num_plans = region_ids.ncol();

    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);
    PlanEnsemble plan_ensemble(
        map_params, num_regions,
        num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;
    Rcpp::NumericMatrix log_unnormalized_region_densities(num_regions, num_plans);

    const int check_int = 50; // check for interrupts every _ iterations
    Rcpp::Rcout << "Computing Log Target Density!" << std::endl;
    RcppThread::ProgressBar bar(num_plans, 1);
    pool.parallelFor(0, num_plans, [&] (int i) {
        static thread_local PlanMultigraph plan_multigraph(map_params);
        static thread_local std::vector<bool> county_component_lookup(
            num_regions * map_params.num_counties, false
        );

        bool hierarchically_valid = plan_multigraph.is_hierarchically_valid(
            *plan_ensemble.plan_ptr_vec[i], county_component_lookup
        );



        // auto county_splits_result = county_components.count_county_splits(plans_vec[i]);
        // check number of counties is valid and no double county intersect region components
        // if too many then log(target) = -Inf

        // TODO only check for each region 

        if(!hierarchically_valid){
            for (size_t region_id = 0; region_id < num_regions; region_id++)
            {
                log_unnormalized_region_densities(region_id, i) = -arma::math::inf();
            }
            ++bar;
            RcppThread::checkUserInterrupt(i % check_int == 0);
            return; // return to break the lambda 
        }
        // Now check hard constraints 
        auto hard_score = scoring_function.compute_hard_plan_constraints_score(*plan_ensemble.plan_ptr_vec[i]);
        if(!hard_score.first){
            for (size_t region_id = 0; region_id < num_regions; region_id++)
            {
                log_unnormalized_region_densities(region_id, i) = -arma::math::inf();
            }
            ++bar;
            RcppThread::checkUserInterrupt(i % check_int == 0);
            return; // return to break the lambda 
        }


        // compute contribution for each region
        for (size_t region_id = 0; region_id < num_regions; region_id++)
        {
            log_unnormalized_region_densities(region_id, i) = 0.0;
            if(rho != 0){
                log_unnormalized_region_densities(region_id, i) += rho * plan_ensemble.plan_ptr_vec[i]->compute_log_region_spanning_trees(
                    map_params, region_id
                );
            }
            log_unnormalized_region_densities(region_id, i) -= scoring_function.compute_region_score(
                *plan_ensemble.plan_ptr_vec[i], region_id, is_final
            );
        }
        
        ++bar;
        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    pool.wait();

    return log_unnormalized_region_densities;
}


arma::vec compute_plans_log_optimal_weights(
    List const &adj_list, arma::uvec const &counties, arma::uvec const &pop,
    List const &constraints, double const pop_temper,  double const rho,
    std::string const &splitting_schedule_str,
    int const ndists, int const total_seats, Rcpp::IntegerVector const &district_seat_sizes,
    int const num_regions,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int num_threads
){
    // Create the plan objects
    int num_plans = region_ids.ncol();

    // sampling space
    SamplingSpace sampling_space = SamplingSpace::GraphSpace;
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;

    // create thread pool
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    RcppThread::ThreadPool pool(num_threads);

    // create the map param object
    MapParams map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper);
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);


    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);
    PlanEnsemble plan_ensemble(
        map_params, num_regions,
        num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );

    // get splitting schedule 
    Rcpp::List control;
    SplittingSizeScheduleType splitting_schedule_type = get_splitting_size_regime(splitting_schedule_str);
    auto splitting_schedule_ptr = get_splitting_schedule(
        1, ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        splitting_schedule_type, control
    );
    splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
        0, num_regions-1
    );
    splitting_schedule_ptr->print_current_step_splitting_info();


    // create the splitter
    NaiveTopKSplitter tree_splitter(map_params.V, 1);

    bool compute_log_splitting_prob = splitting_schedule_ptr->schedule_type != SplittingSizeScheduleType::DistrictOnlySMD &&
    plan_ensemble.plan_ptr_vec[0]->num_regions != ndists;
    bool const counties_on = map_params.num_counties > 1;

    
    arma::vec log_weights(num_plans, arma::fill::none);

    const int nsims = plan_ensemble.nsims;
    const int check_int = 50; // check for interrupts every _ iterations

    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        static thread_local PlanMultigraph plan_multigraph(map_params);

        // build the multigraph 
        plan_multigraph.build_plan_multigraph(*plan_ensemble.plan_ptr_vec[i]);
        plan_multigraph.Rprint_detailed(*plan_ensemble.plan_ptr_vec[i]);
        // plan_multigraph.pair_map.Rprint();
        // // remove invalid hierarchical merges
        // plan_multigraph.remove_invalid_hierarchical_merge_pairs(*plan_ensemble.plan_ptr_vec[i]);
        // plan_multigraph.pair_map.Rprint();
        // // remove invalid size pairs 
        // plan_multigraph.remove_invalid_size_pairs(*plan_ensemble.plan_ptr_vec[i], *splitting_schedule_ptr);
        // plan_multigraph.pair_map.Rprint();
        // now make the output vector 
        std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
        region_pairs_tuple_vec.reserve(plan_multigraph.pair_map.num_hashed_pairs);

        double incremental_weight = 0.0;

        for(auto const key_val_pair: plan_multigraph.pair_map.get_all_values()){
            
            double log_boundary_len;
            if(key_val_pair.second.admin_adjacent){
                // if administratively adjacent only count within same county
                log_boundary_len = std::log(key_val_pair.second.within_county_edges);
            }else{
                // else only count across county 
                log_boundary_len = std::log(key_val_pair.second.across_county_edges);
            }

            incremental_weight += log_boundary_len;
            REprintf("(%u, %u) - %.4f \n", 
            key_val_pair.first.first, 
                key_val_pair.first.second, 
                std::exp(log_boundary_len));
        }

        // REprintf("I=%d\n", i);
        // log_weights(i) = compute_log_optimal_incremental_weights(
        //     *plan_ensemble.plan_ptr_vec[i], plan_multigraph,
        //     *splitting_schedule_ptr, tree_splitter,
        //     sampling_space, scoring_function, 
        //     rho, false, is_final
        // );


        REprintf("%f vs %f \n", 
            1.0 / static_cast<double>(incremental_weight), 
            std::exp(log_weights(i))
        );

        ++bar;

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();
    return log_weights;
}
