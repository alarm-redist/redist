/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes weight calculation functions to R code
********************************************************/



#include "manual_weights.h"

// TODO need to add checks for 
// - population tolerance is ok
// - plans are connected
Rcpp::NumericVector compute_log_unnormalized_plan_target_density(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double const pop_temper,  double const rho,
    int const ndists, int const num_regions,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int num_threads
){
    // create the map param object
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);
    // Create the plan objects
    int num_plans = region_ids.ncol();
    PlanEnsemble plan_ensemble(
        map_params.V, ndists, num_regions,
        pop, num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes,
        pool 
    );
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;
    Rcpp::NumericVector log_unnormalized_density(num_plans);

    const int check_int = 50; // check for interrupts every _ iterations
    Rcpp::Rcout << "Computing Log Target Density!" << std::endl;
    RcppThread::ProgressBar bar(num_plans, 1);
    pool.parallelFor(0, num_plans, [&] (int i) {
        static thread_local CountyComponents county_components(
            map_params, num_regions
        );
        static thread_local Graph county_graph(map_params.num_counties);

        // auto county_splits_result = county_components.count_county_splits(plans_vec[i]);
        bool const hiearhically_valid = county_components.check_valid_hiearchical_plan(*plan_ensemble.plan_ptr_vec[i], county_graph);
        // check number of counties is valid and no double county intersect region components
        // if too many then log(target) = -Inf
        if(!hiearhically_valid){
            // REprintf("Invalid!\n");
            log_unnormalized_density[i] = -arma::math::inf();
        }else{
            log_unnormalized_density[i] = 0.0;
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
    int const ndists, int const num_regions,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int const num_threads
){
    // create the map param object
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);

    // Create the plan objects
    int num_plans = region_ids.ncol();
    PlanEnsemble plan_ensemble(
        map_params.V, ndists, num_regions,
        pop, num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes,
        pool 
    );
    
    // check if final splits (ie don't do pop_temper)
    bool is_final = num_regions == ndists;
    Rcpp::NumericMatrix log_unnormalized_region_densities(num_regions, num_plans);

    const int check_int = 50; // check for interrupts every _ iterations
    Rcpp::Rcout << "Computing Log Target Density!" << std::endl;
    RcppThread::ProgressBar bar(num_plans, 1);
    pool.parallelFor(0, num_plans, [&] (int i) {
        static thread_local CountyComponents county_components(
            map_params, num_regions
        );
        static thread_local Graph county_graph(map_params.num_counties);

        // auto county_splits_result = county_components.count_county_splits(plans_vec[i]);
        bool const hiearhically_valid = county_components.check_valid_hiearchical_plan(*plan_ensemble.plan_ptr_vec[i], county_graph);
        // check number of counties is valid and no double county intersect region components
        // if too many then log(target) = -Inf
        if(!hiearhically_valid){
            for (size_t region_id = 0; region_id < num_regions; region_id++)
            {
                log_unnormalized_region_densities(region_id, i) = -arma::math::inf();
            }
        }else{
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
        }
        ++bar;
        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    pool.wait();

    return log_unnormalized_region_densities;
}


// std::vector<std::pair<std::pair<int,int>, int>>
std::vector<int> new_boundary_counting(
    MapParams const &map_params,
    std::vector<int> &pair_boundary_count_vec,
    PlanVector const &region_ids, int const num_regions,
    CountyComponents const &county_components
){
    bool const check_counties = map_params.num_counties > 1;

    // iterate through the graph 
    for (int v = 0; v < map_params.V; v++) {
        // Find out which region this vertex corresponds to
        int v_region = region_ids[v];
        int v_county = map_params.counties(v);

        // now iterate over its neighbors
        for (int u : map_params.g[v]) {
            // find which region neighbor corresponds to
            int u_region = region_ids[u];
            // ignore if u_region <= v_region to avoid double counting 
            if(u_region <= v_region) continue;
            if(!check_counties || v_county == map_params.counties(u)){
                // if not checking counties we just always count since
                // incrementing is cheaper then checking and incrementing 
                ++pair_boundary_count_vec[index_from_ordered_pair(v_region, u_region, num_regions)];
                // REprintf("Counting (%d, %d) - Regions (%u, %u) \n", u, v, v_region, u_region);
            }else if(pair_boundary_count_vec[index_from_ordered_pair(v_region, u_region, num_regions)] >= 0){
                // else this edge crosses a county boundary of pair we care about
                // We only count if no edges between 
                //      (v_county, v_region), (u_county, v_region)
                //      (u_county, u_region), (v_county, u_region)
                pair_boundary_count_vec[index_from_ordered_pair(v_region, u_region, num_regions)] += county_components.count_county_boundary(
                    v_region, v_county-1,
                    u_region, map_params.counties(u)-1
                );
                
            }
        }
    }
    return pair_boundary_count_vec;
}


arma::vec compute_plans_log_optimal_weights(
    List const &adj_list, arma::uvec const &counties, arma::uvec const &pop,
    List const &constraints, double const pop_temper,  double const rho,
    std::string const &splitting_schedule_str,
    int const ndists, int const num_regions,
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
    if (num_threads == 1) num_threads = 0;
    RcppThread::ThreadPool pool(num_threads);

    // create the map param object
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    // Create the scoring function 
    ScoringFunction scoring_function(map_params, constraints, pop_temper);

    PlanEnsemble plan_ensemble(
        map_params.V, ndists, num_regions,
        pop, num_plans,
        SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes,
        pool 
    );

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
    plan_ensemble.plan_ptr_vec[0]->num_regions != ndists;
    bool const counties_on = map_params.num_counties > 1;

    


    arma::vec log_weights(num_plans, arma::fill::none);

    const int nsims = plan_ensemble.nsims;
    const int check_int = 50; // check for interrupts every _ iterations

    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        static thread_local CountyComponents county_components(
            map_params, num_regions
        );
        static thread_local std::vector<int> pair_eff_bounary_counts(
            (num_regions*(num_regions-1))/2, -1
        );
        static thread_local EffBoundaryMap pair_map(num_regions, 0.0);
        std::fill(
            pair_eff_bounary_counts.begin(),
            pair_eff_bounary_counts.end(),
            -1
        );
        // first make sure if counties are on then the plan doesn't
        // have an illegal number of splits 
        if(counties_on){
            auto county_splits_result = county_components.count_county_splits(*plan_ensemble.plan_ptr_vec[i]);
            // check number of counties is valid and no double county intersect region components
            // if too many then log(0) = -Inf
            if(county_splits_result.second > num_regions - 1){
                REprintf("Plan %d should have %d county splits or less but it actually has %d!\n",
                    i, num_regions-1, county_splits_result.second);
                throw Rcpp::exception("Plan has too many county splits!\n");
            }else if(county_splits_result.first){
                REprintf("Plan %d has at least one region intersect county with more than 1 connected component!\n",
                    i);
                throw Rcpp::exception("Plan has at least one region intersect county with more than 1 connected component!\n");
            }
            county_components.build_component_graph(plan_ensemble.plan_ptr_vec[i]->region_ids);
        }

        auto adj_pairs_ignoring_county = plan_ensemble.plan_ptr_vec[i]->get_or_count_valid_adj_regions(
                map_params, *splitting_schedule_ptr, county_components, false
        ).second;
        REprintf("Size intially was %u\n", adj_pairs_ignoring_county.size());
        REprintf("Adjacent regions: ");
        for(auto const &a_pair: adj_pairs_ignoring_county){
            REprintf("(%d, %d), ", a_pair.first, a_pair.second);
        }
        REprintf("\n");
        // dangerous but removing invalid pairs while iterating over 
        auto it = adj_pairs_ignoring_county.begin();

        while(it != adj_pairs_ignoring_county.end()){
            bool pair_ok = !county_components.counties_on || county_components.check_merging_regions_is_ok(
                it->first, it->second
            );
            REprintf("Regions (%d, %d): ",
                it->first, it->second
            );
            if(pair_ok) {
                REprintf(" is ok!\n");
                // mark this as pair to check 
                pair_eff_bounary_counts[index_from_ordered_pair(it->first, it->second, num_regions)] = 0;
                REprintf("Index %d, value %d\n",index_from_ordered_pair(it->first, it->second, num_regions),
                     pair_eff_bounary_counts[index_from_ordered_pair(it->first, it->second, num_regions)]);
                // increase iterator 
                ++it;
            }else{
                REprintf(" is NOT ok!\n");
                // erase the pair
                it = adj_pairs_ignoring_county.erase(it);
            }
        }
        REprintf("Size now is %u\n", adj_pairs_ignoring_county.size());


        auto new_results_map = new_boundary_counting(
            map_params, pair_eff_bounary_counts, 
            plan_ensemble.plan_ptr_vec[i]->region_ids, num_regions,
             county_components
        );

        double incremental_weight = 0.0;
        REprintf("\nBoundary Length Maps:\n");
        for(const auto &map_pair: adj_pairs_ignoring_county){
            REprintf("(%d, %d) - %d, \n",
                map_pair.first, map_pair.second,
                pair_eff_bounary_counts[index_from_ordered_pair(map_pair.first, map_pair.second, num_regions)]);
            incremental_weight += pair_eff_bounary_counts[index_from_ordered_pair(map_pair.first, map_pair.second, num_regions)];
        }
        REprintf("\n");

        REprintf("I=%d\n", i);
        log_weights(i) = compute_log_optimal_weights(
            map_params, *splitting_schedule_ptr, sampling_space,
            scoring_function, rho,
            *plan_ensemble.plan_ptr_vec[i], 
            tree_splitter, county_components, pair_map,
            false,
            is_final
        );

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
