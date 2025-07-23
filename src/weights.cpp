/********************************************************
* Author: Philip O'Sullivan
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: SMC weight calculation related functions
********************************************************/

constexpr bool DEBUG_WEIGHTS_VERBOSE = false; // Compile-time constant
#include "weights.h"



//' Computes the effective sample size from log incremental weights
//'
//' Takes a vector of log incremental weights and computes the effective sample
//' size which is the sum of the weights squared divided by the sum of squared
//' weights
//'
//'
//' @title Compute Effective Sample Size
//'
//' @param log_wgt vector of log incremental weights
//'
//' @details No modifications to inputs made
//'
//' @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
//'
double compute_n_eff(const arma::subview_col<double> log_wgt) {
    double sum_wgt = 0.0;
    double sum_wgt_squared = 0.0;

    // compute sum of squares and square of sum
    for (const double& log_w : log_wgt) {
        double wgt = std::exp(log_w);
        sum_wgt += wgt;
        sum_wgt_squared += std::exp(2*log_w);
    }


    return std::exp(
        (2 * std::log(sum_wgt)) - std::log(sum_wgt_squared)
    );
}





//' Get the probability the union of two regions was chosen to split
//'
//' Given a plan object and two regions in the plan this returns the probability
//' the union of the two regions was chosen to be split.
//'
//' @title Get Retroactive Split Selection Probability
//'
//' @param plan A plan object
//' @param region1_id The id of the first region to union
//' @param region2_id The id of the second region to union
//'
//' @details No modifications to inputs made
//'
//' @return the log of the probability the union of the two regions would be
//' chosen to split.
//'
double get_log_retroactive_splitting_prob(
        const Plan &plan,
        const std::vector<bool> &valid_region_sizes_to_split,
        const int region1_id, const int region2_id,
        double const selection_alpha = SELECTION_ALPHA
){
    if(DEBUG_WEIGHTS_VERBOSE) Rprintf("Possible options: ");

    // compute weight of the merged region
    auto unioned_region_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];
    double unioned_region_prob = std::pow(unioned_region_size, selection_alpha);
    // get the sum of all regions 
    double prob_sum = unioned_region_prob;

    // add the unioned region
    if (DEBUG_WEIGHTS_VERBOSE) Rprintf(" (unioned) %d | ", unioned_region_size);

    for(int region_id = 0 ; region_id < plan.num_regions; region_id++) {
        auto region_size = plan.region_sizes[region_id];

        // add if valid multidistrict and not the two we started with
        if(valid_region_sizes_to_split[region_size] &&
            region_id != region1_id && 
            region_id != region2_id){
            if (DEBUG_WEIGHTS_VERBOSE) Rprintf(" %d | ", region_size);
            // add the count and label to vector
            prob_sum += std::pow(region_size, selection_alpha);
        }
    }

    // so prob of picking is weight of first over sum
    double log_prob = std::log(unioned_region_prob) - std::log(prob_sum);

    if (DEBUG_WEIGHTS_VERBOSE){
    Rprintf("\nFor regions (%d,%d), size (%u,%u) - Prob is %.5f/%.5f, so log is %.5f\n",
    region1_id, region2_id, 
    plan.region_sizes[region1_id], plan.region_sizes[region2_id],
    unioned_region_prob, prob_sum,
    log_prob);
    }

    return log_prob;

}


// computes the backwards kernel that is uniform in 
// the number of ancestors 
double compute_simple_log_incremental_weight(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, 
    USTSampler &ust_sampler, TreeSplitter const &edge_splitter,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double rho,
    bool compute_log_splitting_prob, bool is_final_plan
){
    // bool for whether we'll need to compute spanning tree count
    bool compute_log_tau = rho != 1;
    // find the index of the two newly split regions 
    auto most_recent_split_regions = plan.get_most_recently_split_regions();
    auto region1_id = most_recent_split_regions.first;
    auto region2_id = most_recent_split_regions.second;

    if(DEBUG_WEIGHTS_VERBOSE){
        Rprintf("The two regions are %d and %d\n", region1_id, region2_id);
        Rcpp::Rcerr << std::flush;
    } 

    bool const use_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;

    // build the plan multigraph and get valid adj pairs 
    auto valid_adj_pairs = plan.get_valid_smc_merge_regions(
        plan_multigraph, splitting_schedule, scoring_function
    );


    // compute forward and backwards kernel term 
    double log_backwards_kernel_term = 0.0; double log_forward_kernel_term = 0.0;
    double log_tau_ratio_term = 0.0;

    if(sampling_space == SamplingSpace::GraphSpace){
        // the number of valid adjacent regions is just the number of elements in the map
        const int num_valid_adj_region_pairs = valid_adj_pairs.size();
        log_backwards_kernel_term -= std::log(
            static_cast<double>(num_valid_adj_region_pairs)
        );
        // taus cancel so just add the boundary length which we have from hash map 
        log_forward_kernel_term = plan.get_log_eff_boundary_len(
            plan_multigraph, splitting_schedule, 
            ust_sampler, edge_splitter,
            region1_id, region2_id
        );
        if(compute_log_tau){
            // do merged region tau if neccesary
            log_tau_ratio_term -= (rho-1)*plan.compute_log_merged_region_spanning_trees(
                plan_multigraph.map_params, region1_id, region2_id
            );
        }
    }else if(sampling_space == SamplingSpace::ForestSpace || sampling_space == SamplingSpace::LinkingEdgeSpace){
        // sum of spanning trees 
        double spanning_tree_sum = 0.0;
        double last_split_merge_log_tau = 0.0;

        for (auto const &region_pair : valid_adj_pairs){
            int pair_region1_id = region_pair.first;
            int pair_region2_id = region_pair.second;
            
            // double log_st_term = plan.compute_log_merged_region_spanning_trees(
            //     map_params, region1_id, region2_id
            // );
            // double st_term = std::exp(log_st_term);
            // Rprintf("Pair (%d, %d) = %f and %f \n", region1_id, region2_id,
            //     log_st_term, st_term);

            // If its the merged region actually split save it for use in forward kernel
            if(pair_region1_id == region1_id && pair_region2_id == region2_id){
                last_split_merge_log_tau = plan.compute_log_merged_region_spanning_trees(
                    plan_multigraph.map_params, region1_id, region2_id
                );
                spanning_tree_sum += std::exp(last_split_merge_log_tau);
            }else{
                spanning_tree_sum += std::exp(plan.compute_log_merged_region_spanning_trees(
                    plan_multigraph.map_params, pair_region1_id, pair_region2_id
                ));
            }
        }

        // Rprintf("The total is %f!\n", spanning_tree_sum);
        log_backwards_kernel_term -= std::log(spanning_tree_sum);

        // One over tau of merged region 
        log_forward_kernel_term -= last_split_merge_log_tau;
        // tree boundary length 
        log_forward_kernel_term += plan.get_log_eff_boundary_len(
            plan_multigraph, splitting_schedule, 
            ust_sampler, edge_splitter,
            region1_id, region2_id
        );

        if(compute_log_tau){
            log_tau_ratio_term -= (rho-1)*last_split_merge_log_tau;
        }
    }
    if(compute_log_tau){
        // if(plan.region_sizes[region1_id] == 1){
        //     log_tau_ratio_term += (rho-1)*plan.compute_log_region_spanning_trees(
        //         plan_multigraph.map_params, region1_id
        //     );
        // }
        // if(plan.region_sizes[region2_id] == 1){
        //     log_tau_ratio_term += (rho-1)*plan.compute_log_region_spanning_trees(
        //         plan_multigraph.map_params, region2_id
        //     );
        // }
        log_tau_ratio_term += (rho-1)*plan.compute_log_region_spanning_trees(
            plan_multigraph.map_params, region1_id
        );
        log_tau_ratio_term += (rho-1)*plan.compute_log_region_spanning_trees(
            plan_multigraph.map_params, region2_id
        );
    }


    double log_splitting_prob = 0;
    // for one district split the probability that region was chosen to be split is always 1
    if(compute_log_splitting_prob){
        // in generalized region split find probability you would have 
        // picked to split the union of the the two regions 
        log_splitting_prob = get_log_retroactive_splitting_prob(plan, 
        splitting_schedule.valid_region_sizes_to_split, region1_id, region2_id);
        if(DEBUG_WEIGHTS_VERBOSE) Rprintf("Computed split prob %f\n", std::exp(log_splitting_prob));
    }

    double region1_score, region2_score, merged_region_score;

    region1_score = region2_score = merged_region_score = 0;

    // compute if any constraints 
    if(scoring_function.any_soft_region_constraints){
        // compute scoring functions
        region1_score = scoring_function.compute_region_score(
            plan, region1_id, is_final_plan
        );
        region2_score = scoring_function.compute_region_score(
            plan, region2_id, is_final_plan
        );
        merged_region_score = scoring_function.compute_merged_region_score(
            plan, region1_id, region2_id, is_final_plan
        );
        if(DEBUG_WEIGHTS_VERBOSE){
            REprintf("Region 1,2 Scores (%f, %f) | Merged Score %f \n",
            region1_score, region2_score, merged_region_score);
        }
    }

    double plan_score, prev_plan_score;
    plan_score = prev_plan_score = 0.0;

    if(scoring_function.any_soft_plan_constraints){
        plan_score = scoring_function.compute_plan_score(plan, is_final_plan);
        prev_plan_score = scoring_function.compute_merged_plan_score(plan, region1_id, region2_id, is_final_plan);
        if(DEBUG_WEIGHTS_VERBOSE){
            REprintf("Entire Plan Score %f | Previous Plan Score %f \n",
            plan_score, prev_plan_score);
        }
    }

    double log_extra_plan_terms = 0.0;
    double log_extra_prev_plan_terms = 0.0;
    // if linking edge space we also need to correct for that
    if(sampling_space == SamplingSpace::LinkingEdgeSpace){
        std::vector<int> merge_index_reshuffle(plan.num_regions);
        // we divide target by number of linking edges so 
        // subtract log linking edges from denominator 
        log_extra_plan_terms -= plan_multigraph.compute_log_multigraph_tau(plan.num_regions, merge_index_reshuffle, scoring_function);

        

            if(plan.num_regions > 2){ 


            // TEMP just rebuild the multigraph 
            std::vector<RegionID> flattened_all_plans(plan_multigraph.map_params.V);
            PlanVector plan_region_ids(flattened_all_plans, 0, plan_multigraph.map_params.V);
            // REprintf("Size is %u!\n", plan_region_ids.size());

            // set merge reindex 
            int const merged_reindex = plan.num_regions-2;
            for (int current_reindex = 0, i = 0; i < plan.num_regions; i++){
                if(i == region1_id || i == region2_id){
                    merge_index_reshuffle[i] = merged_reindex;
                }else{
                    merge_index_reshuffle[i] = current_reindex;
                    ++current_reindex;
                }
                // REprintf("Mapping %d to %d!\n", i, merge_index_reshuffle[i]);
            }

            // REprintf("Merging (%u, %u)\n", region1_id, region2_id);
            for (size_t i = 0; i < plan.region_ids.size(); i++)
            {
                plan_region_ids[i] = merge_index_reshuffle[plan.region_ids[i]];
                // REprintf("Plan %u | Merged %u \n", plan.region_ids[i], plan_region_ids[i]);
            }

            PlanMultigraph temp_multi(plan_multigraph.map_params);
            temp_multi.build_plan_multigraph(plan_region_ids, plan.num_regions -1);

            // double merged_tau = plan_multigraph.compute_merged_log_multigraph_tau(
            //     plan.num_regions, merge_index_reshuffle,
            //     region1_id, region2_id, scoring_function
            // );
            double temp_tau = temp_multi.compute_log_multigraph_tau(plan.num_regions-1, merge_index_reshuffle, scoring_function);

            log_extra_prev_plan_terms -= temp_tau;

            // if(merged_tau != temp_tau){
            //     REprintf("Merged Plan - Old Code - %f, New Code - %f and Equality Check = %s\n",
            //     merged_tau, temp_tau, 
            //     (merged_tau == temp_tau) ? "TRUE" : "FALSE" );
            // }

            }

        // log_extra_prev_plan_terms -= plan_multigraph.compute_merged_log_multigraph_tau(
        //         plan.num_regions, merge_index_reshuffle,
        //         region1_id, region2_id, scoring_function
        //     );

            
    }

    // The weight is 
    //      - Numerator: backwards kernel * e^-(J(new_region1)+J(new_region2)) * e^-J(new plan)
    //      - Denominator: multid_selection_prob * forward kernel term * e^-J(old_region) * e^-J(old plan)
    // So
    //      - log numerator: backwards kernel -J(new_region1)+J(new_region2))
    //      - log Denominator: log(multid_selection_prob) + log(forward kernel) - J(old_region)

    const double log_numerator = log_backwards_kernel_term - (region1_score + region2_score + plan_score) + log_extra_plan_terms;
    const double log_denom =  log_forward_kernel_term + log_splitting_prob - (merged_region_score + prev_plan_score) + log_extra_prev_plan_terms;


    if(DEBUG_WEIGHTS_VERBOSE){
    int region1_size = plan.region_sizes[region1_id];
    int region2_size = plan.region_sizes[region2_id];
    Rprintf("Doing (%d,%d) - sizes (%d, %d): forward %f, backward %f, split prob %f, ratio %f!\n", 
        region1_id, region2_id, region1_size, region2_size, std::exp(log_forward_kernel_term), 
        std::exp(log_backwards_kernel_term), 
        std::exp(log_splitting_prob), std::exp(log_numerator - merged_region_score));

    Rprintf("Numerator - %f * %f * %f\n", 
        std::exp(region1_score), std::exp(region2_score), 
        std::exp(log_backwards_kernel_term)
    );

    Rprintf("Denominator - %f * %f * %f \n", 
        std::exp(log_splitting_prob),
        std::exp(log_forward_kernel_term),
        std::exp(merged_region_score)
    );}

    // now just take the fraction
    double incremental_weight = log_numerator - log_denom + log_tau_ratio_term;
    // DONT KNOW WHY BUT VALIDATION ONLY WORKS IF YOU SUBTRACT LOG TAU
    // THEORY SAYS IT SHOULD BE ADD SO NEED TO INVESTIGATE 

    // Check its not infinity
    if (!std::isfinite(incremental_weight)) {
        plan.Rprint(true);
        throw Rcpp::exception(
            "One of the plan incremental weights is not finite!"
            "Try checking if constraint strength is too large and causing overflow errors.\n"
        );
    }


    return incremental_weight;
}


void compute_all_plans_log_simple_incremental_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    TreeSplitter const &tree_splitter,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    int verbosity
){
    int const M = (int) plans_ptr_vec.size();
    int const num_regions = plans_ptr_vec[0]->num_regions;
    const int check_int = 50; // check for interrupts every _ iterations


    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        static thread_local PlanMultigraph plan_multigraph(map_params);
        static thread_local USTSampler ust_sampler(map_params, splitting_schedule);

        double log_incr_weight = compute_simple_log_incremental_weight(
            *plans_ptr_vec[i], plan_multigraph,
            splitting_schedule, ust_sampler, tree_splitter, 
            sampling_space,
            scoring_function, rho,
            compute_log_splitting_prob, 
            is_final_plans
        );

        log_incremental_weights[i] = log_incr_weight;

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();


    return;
}



// eventually need to modify to allow presaved options
// OLD DOCUMENTATION FROM GRAPH THING NEED TO UPDATE
//' Compute the optimal log incremental weight of a plan
//'
//' Given a plan object this computes the minimum variance weights as derived in
//' <PAPER NAME HERE>. This is equal to the inverse of a sum over all
//' adjacent regions in a plan if using generalized region split and 
//' sum over all districts adajacent to the remainder if using one district
//' split approach.
//'
//' @title Compute Optimal Incremental Weight of a plan
//'
//' @param g The underlying map graph
//' @param plan A plan object
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not.
//' @param target The target population for a single district
//' @param pop_temper The population tempering parameter
//'
//' @details No modifications to inputs made
//'
//' @return the log of the incremental weight of the plan
//'
double compute_log_optimal_incremental_weights(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, 
    USTSampler &ust_sampler, TreeSplitter const &edge_splitter,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double const rho,
    bool compute_log_splitting_prob, bool is_final_plan
){
    // plan.Rprint();
    // bool for whether we'll need to compute spanning tree count
    bool const compute_log_tau = rho != 1;

    // boolean for whether or not to compute the splitting probability of merged regions
    // dont need to do when
    double incremental_weight = 0.0;

    if(DEBUG_WEIGHTS_VERBOSE) Rprintf("Getting Pairs!\n");


    // get region pair to effective boundary length map
    auto region_pair_log_eff_boundary_map = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
        plan_multigraph, splitting_schedule, scoring_function,
        ust_sampler, edge_splitter
    );

    // get region multigraph if using linked edge 
    bool const use_linked_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    std::vector<int> merge_index_reshuffle(
        use_linked_edge_space ? plan.num_regions : 0
    );

    // iterate over the pairs 
    if(DEBUG_WEIGHTS_VERBOSE){
        REprintf("There are %d adjacent pairs!\n",
            region_pair_log_eff_boundary_map.size()
        );
    }

    // compute plan score if needed 
    double plan_score = 0.0;

    if(scoring_function.any_soft_plan_constraints){
        plan_score += scoring_function.compute_plan_score(plan, is_final_plan);
        if(DEBUG_WEIGHTS_VERBOSE){
            REprintf("Entire Plan Score %f \n",
            plan_score);
        }
    }
    

    // Now iterate over adjacent region pairs and add splitting and pop temper
    for (const auto& pair_tuple: region_pair_log_eff_boundary_map){
        double log_of_sum_term = 0.0;

        const int region1_id = std::get<0>(pair_tuple); // get the smaller region id  
        const int region2_id = std::get<1>(pair_tuple); // get the bigger region id  
        const double eff_log_boundary_len = std::get<2>(pair_tuple); // get the effective boundary length

        if(DEBUG_WEIGHTS_VERBOSE){
            REprintf("Pair (%d,%d) - %f\n", region1_id, region2_id, std::exp(eff_log_boundary_len));
        }
        
        // add to term
        log_of_sum_term += eff_log_boundary_len;

        // If not compute_log_splitting_prob then its just log(1) = 0
        if(compute_log_splitting_prob){
            // in generalized region split find probability you would have 
            // picked to split the union of the the two regions 
            double log_splitting_prob = get_log_retroactive_splitting_prob(
                plan, splitting_schedule.valid_region_sizes_to_split, 
                region1_id, region2_id
            );
            log_of_sum_term += log_splitting_prob;
            if(DEBUG_WEIGHTS_VERBOSE) REprintf("\n");
        }


        // If using linked edge add multigraph tau
        if(use_linked_edge_space){
            // log_of_sum_term -= plan_multigraph.compute_merged_log_multigraph_tau(
            //     plan.num_regions, merge_index_reshuffle,
            //     region1_id, region2_id, scoring_function
            // );

            if(plan.num_regions > 2){ 


            // TEMP just rebuild the multigraph 
            std::vector<RegionID> flattened_all_plans(plan_multigraph.map_params.V);
            PlanVector plan_region_ids(flattened_all_plans, 0, plan_multigraph.map_params.V);
            // REprintf("Size is %u!\n", plan_region_ids.size());

            // set merge reindex 
            int const merged_reindex = plan.num_regions-2;
            for (int current_reindex = 0, i = 0; i < plan.num_regions; i++){
                if(i == region1_id || i == region2_id){
                    merge_index_reshuffle[i] = merged_reindex;
                }else{
                    merge_index_reshuffle[i] = current_reindex;
                    ++current_reindex;
                }
                // REprintf("Mapping %d to %d!\n", i, merge_index_reshuffle[i]);
            }

            // REprintf("Merging (%u, %u)\n", region1_id, region2_id);
            for (size_t i = 0; i < plan.region_ids.size(); i++)
            {
                plan_region_ids[i] = merge_index_reshuffle[plan.region_ids[i]];
                // REprintf("Plan %u | Merged %u \n", plan.region_ids[i], plan_region_ids[i]);
            }

            PlanMultigraph temp_multi(plan_multigraph.map_params);
            temp_multi.build_plan_multigraph(plan_region_ids, plan.num_regions -1);

            // double merged_tau = plan_multigraph.compute_merged_log_multigraph_tau(
            //     plan.num_regions, merge_index_reshuffle,
            //     region1_id, region2_id, scoring_function
            // );
            double temp_tau = temp_multi.compute_log_multigraph_tau(plan.num_regions-1, merge_index_reshuffle, scoring_function);

            log_of_sum_term -= temp_tau;

            // if(merged_tau != temp_tau){
            //     REprintf("Merged Plan - Old Code - %f, New Code - %f and Equality Check = %s\n",
            //     merged_tau, temp_tau, 
            //     (merged_tau == temp_tau) ? "TRUE" : "FALSE" );
            // }

            }
            
        }


        // compute score ratio if any constraints 
        if(scoring_function.any_soft_region_constraints){
            // compute scoring functions
            const double region1_score = scoring_function.compute_region_score(
                plan, region1_id, is_final_plan
            );
            const double region2_score = scoring_function.compute_region_score(
                plan, region2_id, is_final_plan
            );
            const double merged_region_score = scoring_function.compute_merged_region_score(
                plan, region1_id, region2_id, is_final_plan
            );
            if(DEBUG_WEIGHTS_VERBOSE){
                REprintf("Region 1,2 Scores (%f, %f) | Merged Score %f \n",
                region1_score, region2_score, merged_region_score);
            }
            // log ratio is (log region1 + log region2) - log score merged
            double const score_ratio = region1_score + region2_score - merged_region_score;
            log_of_sum_term += score_ratio;
        }

        
        if(scoring_function.any_soft_plan_constraints){
            double merged_plan_score = scoring_function.compute_merged_plan_score(plan, region1_id, region2_id, is_final_plan);
            if(DEBUG_WEIGHTS_VERBOSE){
                REprintf(
                    "For Regions (%u, %u) Merged Entire Plan Score %f \n",
                    region1_id, region2_id, merged_plan_score
                );
            }
            // term is e^- score(merged plan) so we subtract the log score 
            log_of_sum_term -= merged_plan_score;
        }

        if(DEBUG_WEIGHTS_VERBOSE){
        int region1_size = plan.region_sizes[region1_id];
        int region2_size = plan.region_sizes[region2_id];
        Rprintf("Adding (%d,%d) - sizes (%d, %d): len %f, split prob %f, ratio %f!\n", 
            region1_id, region2_id, region1_size, region2_size, std::exp(eff_log_boundary_len), 
            std::exp(log_of_sum_term), log_of_sum_term);
        }

        

        // Do taus if neccesary 
        if(compute_log_tau){
            double log_tau_ratio = 0.0;
            // add merged region
            log_tau_ratio += (rho-1)*plan.compute_log_merged_region_spanning_trees(
                plan_multigraph.map_params, region1_id, region2_id
            );
            // subtract split regions
            log_tau_ratio -= (rho-1)*plan.compute_log_region_spanning_trees(
                plan_multigraph.map_params, region1_id
            );
            log_tau_ratio -= (rho-1)*plan.compute_log_region_spanning_trees(
                plan_multigraph.map_params, region2_id
            );
            log_of_sum_term += log_tau_ratio;
        }

        // Now exponentiate and add to the sum
        incremental_weight += std::exp(log_of_sum_term);

    }

    // Check its not infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        Rcpp::stop("Error! weight is negative infinity for some reason \n");
    }

    // Extra term if needed
    double extra_log_terms = 0.0;
    if(use_linked_edge_space){
        // need number of linking edges for current plan
        extra_log_terms -= plan_multigraph.compute_log_multigraph_tau(
            plan.num_regions, merge_index_reshuffle, scoring_function
        );
    }

    // Now add the log extra terms and subtract the plan score 
    incremental_weight = extra_log_terms - std::log(incremental_weight) - plan_score;

    if(DEBUG_WEIGHTS_VERBOSE){
        REprintf("Weight=%f, log weight = %f\n", std::exp(incremental_weight), incremental_weight);
    }

    
    if (!std::isfinite(incremental_weight)) {
        plan.Rprint(true);
        throw Rcpp::exception(
            "One of the plan incremental weights is not finite!"
            "Try checking if constraint strength is too large and causing overflow errors.\n"
        );
    }

    // now return the log of the inverse of the sum
    return incremental_weight;

}




//' NEED TO UPDATE THIS IS OLD DOCUMENTATION FOR GRAPH STUFF
//' Computes log unnormalized weights for vector of plans
//'
//' Using the procedure outlined in <PAPER HERE> this function computes the log
//' incremental weights and the unnormalized weights for a vector of plans (which
//' may or may not be the same depending on the parameters).
//'
//' @title Compute Log Unnormalized Weights
//'
//' @param pool A threadpool for multithreading
//' @param g A graph (adjacency list) passed by reference
//' @param plans_ptr_vec A vector of plans to compute the log unnormalized weights
//' of
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not. If `split_district_only` is true
//' then uses optimal weights from one-district split scheme.
//' @param log_incremental_weights A vector of the log incremental weights
//' computed for the plans. The value of `log_incremental_weights[i]` is
//' the log incremental weight for `plans_ptr_vec[i]`
//' @param unnormalized_sampling_weights A vector of the unnormalized sampling
//' weights to be used with sampling the `plans_ptr_vec` in the next iteration of the
//' algorithm. Depending on the other hyperparameters this may or may not be the
//' same as `exp(log_incremental_weights)`
//' @param target Target population of a single district
//' @param pop_temper <DETAILS NEEDED>
//'
//' @details Modifications
//'    - The `log_incremental_weights` is updated to contain the incremental
//'    weights of the plans
//'    - The `unnormalized_sampling_weights` is updated to contain the unnormalized
//'    sampling weights of the plans for the next round
void compute_all_plans_log_optimal_incremental_weights(
    RcppThread::ThreadPool &pool,
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function,
    double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    TreeSplitter const &tree_splitter,
    bool compute_log_splitting_prob, bool is_final_plans,
    arma::subview_col<double> log_incremental_weights,
    int verbosity
){
    const int nsims = static_cast<int>(plans_ptr_vec.size());
    const int check_int = 50; // check for interrupts every _ iterations
    const int num_regions = plans_ptr_vec[0]->num_regions;
    if(DEBUG_WEIGHTS_VERBOSE) Rprintf("About to start computing weights!\n");

    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        static thread_local PlanMultigraph plan_multigraph(map_params);
        static thread_local USTSampler ust_sampler(map_params, splitting_schedule);


        // REprintf("I=%d\n", i);
        double log_incr_weight = compute_log_optimal_incremental_weights(
            *plans_ptr_vec[i], plan_multigraph, 
            splitting_schedule, ust_sampler, tree_splitter,
            sampling_space, scoring_function, 
            rho, compute_log_splitting_prob, 
            is_final_plans
        );

        log_incremental_weights[i] = log_incr_weight;

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}