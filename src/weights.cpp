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






// Returns a sampler over a vector of adjacent pairs where the probability 
// of a pair is decided according to `selection_type`
//'
//' Current supported options are
//'     - uniform - Every pair has equal probability
//'     - district_pair - double district pairs have weight 1000, one district is 10,
//'         and two multidistricts have 1/(1+sum of their dvals)
//'
//' @title Get Sampler over Adj Regions List
//'
//' @param plan A plan object
//' @param adj_pairs_and_boundary_lens A vector where each pair is 
//' (adj region1, adj region2, boundary length between 2 regions)
//' @param selection_type A string controlling the function to use
//' in assigning the unnormalized weight to each pair
//'
//' @details No modifications to inputs made
//'
//' @return A sampler where index i has probability proportional to the weight 
//' given to that pair 
//'
arma::vec get_adj_pair_sampler(
    Plan const &plan,
    std::vector<std::tuple<int, int, double>> const &adj_pairs_and_boundary_lens,
    std::string const &selection_type
){
    // make a vector for the unnormalized weights 
    arma::vec unnormalized_sampling_weights(adj_pairs_and_boundary_lens.size());

    // if uniform then just make all the weights 1 over the size 
    if(selection_type == "uniform"){
        std::fill(
            unnormalized_sampling_weights.begin(), 
            unnormalized_sampling_weights.end(), 
            static_cast<double>(adj_pairs_and_boundary_lens.size())
            );
    }else if(selection_type == "district_pair"){
        // if district pair then give huge weight to pairs of districts 
        for (size_t i = 0; i < adj_pairs_and_boundary_lens.size(); i++)
        {
            int region1_id = std::get<0>(adj_pairs_and_boundary_lens[i]);
            int region2_id = std::get<1>(adj_pairs_and_boundary_lens[i]);
            int region1_dval = plan.region_sizes[region1_id]; int region2_dval = plan.region_sizes[region2_id];

            // check if both are districts
            if(region1_dval == 1 && region2_dval == 1){
                unnormalized_sampling_weights.at(i) = 1000.0;
            }else if(region1_dval == 1 || region2_dval == 1){
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 10.0;
            }else{
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs 
                unnormalized_sampling_weights.at(i) = 1/(1+ static_cast<double>(region1_dval+region2_dval));
            }
        }
        
    }else if(selection_type == "multidistrict_pair"){
        // if multidistrict pair then give huge weight to pairs of districts 
        for (size_t i = 0; i < adj_pairs_and_boundary_lens.size(); i++)
        {
            int region1_id = std::get<0>(adj_pairs_and_boundary_lens[i]);
            int region2_id = std::get<1>(adj_pairs_and_boundary_lens[i]);
            int region1_dval = plan.region_sizes[region1_id]; int region2_dval = plan.region_sizes[region2_id];

            // check if both are districts
            if(region1_dval == 1 && region2_dval == 1){
                unnormalized_sampling_weights.at(i) = .5;
            }else if(region1_dval == 1 || region2_dval == 1){
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 1.0;
            }else{
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs 
                unnormalized_sampling_weights.at(i) = 1000.0;
            }
        }
        
    }else{
        throw Rcpp::exception("No valid adjacent pair sampler provided!");
    }
    // now return normalized weights 
    arma::vec cum_wgts = arma::cumsum(unnormalized_sampling_weights); 
    // now normalize them
    cum_wgts = cum_wgts / cum_wgts(cum_wgts.size()-1);

    // Return 1 over the log of the size 
    return cum_wgts;
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
    if (DEBUG_WEIGHTS_VERBOSE) Rprintf(" (unioned) %d\n", unioned_region_size);

    for(int region_id = 0 ; region_id < plan.num_regions; region_id++) {
        auto region_size = plan.region_sizes[region_id];

        // add if valid multidistrict and not the two we started with
        if(valid_region_sizes_to_split[region_size] &&
            region_size > 1 && region_id != region1_id && 
            region_id != region2_id){
            if (DEBUG_WEIGHTS_VERBOSE) Rprintf(" %d ", region_size);
            // add the count and label to vector
            prob_sum += std::pow(region_size, selection_alpha);
        }
    }

    // so prob of picking is weight of first over sum
    double log_prob = std::log(unioned_region_prob) - std::log(prob_sum);

    // Rprintf("For regions (%d,%d), size (%d,%d) - Prob is %d/%d, so log is %.5f\n",
    // region1_id, region2_id, 
    // (int) plan.region_sizes[region1_id], (int) plan.region_sizes[region2_id],
    // unioned_region_dnk, total_multi_ds,
    // log_prob);

    return log_prob;

}


// computes the backwards kernel that is uniform in 
// the number of ancestors 
double compute_simple_log_incremental_weight(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, const TreeSplitter &edge_splitter,
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
    auto valid_adj_pairs = plan.get_valid_smc_merge_regions(plan_multigraph, splitting_schedule);


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
            plan_multigraph, splitting_schedule, edge_splitter,
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
            plan_multigraph, splitting_schedule, edge_splitter,
            region1_id, region2_id
        );

        if(compute_log_tau){
            log_tau_ratio_term -= (rho-1)*last_split_merge_log_tau;
        }
    }
    if(compute_log_tau){
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
    if(scoring_function.any_soft_constraints){
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
    }

    double log_extra_plan_terms = 0.0;
    double log_extra_prev_plan_terms = 0.0;
    // if linking edge space we also need to correct for that
    if(sampling_space == SamplingSpace::LinkingEdgeSpace){
        auto region_multigraph = plan_multigraph.pair_map.get_multigraph_counts(
            plan.num_regions
        );
        // we divide target by number of linking edges so 
        // subtract log linking edges from denominator 
        log_extra_plan_terms -= compute_log_region_multigraph_spanning_tree(region_multigraph);
        std::vector<int> merge_index_reshuffle(plan.num_regions);
        log_extra_prev_plan_terms -= get_log_merged_region_multigraph_spanning_tree(
                region_multigraph,
                merge_index_reshuffle,
                region1_id, region2_id
            );
    }

    // The weight is 
    //      - Numerator: backwards kernel * e^-(J(new_region1)+J(new_region2))
    //      - Denominator: multid_selection_prob * forward kernel term * e^-J(old_region)
    // So
    //      - log numerator: backwards kernel -J(new_region1)+J(new_region2))
    //      - log Denominator: log(multid_selection_prob) + log(forward kernel) - J(old_region)

    const double log_numerator = log_backwards_kernel_term - (region1_score + region2_score) + log_extra_plan_terms;
    const double log_denom =  log_forward_kernel_term + log_splitting_prob - merged_region_score + log_extra_prev_plan_terms;


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

    // Check its not - infinity
    if(incremental_weight == -std::numeric_limits<double>::infinity()){
        throw Rcpp::exception("Error! weight is negative infinity for some reason \n");
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
        PlanMultigraph plan_multigraph(map_params);

        double log_incr_weight = compute_simple_log_incremental_weight(
            *plans_ptr_vec[i], plan_multigraph,
            splitting_schedule, tree_splitter, 
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




// NOT FULLY TESTED but taken from district_graph function which was tested
//' Creates the region level graph of a plan
//'
//' Given a plan object this returns a graph of the regions in the plan using
//' the region ids as indices
//'
//' @title Get Region-Level Graph
//'
//' @param g The graph of the entire map
//' @param plan A plan object
//'
//' @details No modifications to inputs made
//'
//' @return the log of the probability the specific value of `region_to_split` was chosen
//'
Graph get_region_graph(const Graph &g, const Plan &plan) {
    int V = g.size();
    Graph out;

    // make a matrix where entry i,j represents if region i and j are adjacent
    std::vector<std::vector<bool>> gr_bool(
            plan.num_regions, std::vector<bool>(plan.num_regions, false)
        );

    // iterate over all vertices in g
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        // Find out which region this vertex corresponds to
        int region_num_i = plan.region_ids[i];

        // now iterate over its neighbors
        for (int nbor : nbors) {
            // find which region neighbor corresponds to
            int region_num_j = plan.region_ids[nbor];
            // if they are different regions mark matrix true since region i
            // and region j are adjacent as they share an edge across
            if (region_num_i != region_num_j) {
                gr_bool.at(region_num_i).at(region_num_j) = true;
            }
        }
    }


    // Now build the region level graph
    for (int i = 0; i < plan.num_regions; i++) {
        // create vector of i's neighbors
        std::vector<int> tmp;
        for (int j = 0; j < plan.num_regions; j++) {
            // check if i and j are adjacent, if so add j
            if (gr_bool.at(i).at(j)) {
                tmp.push_back(j);
            }
        }
        out.push_back(tmp);
    }

    return out;
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
    const SplittingSchedule &splitting_schedule, const TreeSplitter &edge_splitter,
    SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double const rho,
    bool compute_log_splitting_prob, bool is_final_plan
){

    // bool for whether we'll need to compute spanning tree count
    bool const compute_log_tau = rho != 1;

    // boolean for whether or not to compute the splitting probability of merged regions
    // dont need to do when
    double incremental_weight = 0.0;

    if(DEBUG_WEIGHTS_VERBOSE) Rprintf("Getting Pairs!");


    // get region pair to effective boundary length map
    auto region_pair_log_eff_boundary_map = plan.get_valid_adj_regions_and_eff_log_boundary_lens(
        plan_multigraph, splitting_schedule, edge_splitter
    );

    // get region multigraph if using linked edge 
    bool const use_linked_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    RegionMultigraphCount region_multigraph(
        use_linked_edge_space ? plan.num_regions : 0
    );
    std::vector<int> merge_index_reshuffle(
        use_linked_edge_space ? plan.num_regions : 0
    );
    if(use_linked_edge_space){
        region_multigraph = plan_multigraph.pair_map.get_multigraph_counts(plan.num_regions);
    }

    // iterate over the pairs 
    if(DEBUG_WEIGHTS_VERBOSE){
        REprintf("There are %d adjacent pairs!\n",
            region_pair_log_eff_boundary_map.size()
        );
    }
    

    // Now iterate over adjacent region pairs and add splitting and pop temper
    for (const auto& pair_tuple: region_pair_log_eff_boundary_map){
        double log_of_sum_term = 0.0;

        const int region1_id = std::get<0>(pair_tuple); // get the smaller region id  
        const int region2_id = std::get<1>(pair_tuple); // get the bigger region id  
        const double eff_log_boundary_len = std::get<2>(pair_tuple); // get the effective boundary length

        if(DEBUG_WEIGHTS_VERBOSE){
            REprintf("(%d,%d) - %f\n", region1_id, region2_id, std::exp(eff_log_boundary_len));
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
        }


        // If using linked edge add multigraph tau
        if(use_linked_edge_space){
            log_of_sum_term -= get_log_merged_region_multigraph_spanning_tree(
                region_multigraph,
                merge_index_reshuffle,
                region1_id, region2_id
            );
        }


        // compute score ratio if any constraints 
        if(scoring_function.any_soft_constraints){
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
            // log ratio is (log region1 + log region2) - log score merged
            double const score_ratio = region1_score + region2_score - merged_region_score;
            log_of_sum_term += score_ratio;
        }

        // int region1_size = plan.region_sizes[region1_id];
        // int region2_size = plan.region_sizes[region2_id];
        // Rprintf("Adding (%d,%d) - sizes (%d, %d): len %f, split prob %f, ratio %f!\n", 
        //     region1_id, region2_id, region1_size, region2_size, std::exp(eff_log_boundary_len), 
        //     std::exp(log_of_sum_term), std::exp(log_of_sum_term));

        

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
        extra_log_terms -= compute_log_region_multigraph_spanning_tree(region_multigraph);
    }

    if(DEBUG_WEIGHTS_VERBOSE){
        REprintf("Weight=%f, log weight = %f\n", incremental_weight, extra_log_terms-std::log(incremental_weight));
    }
    
    if (!std::isfinite(extra_log_terms-std::log(incremental_weight))) {
        plan.Rprint(true);
        throw Rcpp::exception("log_of_sum_term is not finite!");
    }

    // now return the log of the inverse of the sum
    return extra_log_terms-std::log(incremental_weight);

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
    if(DEBUG_WEIGHTS_VERBOSE) Rprintf("About to start computing weights!");


    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        PlanMultigraph plan_multigraph(map_params);


        // REprintf("I=%d\n", i);
        double log_incr_weight = compute_log_optimal_incremental_weights(
            *plans_ptr_vec[i], plan_multigraph, 
            splitting_schedule, tree_splitter,
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