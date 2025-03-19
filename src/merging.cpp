/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Merging Plans
********************************************************/

#include "merging.h"

constexpr bool DEBUG_MERGING_VERBOSE = false; // Compile-time constant

// computes log metropolis hastings ratio
double get_log_mh_ratio(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    const int region1_id, const int region2_id,
    double const current_log_eff_boundary_len, 
    double const proposed_log_eff_boundary_len, 
    const double log_current_pair_merge_prob, 
    const double log_new_pair_merge_prob, 
    Plan const &current_plan, Plan const &proposed_plan,
    double const rho, bool const is_final
){
    double log_mh_ratio = 0.0;
    // We add the scores of the current regions 
    log_mh_ratio += scoring_function.compute_region_score(
        current_plan, region1_id, is_final
    );
    log_mh_ratio += scoring_function.compute_region_score(
        current_plan, region2_id, is_final
    );
    // We subtract the scores of the proposed region
    log_mh_ratio -= scoring_function.compute_region_score(
        proposed_plan, region1_id, is_final
    );
    log_mh_ratio -= scoring_function.compute_region_score(
        proposed_plan, region2_id, is_final
    );
    // we compute taus if neccesary 
    if(rho != 1){
        // we add scorse of proposed region
        log_mh_ratio += (rho-1) * proposed_plan.compute_log_region_spanning_trees(
            map_params, region1_id
        );
        log_mh_ratio += (rho-1) * proposed_plan.compute_log_region_spanning_trees(
            map_params, region2_id
        );
        // we subtract scores of current region
        log_mh_ratio -= (rho-1) * current_plan.compute_log_region_spanning_trees(
            map_params, region1_id
        );
        log_mh_ratio -= (rho-1) * current_plan.compute_log_region_spanning_trees(
            map_params, region2_id
        );
    }
    // we add forward kernel proposed to current terms 
    log_mh_ratio += current_log_eff_boundary_len;
    log_mh_ratio += log_new_pair_merge_prob;
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("We added %f and %f, ", 
            current_log_eff_boundary_len, log_new_pair_merge_prob
        );
    }
    // we subtract forward kernel current to proposed terms 
    log_mh_ratio -= proposed_log_eff_boundary_len;
    log_mh_ratio -= log_current_pair_merge_prob;
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("We subtracted %f and %f, ", 
            proposed_log_eff_boundary_len, log_current_pair_merge_prob
        );
        Rprintf("Ratio is now %f\n", std::exp(log_mh_ratio));
    }

    return log_mh_ratio;
}

// runs a mergesplit Metropolis Hastings Step
std::tuple<bool, bool, double> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state,
    Plan &plan, Plan &new_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    std::string const merge_prob_type, bool save_edge_selection_prob,
    std::vector<std::pair<int,int>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts,
    double const rho, bool const is_final
){
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("Doing regions %d, %d!\n", adj_region_pairs[0].first, adj_region_pairs[0].second);
    }
    // sample a pair 
    int sampled_pair_index = rng_state.r_int_unnormalized_wgt(unnormalized_pair_wgts);
    std::pair<int, int> merge_pair = adj_region_pairs[sampled_pair_index];

    int region1_id = merge_pair.first; 
    int region2_id = merge_pair.second;

    if(DEBUG_MERGING_VERBOSE){
        Rprintf("Picked pair (%d, %d)\n", region1_id, region2_id);
    }

    // try to draw a region 
    std::tuple<bool, EdgeCut, double> edge_search_result = ust_sampler.attempt_to_find_valid_tree_mergesplit(
        map_params, splitting_schedule,
        rng_state, tree_splitter,
        plan, region1_id, region2_id,
        save_edge_selection_prob
    );
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("A Splitting Checkpoint 1.\n");
    }
    
    // If nothing drawn immediately return
    if(!std::get<0>(edge_search_result)){
        if(DEBUG_MERGING_VERBOSE){
            Rprintf("Failed!\n");
        }
        return std::make_tuple(false, false, -1*std::log(0.0));
    }


    // IN THE FUTURE CAN AVOID THE COPYING BY JUST TRAVERSING THE TREE
    // Just traverse tree and check if not in merged region or something

    // shallow copy the new plan to be the old one 
    new_plan = plan;
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("A Splitting Checkpoint 1.5!\n");
    }
    // now split that region we found on the old one
    new_plan.update_from_successful_split(
        ust_sampler.ust, std::get<1>(edge_search_result),
        region1_id, region2_id,
        std::get<2>(edge_search_result), 
        false
    );
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("A Splitting Checkpoint 2.\n");
    }
    // get adj pairs 
    auto new_valid_adj_region_pairs = new_plan.get_valid_adj_regions(
        map_params, splitting_schedule
    );
    // get the weights
    auto new_valid_pair_weights = get_adj_pair_unnormalized_weights(
        new_plan,
        new_valid_adj_region_pairs,
        merge_prob_type
    );
    int region_pair_proposal_index = -1;
    // Find the index of the pair in the proposed plan
    for (size_t i = 0; i < new_valid_adj_region_pairs.size(); i++)
    {

        if(new_valid_adj_region_pairs[i].first == region1_id 
        && new_valid_adj_region_pairs[i].second == region2_id){
            region_pair_proposal_index = i;
            break;
        }
    }
    // means disconnected thing glitch. If encountered just print and ignore
    if(region_pair_proposal_index == -1){
        REprintf("Merged and Split Regions %d and %d no longer adjacent!\n",
            region1_id, region2_id);
        throw Rcpp::exception("Error in merge split!\n");
    }
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("selected new pair index is %d!\n", region_pair_proposal_index);
    }
    // compute the boundary length 
    double current_log_eff_boundary = plan.get_log_eff_boundary_len(
        map_params, splitting_schedule, tree_splitter,
        region1_id, region2_id
    );
    double proposed_log_eff_boundary = new_plan.get_log_eff_boundary_len(
        map_params, splitting_schedule, tree_splitter,
        region1_id, region2_id
    );

    if(DEBUG_MERGING_VERBOSE){
        Rprintf("Doing regions %d, %d!\n", region1_id, region2_id);
        Rprintf(
            "Current Plan: %d Adjacent Regions and I picked index %d ",
            (int) adj_region_pairs.size(), sampled_pair_index
        );
        Rprintf(
            "and %f\n",
            std::exp(std::log(unnormalized_pair_wgts(sampled_pair_index)) - std::log(arma::sum(unnormalized_pair_wgts)))
        );
        Rprintf(
            "Proposed Plan: %d Adjacent Regions and I picked index %d ",
            (int) new_valid_adj_region_pairs.size(), region_pair_proposal_index
        );
        Rprintf(
            "and %f\n",
            std::exp(std::log(new_valid_pair_weights(region_pair_proposal_index)) - std::log(arma::sum(new_valid_pair_weights)))
        );
    }

    double log_mh_ratio = get_log_mh_ratio(
        map_params, scoring_function,
        region1_id, region2_id,
        current_log_eff_boundary, 
        proposed_log_eff_boundary, 
        std::log(unnormalized_pair_wgts(sampled_pair_index)) - std::log(arma::sum(unnormalized_pair_wgts)), 
        std::log(new_valid_pair_weights(region_pair_proposal_index)) - std::log(arma::sum(new_valid_pair_weights)), 
        plan, new_plan,
        rho, is_final
    );
    bool proposal_accepted = rng_state.r_unif() <= std::exp(log_mh_ratio);
    if(proposal_accepted){
        // if successful then actually update
        plan.update_from_successful_split(
            ust_sampler.ust, std::get<1>(edge_search_result),
            region1_id, region2_id,
            std::get<2>(edge_search_result), 
            false
        );
        // update pairs and weights to be current one
        unnormalized_pair_wgts = new_valid_pair_weights;
        adj_region_pairs = new_valid_adj_region_pairs;
        return std::make_tuple(true, true, log_mh_ratio);
    }else{
        // else reject and do nothing 
        return std::make_tuple(true, false, log_mh_ratio);
    }
}


// runs num_steps_to_run merge split stesp
int run_merge_split_steps(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &dummy_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    std::string const merge_prob_type,
    double const rho, bool const is_final, 
    int num_steps_to_run
){
    int num_succesful_steps = 0;
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // Get pairs of adj districts
    std::vector<std::pair<int,int>> current_plan_adj_region_pairs = plan.get_valid_adj_regions(
        map_params, splitting_schedule
    );
    arma::vec current_plan_pair_unnoramalized_wgts = get_adj_pair_unnormalized_weights(
        plan,
        current_plan_adj_region_pairs,
        merge_prob_type
    );
    // run the merge split steps and count success 
    for (size_t i = 0; i < num_steps_to_run; i++)
    {
        std::tuple<bool, bool, double> mergesplit_result = attempt_mergesplit_step(
            map_params, splitting_schedule, scoring_function,
            rng_state,
            plan, dummy_plan, 
            ust_sampler, tree_splitter,
            merge_prob_type, save_edge_selection_prob,
            current_plan_adj_region_pairs,
            current_plan_pair_unnoramalized_wgts,
            rho, is_final
        );
        // increase count if successful
        if(std::get<1>(mergesplit_result)){
            ++num_succesful_steps;
        }
    }
    return num_succesful_steps;
    
}



void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec, 
    TreeSplitter const &tree_splitter,
    std::string const merge_prob_type, 
    double const rho, bool const is_final, 
    int const nsteps_to_run,
    Rcpp::IntegerMatrix::Column success_count_vec
){
    int M = (int) plan_ptrs_vec.size();

    // thread safe id counter for seeding RNG generator 
    std::atomic<int> thread_id_counter{0};
    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        static thread_local int thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
        static thread_local USTSampler ust_sampler(map_params.V);
        // Create variables needed for each 

        // store the number of succesful runs
        success_count_vec[i] = run_merge_split_steps(
            map_params, splitting_schedule, scoring_function,
            rng_states[thread_id], sampling_space,
            *plan_ptrs_vec.at(i), *new_plan_ptrs_vec.at(i), 
            ust_sampler, tree_splitter,
            merge_prob_type, 
            rho, is_final, 
            nsteps_to_run
        );

    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}

