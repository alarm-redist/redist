/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Merging Plans
********************************************************/

#include "merging.h"

constexpr bool DEBUG_MERGING_VERBOSE = false; // Compile-time constant




arma::vec get_adj_pair_unnormalized_weights(
    Plan const &plan,
    std::vector<std::pair<RegionID, RegionID>> const &valid_region_adj_pairs,
    std::string const &selection_type
){
    // make a vector for the unnormalized weights 
    arma::vec unnormalized_sampling_weights(valid_region_adj_pairs.size());

    // if uniform then just make all the weights 1 over the size 
    if(selection_type == "uniform"){
        std::fill(
            unnormalized_sampling_weights.begin(), 
            unnormalized_sampling_weights.end(), 
            static_cast<double>(valid_region_adj_pairs.size())
            );
    }else if(selection_type == "district_pair"){
        // if district pair then give huge weight to pairs of districts 
        for (size_t i = 0; i < valid_region_adj_pairs.size(); i++)
        {
            int region1_id = valid_region_adj_pairs[i].first;
            int region2_id = valid_region_adj_pairs[i].second;
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
        for (size_t i = 0; i < valid_region_adj_pairs.size(); i++)
        {
            int region1_id = valid_region_adj_pairs[i].first;
            int region2_id = valid_region_adj_pairs[i].second;
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
    // Return weights 
    return unnormalized_sampling_weights;
}

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
std::tuple<bool, bool, double, int> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &new_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph,
    std::string const merge_prob_type, bool save_edge_selection_prob,
    std::vector<std::pair<RegionID, RegionID>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts,
    double const rho, bool const is_final
){
    // sample a pair 
    int sampled_pair_index = rng_state.r_int_unnormalized_wgt(unnormalized_pair_wgts);
    std::pair<int, int> merge_pair = adj_region_pairs[sampled_pair_index];

    int region1_id = merge_pair.first; 
    int region2_id = merge_pair.second;
    int merged_region_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];

    if(DEBUG_MERGING_VERBOSE){
        Rprintf("Picked pair (%d, %d)\n", region1_id, region2_id);
    }

    // try to draw a region 
    std::tuple<bool, EdgeCut> edge_search_result = ust_sampler.attempt_to_find_valid_tree_mergesplit(
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
        return std::make_tuple(false, false, -1*std::log(0.0), merged_region_size);
    }


    // IN THE FUTURE CAN AVOID THE COPYING BY JUST TRAVERSING THE TREE
    // Just traverse tree and check if not in merged region or something

    // copy the new plan to be the old one 
    new_plan.shallow_copy(plan);
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("A Splitting Checkpoint 1.5!\n");
    }
    // now split that region we found on the old one
    new_plan.update_from_successful_split(
        tree_splitter,
        ust_sampler, std::get<1>(edge_search_result),
        region1_id, region2_id,
        false
    );
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("A Splitting Checkpoint 2.\n");
    }
    // check new plan is hierarchically valid if needed
    auto build_attempt = new_plan.attempt_to_get_valid_mergesplit_pairs(
        proposed_plan_multigraph, splitting_schedule
    );
    bool new_plan_valid = build_attempt.first;

    if(DEBUG_MERGING_VERBOSE){
        Rprintf("%d county splits!\n", proposed_plan_multigraph.num_county_region_components);
    }
    if(!new_plan_valid){
        if(DEBUG_MERGING_VERBOSE){
            REprintf("%d splits for %d regions\n", 
                proposed_plan_multigraph.num_county_region_components, plan.num_regions);
        }
        // return failure
        return std::make_tuple(true, false, -1*std::log(0.0), merged_region_size);
    }
    // get adj pairs 
    auto new_valid_adj_region_pairs = build_attempt.second;
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
        REprintf("Current adj pairs: ");
        for(const auto &a_pair: adj_region_pairs){
            REprintf("(%d, %d), ", a_pair.first, a_pair.second);
        }
        REprintf("\nProposed adj pairs: ");
        for(const auto &a_pair: new_valid_adj_region_pairs){
            REprintf("(%d, %d), ", a_pair.first, a_pair.second);
        }
        REprintf("\nCurrent Plan is:\n");
        for(const auto &id: plan.region_ids){
            REprintf("%d;", id);
        }
        REprintf("\n\nProposed Plan is:\n");
        for(const auto &id: new_plan.region_ids){
            REprintf("%d;", id);
        }
        REprintf("\n\n");
        throw Rcpp::exception("Error in merge split!\n");
    }
    if(DEBUG_MERGING_VERBOSE){
        Rprintf("selected new pair index is %d!\n", region_pair_proposal_index);
    }
    bool using_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // compute the boundary length 
    double current_log_eff_boundary = plan.get_log_eff_boundary_len(
        current_plan_multigraph, splitting_schedule, tree_splitter, 
        region1_id, region2_id
    );
    double proposed_log_eff_boundary = new_plan.get_log_eff_boundary_len(
        proposed_plan_multigraph, splitting_schedule, tree_splitter,
        region1_id, region2_id
    );
    // If linking edge space we need to subtract linking edge correction term
    if(using_linking_edge_space){
        // add instead of subtract bc its flipped in MH ratio
        current_log_eff_boundary += plan.compute_log_linking_edge_count(current_plan_multigraph);
        proposed_log_eff_boundary += new_plan.compute_log_linking_edge_count(proposed_plan_multigraph);
    }

    if(DEBUG_MERGING_VERBOSE){
        Rprintf("Finding Adjacent regions %d, %d!\n", region1_id, region2_id);
        Rprintf(
            "Current Plan: %d Adjacent Regions and I picked index %d ",
            (int) adj_region_pairs.size(), sampled_pair_index
        );
        Rprintf(
            "with probability %f\n",
            std::exp(std::log(unnormalized_pair_wgts(sampled_pair_index)) - std::log(arma::sum(unnormalized_pair_wgts)))
        );
        Rprintf(
            "Proposed Plan: %d Adjacent Regions and I picked index %d ",
            (int) new_valid_adj_region_pairs.size(), region_pair_proposal_index
        );
        Rprintf(
            "with probability %f\n",
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
    if(DEBUG_MERGING_VERBOSE) Rprintf("Ratio is %f and it is", std::exp(log_mh_ratio));
    if(proposal_accepted){
        if(DEBUG_MERGING_VERBOSE) Rprintf("ACCEPTED!! Now updating plans\n", std::exp(log_mh_ratio));
        // if successful then actually update
        plan.update_from_successful_split(
            tree_splitter,
            ust_sampler, std::get<1>(edge_search_result),
            region1_id, region2_id,
            false
        );
        // swap the plan multigraphs
        swap_plan_multigraphs(current_plan_multigraph, proposed_plan_multigraph);
        // update pairs and weights to be current one
        unnormalized_pair_wgts = new_valid_pair_weights;
        adj_region_pairs = new_valid_adj_region_pairs;
        return std::make_tuple(true, true, log_mh_ratio, merged_region_size);
    }else{
        if(DEBUG_MERGING_VERBOSE) Rprintf("REJECTED!!\n", std::exp(log_mh_ratio));
        // else reject and do nothing 
        return std::make_tuple(true, false, log_mh_ratio, merged_region_size);
    }
}


// runs num_steps_to_run merge split stesp
int run_merge_split_steps(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    RNGState &rng_state, SamplingSpace const sampling_space,
    Plan &plan, Plan &dummy_plan, 
    USTSampler &ust_sampler, TreeSplitter const &tree_splitter,
    PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph,
    std::string const merge_prob_type,
    double const rho, bool const is_final, 
    int num_steps_to_run,
    std::vector<int> &tree_sizes, std::vector<int> &successful_tree_sizes
){
    int num_succesful_steps = 0;
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;

    // Build the multigraph and get pairs of adj districts
    auto current_plan_adj_region_pairs = plan.attempt_to_get_valid_mergesplit_pairs(
        current_plan_multigraph, splitting_schedule
    ).second;
    arma::vec current_plan_pair_unnoramalized_wgts = get_adj_pair_unnormalized_weights(
        plan,
        current_plan_adj_region_pairs,
        merge_prob_type
    );
    // run the merge split steps and count success 
    for (size_t i = 0; i < num_steps_to_run; i++)
    {
        std::tuple<bool, bool, double, int> mergesplit_result = attempt_mergesplit_step(
            map_params, splitting_schedule, scoring_function,
            rng_state, sampling_space,
            plan, dummy_plan, 
            ust_sampler, tree_splitter,
            current_plan_multigraph,
            proposed_plan_multigraph,
            merge_prob_type, save_edge_selection_prob,
            current_plan_adj_region_pairs,
            current_plan_pair_unnoramalized_wgts,
            rho, is_final
        );
        // count tree size
        ++tree_sizes[std::get<3>(mergesplit_result)-1];
        // increase count if successful
        if(std::get<1>(mergesplit_result)){
            ++num_succesful_steps;
            ++successful_tree_sizes[std::get<3>(mergesplit_result)-1];
        }else if(std::get<0>(mergesplit_result)){
            ++successful_tree_sizes[std::get<3>(mergesplit_result)-1];
        }
    }
    if(DEBUG_MERGING_VERBOSE) Rprintf("Total success is %d\n", num_succesful_steps);
    return num_succesful_steps;
    
}





