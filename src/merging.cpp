/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2024/11
 * Purpose: Functions for Merging Plans
 ********************************************************/

#include "merging.h"

constexpr bool DEBUG_MERGING_VERBOSE = false; // Compile-time constant

/*
 *  Returns a sampler over a vector of adjacent pairs where the probability
 *  of a pair is decided according to `selection_type`
 *
 *  Current supported options are
 *      - uniform - Every pair has equal probability
 *      - district_pair - double district pairs have weight 1000, one district is 10,
 *          and two multidistricts have 1/(1+sum of their dvals)
 *
 *  @title Get Sampler over Adj Regions List
 *
 *  @param plan A plan object
 *  @param adj_pairs_and_boundary_lens A vector where each pair is
 *  (adj region1, adj region2, boundary length between 2 regions)
 *  @param selection_type A string controlling the function to use
 *  in assigning the unnormalized weight to each pair
 *
 *  @details No modifications to inputs made
 *
 *  @return A sampler where index i has probability proportional to the weight
 *  given to that pair
 */
arma::vec get_adj_pair_unnormalized_weights(
    Plan const &plan, std::vector<std::pair<RegionID, RegionID>> const &valid_region_adj_pairs,
    std::string const &selection_type) {
    // make a vector for the unnormalized weights
    arma::vec unnormalized_sampling_weights(valid_region_adj_pairs.size());

    // if uniform then just make all the weights 1 over the size
    if (selection_type == "uniform") {
        std::fill(unnormalized_sampling_weights.begin(), unnormalized_sampling_weights.end(),
                  static_cast<double>(valid_region_adj_pairs.size()));
    } else if (selection_type == "district_pair") {
        // if district pair then give huge weight to pairs of districts
        for (size_t i = 0; i < valid_region_adj_pairs.size(); i++) {
            int region1_id = valid_region_adj_pairs[i].first;
            int region2_id = valid_region_adj_pairs[i].second;
            int region1_dval = plan.region_sizes[region1_id];
            int region2_dval = plan.region_sizes[region2_id];

            // check if both are districts
            if (region1_dval == 1 && region2_dval == 1) {
                unnormalized_sampling_weights.at(i) = 1000.0;
            } else if (region1_dval == 1 || region2_dval == 1) {
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 10.0;
            } else {
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs
                unnormalized_sampling_weights.at(i) =
                    1 / (1 + static_cast<double>(region1_dval + region2_dval));
            }
        }

    } else if (selection_type == "multidistrict_pair") {
        // if multidistrict pair then give huge weight to pairs of districts
        for (size_t i = 0; i < valid_region_adj_pairs.size(); i++) {
            int region1_id = valid_region_adj_pairs[i].first;
            int region2_id = valid_region_adj_pairs[i].second;
            int region1_dval = plan.region_sizes[region1_id];
            int region2_dval = plan.region_sizes[region2_id];

            // check if both are districts
            if (region1_dval == 1 && region2_dval == 1) {
                unnormalized_sampling_weights.at(i) = .5;
            } else if (region1_dval == 1 || region2_dval == 1) {
                // else check if at least one is a district
                unnormalized_sampling_weights.at(i) = 1.0;
            } else {
                // if both multidistricts then do 1/1+(sum of two dvals)
                // This penalizes bigger pairs
                unnormalized_sampling_weights.at(i) = 1000.0;
            }
        }

    } else {
        throw Rcpp::exception("No valid adjacent pair sampler provided!");
    }
    // Return weights
    return unnormalized_sampling_weights;
}

// computes log metropolis hastings ratio
// new_region1_log_compactness, new_region2_log_compactness
// are pass by reference to save the log tau of the new regions
// if the split was successful
double get_log_mh_ratio(MapParams const &map_params, ScoringFunction const &scoring_function,
                        const int region1_id, const int region2_id,
                        double const current_log_eff_boundary_len,
                        double const proposed_log_eff_boundary_len,
                        const double log_current_pair_merge_prob,
                        const double log_new_pair_merge_prob,
                        double &new_region1_log_compactness,
                        double &new_region2_log_compactness, Plan const &current_plan,
                        Plan const &proposed_plan, double const rho, bool const is_final,
                        bool const using_caching, WeightCache *weight_cache,
                        GranularMCMCTimes &granular_times) {
    double log_mh_ratio = 0.0;
    // We add the scores of the current regions and plan and subtract the scores of the 
    // proposed regions and plan
    // first we'll start with the regions 
    auto region_score_time = maybe_now(); // optional timing 
    log_mh_ratio +=
        scoring_function.compute_region_soft_score(current_plan, region1_id, is_final);
    log_mh_ratio +=
        scoring_function.compute_region_soft_score(current_plan, region2_id, is_final);
    log_mh_ratio -=
        scoring_function.compute_region_soft_score(proposed_plan, region1_id, is_final);
    log_mh_ratio -=
        scoring_function.compute_region_soft_score(proposed_plan, region2_id, is_final);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.region_scores, region_score_time);
    }


    // now we do the plan scores 
    auto plan_score_time = maybe_now();
    log_mh_ratio += scoring_function.compute_plan_score(current_plan).second;
    log_mh_ratio -= scoring_function.compute_plan_score(proposed_plan).second;
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.plan_scores, plan_score_time);
    }
    
    // we compute taus if neccesary
    if (rho != 1) {
        auto log_st_time = maybe_now();
        // we add scorse of proposed region
        new_region1_log_compactness =
            (rho - 1) * proposed_plan.compute_log_region_spanning_trees(map_params, region1_id);
        new_region2_log_compactness =
            (rho - 1) * proposed_plan.compute_log_region_spanning_trees(map_params, region2_id);
        log_mh_ratio += new_region1_log_compactness;
        log_mh_ratio += new_region2_log_compactness;
        // we subtract scores of current region
        // This function return (rho-1)*log compactness
        log_mh_ratio -= compute_or_fetch_log_region_compactness(
            current_plan, region1_id, map_params, rho, using_caching, weight_cache);
        log_mh_ratio -= compute_or_fetch_log_region_compactness(
            current_plan, region2_id, map_params, rho, using_caching, weight_cache);
        if constexpr (perf_config::track_granular_times){
            add_elapsed(granular_times.tau_terms, log_st_time);
        }
    }
    // we add forward kernel proposed to current terms
    log_mh_ratio += current_log_eff_boundary_len;
    log_mh_ratio += log_new_pair_merge_prob;
    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("We added %f and %f, ", current_log_eff_boundary_len, log_new_pair_merge_prob);
    }
    // we subtract forward kernel current to proposed terms
    log_mh_ratio -= proposed_log_eff_boundary_len;
    log_mh_ratio -= log_current_pair_merge_prob;
    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("We subtracted %f and %f, ", proposed_log_eff_boundary_len,
                log_current_pair_merge_prob);
        Rprintf("Ratio is now %f\n", std::exp(log_mh_ratio));
    }

    return log_mh_ratio;
}

// runs a mergesplit Metropolis Hastings Step
std::tuple<bool, bool, double, int> attempt_mergesplit_step(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, RNGState &rng_state,
    SamplingSpace const sampling_space, Plan &plan, Plan &new_plan, USTSampler &ust_sampler,
    TreeSplitter &tree_splitter, PlanMultigraph &current_plan_multigraph,
    PlanMultigraph &proposed_plan_multigraph, std::string const merge_prob_type,
    bool save_edge_selection_prob, std::vector<std::pair<RegionID, RegionID>> &adj_region_pairs,
    arma::vec &unnormalized_pair_wgts, double const rho, bool const is_final, bool const do_mh,
    bool const using_caching, WeightCache *weight_cache,
    GranularMCMCTimes &granular_times
) {
    // sample a pair
    auto pair_selection_time = maybe_now(); // optional timing 
    int sampled_pair_index = rng_state.r_int_unnormalized_wgt(unnormalized_pair_wgts);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.selecting_merge_pair, pair_selection_time); // optional timing 
    }

    std::pair<int, int> merge_pair = adj_region_pairs[sampled_pair_index];

    int region1_id = merge_pair.first;
    int region2_id = merge_pair.second;
    int merged_region_size = plan.region_sizes[region1_id] + plan.region_sizes[region2_id];

    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("Picked pair (%d, %d)\n", region1_id, region2_id);
    }

    plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In attempt_mergesplit_step, calling on `plan` before wilson call"
    );
    std::string diff_msg;
    Tree const forest_graph_before = plan.get_forest_adj();
    // try to draw a region
    auto wilson_time = maybe_now(); // optional timing 
    std::tuple<bool, EdgeCut> edge_search_result =
        ust_sampler.attempt_to_find_valid_tree_mergesplit(rng_state, scoring_function,
                                                          tree_splitter, plan, region1_id,
                                                          region2_id, save_edge_selection_prob);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.wilson_time, wilson_time); // optional timing
    }
    if (!plan.forest_graph_equals_order_insensitive(forest_graph_before, diff_msg)) {
        throw std::runtime_error(
            "Current plan forest_graph changed wcalling attempt_to_find_valid_tree_mergesplit on plan.\n" +
            diff_msg
        );
    }
     plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In attempt_mergesplit_step, calling on `plan` after wilson call"
    );

    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("A Splitting Checkpoint 1.\n");
    }

    // If nothing drawn immediately return
    if (!std::get<0>(edge_search_result)) {
        if constexpr (DEBUG_MERGING_VERBOSE) {
            Rprintf("Failed!\n");
        }
        return std::make_tuple(false, false, std::log(0.0), merged_region_size);
    }

    // IN THE FUTURE CAN AVOID THE COPYING BY JUST TRAVERSING THE TREE
    // Just traverse tree and check if not in merged region or something
    // copy the new plan to be the old one
    auto initial_copy_time = maybe_now(); // optional timing 
    new_plan.shallow_copy(plan);
    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("A Splitting Checkpoint 1.5!\n");
    }
    new_plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In attempt_mergesplit_step, calling on `new_plan` after copying"
    );
    // now split that region we found on the old one
    new_plan.update_from_successful_split(tree_splitter, ust_sampler,
                                          std::get<1>(edge_search_result), region1_id,
                                          region2_id, false);
    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("A Splitting Checkpoint 2.\n");
    }
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.plan_copying, initial_copy_time); // optional timing 
    }
    new_plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In attempt_mergesplit_step, calling on `new_plan` after successful update"
    );

    if (!plan.forest_graph_equals_order_insensitive(forest_graph_before, diff_msg)) {
        throw std::runtime_error(
            "Current plan forest_graph changed by update_from_successful_split on NEW plan.\n" +
            diff_msg
        );
    }

    // check new plan is hierarchically valid if needed
    auto new_pair_building_time = maybe_now(); // optional timing 
    auto build_attempt = new_plan.attempt_to_get_valid_mergesplit_pairs(
        proposed_plan_multigraph, splitting_schedule, scoring_function, is_final);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.get_valid_pairs, new_pair_building_time); // optional timing 
    }

    new_plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In attempt_mergesplit_step, calling on `new_plan` after attempt_to_get_valid_mergesplit_pairs"
    );

    // new plan is valid if build attempt successful and passes any hard constraints
    auto hard_constraint_time = maybe_now();
    bool new_plan_valid =
        build_attempt.first &&
        scoring_function.new_split_ok(new_plan, region1_id, region2_id, is_final);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.hard_constraint_time, hard_constraint_time); // optional timing 
    }

    

    if constexpr (DEBUG_MERGING_VERBOSE) {
        Rprintf("%d county splits!\n", proposed_plan_multigraph.num_county_region_components);
    }
    if (!new_plan_valid) {
        if constexpr (DEBUG_MERGING_VERBOSE) {
            REprintf("%d splits for %d regions\n",
                     proposed_plan_multigraph.num_county_region_components, plan.num_regions);
        }
        // return failure
        return std::make_tuple(true, false, std::log(0.0), merged_region_size);
    }
    // get adj pairs
    auto new_valid_adj_region_pairs = build_attempt.second;
    // get the weights
    auto pair_weight_time = maybe_now();
    auto new_valid_pair_weights = get_adj_pair_unnormalized_weights(
        new_plan, new_valid_adj_region_pairs, merge_prob_type);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.selecting_merge_pair, pair_weight_time); // optional timing 
    }


    // compute mh ratio if needed
    double log_mh_ratio;
    bool proposal_accepted;

    if (do_mh) {
        int region_pair_proposal_index = -1;
        // Find the index of the pair in the proposed plan
        for (size_t i = 0; i < new_valid_adj_region_pairs.size(); i++) {

            if (new_valid_adj_region_pairs[i].first == region1_id &&
                new_valid_adj_region_pairs[i].second == region2_id) {
                region_pair_proposal_index = i;
                break;
            }
        }
        // means disconnected thing glitch. If encountered just print and ignore
        if (region_pair_proposal_index == -1)
        {
            std::ostringstream oss;

            oss << "Merged and split regions "
                << region1_id << " and " << region2_id
                << " no longer adjacent.\n";

            oss << "Current adj pairs: ";
            for (auto const &a_pair : adj_region_pairs) {
                oss << "(" << static_cast<int>(a_pair.first)
                    << ", " << static_cast<int>(a_pair.second) << "), ";
            }
            oss << "\n";

            oss << "Proposed adj pairs: ";
            for (auto const &a_pair : new_valid_adj_region_pairs) {
                oss << "(" << static_cast<int>(a_pair.first)
                    << ", " << static_cast<int>(a_pair.second) << "), ";
            }
            oss << "\n";

            oss << "Current Plan IDs:\n";
            for (auto const id : plan.region_ids) {
                oss << static_cast<int>(id) << ";";
            }
            oss << "\n";

            oss << "Proposed Plan IDs:\n";
            for (auto const id : new_plan.region_ids) {
                oss << static_cast<int>(id) << ";";
            }
            oss << "\n";

            oss << "Current plan debug:\n";
            oss << plan.debug_string(true);

            oss << "Proposed plan debug:\n";
            oss << new_plan.debug_string(true);

            throw std::runtime_error(oss.str());
        }
        if constexpr (DEBUG_MERGING_VERBOSE) {
            Rprintf("selected new pair index is %d!\n", region_pair_proposal_index);
        }
        bool const using_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
        // compute the boundary length
        auto eff_boundary_time = maybe_now();
        double current_log_eff_boundary = plan.get_log_eff_boundary_len(
            current_plan_multigraph, splitting_schedule, ust_sampler, tree_splitter,
            scoring_function, region1_id, region2_id);
        double proposed_log_eff_boundary = new_plan.get_log_eff_boundary_len(
            proposed_plan_multigraph, splitting_schedule, ust_sampler, tree_splitter,
            scoring_function, region1_id, region2_id);
        // If linking edge space we need to subtract linking edge correction term
        if (using_linking_edge_space) {
            // add instead of subtract bc its flipped in MH ratio
            current_log_eff_boundary += current_plan_multigraph.compute_log_multigraph_tau(
                plan.num_regions, scoring_function);
            proposed_log_eff_boundary += proposed_plan_multigraph.compute_log_multigraph_tau(
                new_plan.num_regions, scoring_function);
        }
        if constexpr (perf_config::track_granular_times){
            add_elapsed(granular_times.eff_boundary_length, eff_boundary_time); // optional timing
        }

        

        if constexpr (DEBUG_MERGING_VERBOSE) {
            Rprintf("Finding Adjacent regions %d, %d!\n", region1_id, region2_id);
            Rprintf("Current Plan: %d Adjacent Regions and I picked index %d ",
                    (int)adj_region_pairs.size(), sampled_pair_index);
            Rprintf("with probability %f\n",
                    std::exp(std::log(unnormalized_pair_wgts(sampled_pair_index)) -
                             std::log(arma::sum(unnormalized_pair_wgts))));
            Rprintf("Proposed Plan: %d Adjacent Regions and I picked index %d ",
                    (int)new_valid_adj_region_pairs.size(), region_pair_proposal_index);
            Rprintf("with probability %f\n",
                    std::exp(std::log(new_valid_pair_weights(region_pair_proposal_index)) -
                             std::log(arma::sum(new_valid_pair_weights))));
        }

        double new_region1_log_compactness, new_region2_log_compactness;

        log_mh_ratio =
            get_log_mh_ratio(map_params, scoring_function, region1_id, region2_id,
                             current_log_eff_boundary, proposed_log_eff_boundary,
                             std::log(unnormalized_pair_wgts(sampled_pair_index)) -
                                 std::log(arma::sum(unnormalized_pair_wgts)),
                             std::log(new_valid_pair_weights(region_pair_proposal_index)) -
                                 std::log(arma::sum(new_valid_pair_weights)),
                             new_region1_log_compactness, new_region2_log_compactness, plan,
                             new_plan, rho, is_final, using_caching, weight_cache,
                             granular_times);
        proposal_accepted = rng_state.r_unif() <= std::exp(log_mh_ratio);
        if constexpr (DEBUG_MERGING_VERBOSE)
            Rprintf("Ratio is %f and it is ", std::exp(log_mh_ratio));

        // update the cached compactness if needed
        if (proposal_accepted && using_caching) {
            weight_cache->region_cache_values[region1_id] = new_region1_log_compactness;
            weight_cache->this_plan_order_added[region1_id] =
                new_plan.region_added_order[region1_id];
            weight_cache->region_cache_values[region2_id] = new_region2_log_compactness;
            weight_cache->this_plan_order_added[region2_id] =
                new_plan.region_added_order[region2_id];
        }

    } else {
        // always accept
        proposal_accepted = true;
        log_mh_ratio = 0;
    }

    if (proposal_accepted) {
        if constexpr (DEBUG_MERGING_VERBOSE)
            Rprintf("ACCEPTED!! Now updating plans\n");
        new_plan.check_forest_integrity(
            map_params.graph_edge_index,
            "In attempt_mergesplit_step, calling on `new_plan` after proposal accepted"
        );
        // if successful then actually update
        auto final_copy_time = maybe_now(); // optional timing 
        plan.shallow_copy(new_plan);
        if constexpr (perf_config::track_granular_times){
            add_elapsed(granular_times.plan_copying, final_copy_time); // optional timing
        }
        if constexpr (DEBUG_MERGING_VERBOSE)
            Rprintf("Plan updated, now swapping multigraphs\n");

        plan.check_forest_integrity(
            map_params.graph_edge_index,
            "In attempt_mergesplit_step, calling on `plan` after proposal accepted and new plan copied"
        );

        // swap the plan multigraphs
        swap_plan_multigraphs(current_plan_multigraph, proposed_plan_multigraph);
        // update pairs and weights to be current one
        unnormalized_pair_wgts = new_valid_pair_weights;
        adj_region_pairs = new_valid_adj_region_pairs;
        return std::make_tuple(true, true, log_mh_ratio, merged_region_size);
    } else {
        if constexpr (DEBUG_MERGING_VERBOSE)
            Rprintf("REJECTED!! (mh_ratio=%g)\n", std::exp(log_mh_ratio));
        // else reject and do nothing
        return std::make_tuple(true, false, log_mh_ratio, merged_region_size);
    }
}

// runs num_steps_to_run merge split stesp
int run_merge_split_steps(MapParams const &map_params,
                          const SplittingSchedule &splitting_schedule,
                          ScoringFunction const &scoring_function, RNGState &rng_state,
                          SamplingSpace const sampling_space, Plan &plan, Plan &dummy_plan,
                          USTSampler &ust_sampler, TreeSplitter &tree_splitter,
                          PlanMultigraph &current_plan_multigraph,
                          PlanMultigraph &proposed_plan_multigraph,
                          std::string const merge_prob_type, double const rho,
                          bool const is_final, int num_steps_to_run,
                          std::vector<int> &tree_sizes, std::vector<int> &successful_tree_sizes,
                          bool const using_caching, WeightCache *weight_cache,
                          GranularMCMCTimes &granular_times) {
    plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In run_merge_split_steps, calling on `plan` before any code"
    );
    int num_succesful_steps = 0;
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;

    // Build the multigraph and get pairs of adj districts
    auto pair_building_time = maybe_now(); // optional timing 
    auto current_plan_adj_region_pairs =
        plan.attempt_to_get_valid_mergesplit_pairs(current_plan_multigraph, splitting_schedule,
                                                   scoring_function, is_final)
            .second;
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.get_valid_pairs, pair_building_time); // optional timing
    }
    plan.check_forest_integrity(
        map_params.graph_edge_index,
        "In run_merge_split_steps, calling on `plan` after initial multigraph built"
    );

    auto pair_weight_time = maybe_now();
    arma::vec current_plan_pair_unnoramalized_wgts =
        get_adj_pair_unnormalized_weights(plan, current_plan_adj_region_pairs, merge_prob_type);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.selecting_merge_pair, pair_weight_time); // optional timing 
    }

    // run the merge split steps and count success
    for (size_t i = 0; i < num_steps_to_run; i++) {


    std::string const plan_msg =
        "In run_merge_split_steps, for loop i = " + std::to_string(i) +
        ", calling on `plan`";
        plan.check_forest_integrity(
            map_params.graph_edge_index,
            plan_msg
        );


        std::tuple<bool, bool, double, int> mergesplit_result = attempt_mergesplit_step(
            map_params, splitting_schedule, scoring_function, rng_state, sampling_space, plan,
            dummy_plan, ust_sampler, tree_splitter, current_plan_multigraph,
            proposed_plan_multigraph, merge_prob_type, save_edge_selection_prob,
            current_plan_adj_region_pairs, current_plan_pair_unnoramalized_wgts, rho, is_final,
            true, using_caching, weight_cache, granular_times);
        // count tree size
        ++tree_sizes[std::get<3>(mergesplit_result) - 1];
        // increase count if successful
        if (std::get<1>(mergesplit_result)) {
            ++num_succesful_steps;
            ++successful_tree_sizes[std::get<3>(mergesplit_result) - 1];
        } else if (std::get<0>(mergesplit_result)) {
            ++successful_tree_sizes[std::get<3>(mergesplit_result) - 1];
        }
    }
    if constexpr (DEBUG_MERGING_VERBOSE)
        Rprintf("Total success is %d\n", num_succesful_steps);

    std::string const after_plan_msg =
        "In run_merge_split_steps, after all steps, calling on `dummy_plan`";
        plan.check_forest_integrity(
            map_params.graph_edge_index,
            after_plan_msg
    );

    return num_succesful_steps;
}
