/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2024/11
 * Purpose: SMC weight calculation related functions
 ********************************************************/

constexpr bool DEBUG_WEIGHTS_VERBOSE = false; // Compile-time constant
#include "weights.h"

#include "advanced_types.h"          // MapParams, SamplingSpace
#include "base_plan_type.h"          // Plan, PlanMultigraph
#include "redist_constants.h"        // perf_config
#include "scoring.h"                 // ScoringFunction
#include "splitting_schedule_types.h" // SplittingSchedule
#include "tree_splitting.h"          // TreeSplitter
#include "utils.h"                   // maybe_now, add_elapsed
#include "weight_caching.h"          // WeightCache, compute_or_fetch_log_region_compactness
#include "wilson.h"                  // USTSampler

#include <algorithm> // std::max_element, std::max
#include <cmath>     // std::exp, std::log, std::isfinite
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <tuple>     // std::get
#include <utility>   // std::pair, std::make_pair
#include <vector>





/*
 *  log(sum(exp(log_terms))), computed by factoring out the largest term so
 *  that every exponent is <= 0. Summing exp(t) directly overflows past ~709
 *  and flushes to zero below ~-745; after the shift the terms lie in (0, 1]
 *  and the sum in [1, n], so neither can happen.
 *  
 *  Math is 
 *  Σ exp(t_i)  =  exp(M) · Σ exp(t_i − M)          where M = max t_i
 *  log Σ exp(t_i)  =  M + log Σ exp(t_i − M)
 */
static double log_sum_exp(std::vector<double> const &log_terms) {
    // An empty sum is 0 and log(0) is -Inf. Guarding here also keeps
    // `max_element` from dereferencing `end()`, which for an empty vector is
    // a null pointer whenever nothing was ever allocated into it.
    if (log_terms.empty()) return -std::numeric_limits<double>::infinity();

    double const max_term = *std::max_element(log_terms.begin(), log_terms.end());
    // If its not finite then just return
    if (!std::isfinite(max_term)) return max_term;

    double sum = 0.0;
    for (double const term : log_terms) sum += std::exp(term - max_term);
    return max_term + std::log(sum);
}


// computes the backwards kernel that is uniform in
// the number of ancestors
double compute_simple_log_incremental_weight(Plan const &plan, PlanMultigraph &plan_multigraph,
                                             const SplittingSchedule &splitting_schedule,
                                             USTSampler &ust_sampler,
                                             TreeSplitter &edge_splitter,
                                             SamplingSpace const sampling_space,
                                             ScoringFunction const &scoring_function,
                                             double rho, bool compute_log_splitting_prob,
                                             double const multidistrict_selection_alpha,
                                             bool is_final_split, bool const using_caching,
                                             WeightCache *weight_cache,
                                             GranularWeightTimes &granular_times) {
    // bool for whether we'll need to compute spanning tree count
    bool compute_log_tau = rho != 1;
    // find the index of the two newly split regions
    auto most_recent_split_regions = plan.get_most_recently_split_regions();
    auto region1_id = most_recent_split_regions.first;
    auto region2_id = most_recent_split_regions.second;

    if constexpr (DEBUG_WEIGHTS_VERBOSE) {
        Rprintf("The two regions are %d and %d\n", region1_id, region2_id);
        Rcpp::Rcerr << std::flush;
    }

    // build the schedule context once for this plan
    ScheduleContext schedule_ctx;
    splitting_schedule.build_context(plan, schedule_ctx);

    // build the plan multigraph and get valid adj pairs
    auto time_get_pairs = maybe_now(); // optional timing
    auto valid_adj_pairs = plan.get_valid_smc_merge_regions(plan_multigraph, splitting_schedule,
                                                            schedule_ctx,
                                                            scoring_function, is_final_split);
    if constexpr (perf_config::track_granular_times) {
        add_elapsed(granular_times.get_valid_pairs, time_get_pairs); // optional timing
    }

    // compute forward and backwards kernel term
    double log_backwards_kernel_term = 0.0;
    double log_forward_kernel_term = 0.0;
    double log_tau_ratio_term = 0.0;

    auto time_log_tau = maybe_now(); // optional timing

    if (sampling_space == SamplingSpace::GraphSpace) {
        // the number of valid adjacent regions is just the number of elements in the map
        const int num_valid_adj_region_pairs = valid_adj_pairs.size();
        log_backwards_kernel_term -= std::log(static_cast<double>(num_valid_adj_region_pairs));
        // taus cancel so just add the boundary length which we have from hash map
        log_forward_kernel_term = plan.get_log_eff_boundary_len(
            plan_multigraph, splitting_schedule, schedule_ctx, ust_sampler, edge_splitter,
            scoring_function,
            region1_id, region2_id);
        if (compute_log_tau) {
            // do merged region tau if neccesary
            log_tau_ratio_term -=
                (rho - 1) * plan.compute_log_merged_region_spanning_trees(
                                plan_multigraph.map_params, region1_id, region2_id);
        }
    } else if (sampling_space == SamplingSpace::ForestSpace ||
               sampling_space == SamplingSpace::LinkingEdgeSpace) {
        // Sum of spanning trees, collected in log space and combined with
        // `log_sum_exp` for numerical stability.
        std::vector<double> log_spanning_tree_terms;
        log_spanning_tree_terms.reserve(valid_adj_pairs.size());
        double last_split_merge_log_tau = 0.0;

        for (auto const &region_pair : valid_adj_pairs) {
            int pair_region1_id = region_pair.first;
            int pair_region2_id = region_pair.second;

            double const log_tau = plan.compute_log_merged_region_spanning_trees(
                plan_multigraph.map_params, pair_region1_id, pair_region2_id);

            // If its the merged region actually split save it for use in forward kernel
            if (pair_region1_id == region1_id && pair_region2_id == region2_id) {
                last_split_merge_log_tau = log_tau;
            }

            log_spanning_tree_terms.push_back(log_tau);
        }

        log_backwards_kernel_term -= log_sum_exp(log_spanning_tree_terms);

        // One over tau of merged region
        log_forward_kernel_term -= last_split_merge_log_tau;
        // tree boundary length
        log_forward_kernel_term += plan.get_log_eff_boundary_len(
            plan_multigraph, splitting_schedule, schedule_ctx, ust_sampler, edge_splitter,
            scoring_function,
            region1_id, region2_id);

        if (compute_log_tau) {
            log_tau_ratio_term -= (rho - 1) * last_split_merge_log_tau;
        }
    }
    if (compute_log_tau) {
        // These are the two regions the last split created, so their order added
        // counters were both just bumped and the cache can never be fresh for
        // them. The fetch is still worth making: it writes the values back, so a
        // following mergesplit round reads them instead of recomputing.
        log_tau_ratio_term += compute_or_fetch_log_region_compactness(
            plan, region1_id, plan_multigraph.map_params, rho, using_caching, weight_cache);
        log_tau_ratio_term += compute_or_fetch_log_region_compactness(
            plan, region2_id, plan_multigraph.map_params, rho, using_caching, weight_cache);
    }

    if constexpr (perf_config::track_granular_times) {
        add_elapsed(granular_times.tau_terms, time_log_tau); // optional timing
    }

    double log_splitting_prob = 0;
    // for one district split the probability that region was chosen to be split is always 1
    if (compute_log_splitting_prob) {
        // in generalized region split find probability you would have
        // picked to split the union of the the two regions
        auto time_split_prob = maybe_now(); // optional timing
        log_splitting_prob = splitting_schedule.log_merged_region_selection_prob(
            plan, schedule_ctx, region1_id, region2_id, multidistrict_selection_alpha);
        if constexpr (perf_config::track_granular_times) {
            add_elapsed(granular_times.splitting_prob, time_split_prob); // optional timing
        }
        if constexpr (DEBUG_WEIGHTS_VERBOSE)
            Rprintf("Computed split prob %f\n", std::exp(log_splitting_prob));
    }

    double region1_score, region2_score, merged_region_score;

    region1_score = region2_score = merged_region_score = 0;

    // compute if any constraints
    if (scoring_function.any_soft_region_constraints) {
        auto time_region_score = maybe_now(); // optional timing
        // compute scoring functions
        region1_score =
            scoring_function.compute_region_soft_score(plan, region1_id);
        region2_score =
            scoring_function.compute_region_soft_score(plan, region2_id);
        merged_region_score =
            scoring_function
                .compute_merged_region_full_score(plan, region1_id, region2_id)
                .second;
        if constexpr (perf_config::track_granular_times) {
            add_elapsed(granular_times.region_scores, time_region_score); // optional timing
        }
        if constexpr (DEBUG_WEIGHTS_VERBOSE) {
            REprintf("Region (%d,%d) Scores (%f, %f) | Merged Score %f \n", region1_id,
                     region2_id, region1_score, region2_score, merged_region_score);
        }
    }

    double plan_score, prev_plan_score;
    plan_score = prev_plan_score = 0.0;

    if (scoring_function.any_soft_plan_constraints) {
        auto time_plan_score = maybe_now(); // optional timing
        plan_score = scoring_function.compute_plan_score(plan).second;
        prev_plan_score =
            scoring_function
                .compute_merged_plan_score(plan, region1_id, region2_id)
                .second;
        if constexpr (perf_config::track_granular_times) {
            add_elapsed(granular_times.plan_scores, time_plan_score); // optional timing
        }
        if constexpr (DEBUG_WEIGHTS_VERBOSE) {
            REprintf("Entire Plan Score %f | Previous Plan Score %f \n", plan_score,
                     prev_plan_score);
        }
    }

    double log_extra_plan_terms = 0.0;
    double log_extra_prev_plan_terms = 0.0;
    // if linking edge space we also need to correct for that
    auto time_multigraph_tau = maybe_now(); // optional timing
    if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
        // std::vector<int> merge_index_reshuffle(plan.num_regions);
        // we divide target by number of linking edges so
        // subtract log linking edges from denominator
        log_extra_plan_terms -=
            plan_multigraph.compute_log_multigraph_tau(plan.num_regions, scoring_function);

        if (plan.num_regions > 2) {
            double merged_tau = plan_multigraph.compute_merged_log_multigraph_tau(
                plan.num_regions, region1_id, region2_id, scoring_function);
            log_extra_prev_plan_terms -= merged_tau;
        }
    }
    // The optimal weights get their multigraph taus inside
    // `get_valid_adj_regions_and_eff_log_boundary_lens`, so they are counted as
    // pair time rather than tau time. Do the same here so the two weight types
    // stay comparable.
    if constexpr (perf_config::track_granular_times) {
        add_elapsed(granular_times.get_valid_pairs, time_multigraph_tau); // optional timing
    }

    // The weight is
    //      - Numerator: backwards kernel * e^-(J(new_region1)+J(new_region2)) * e^-J(new plan)
    //      - Denominator: multid_selection_prob * forward kernel term * e^-J(old_region) *
    //      e^-J(old plan)
    // So
    //      - log numerator: backwards kernel -J(new_region1)+J(new_region2))
    //      - log Denominator: log(multid_selection_prob) + log(forward kernel) - J(old_region)

    const double log_numerator = log_backwards_kernel_term -
                                 (region1_score + region2_score + plan_score) +
                                 log_extra_plan_terms;
    const double log_denom = log_forward_kernel_term + log_splitting_prob -
                             (merged_region_score + prev_plan_score) +
                             log_extra_prev_plan_terms;

    if constexpr (DEBUG_WEIGHTS_VERBOSE) {
        int region1_size = plan.region_sizes[region1_id];
        int region2_size = plan.region_sizes[region2_id];
        Rprintf("Doing (%d,%d) - sizes (%d, %d): forward %f, backward %f, split prob %f, ratio "
                "%f!\n",
                region1_id, region2_id, region1_size, region2_size,
                std::exp(log_forward_kernel_term), std::exp(log_backwards_kernel_term),
                std::exp(log_splitting_prob), std::exp(log_numerator - merged_region_score));

        Rprintf("Numerator - %f * %f * %f\n", std::exp(region1_score), std::exp(region2_score),
                std::exp(log_backwards_kernel_term));

        Rprintf("Denominator - %f * %f * %f \n", std::exp(log_splitting_prob),
                std::exp(log_forward_kernel_term), std::exp(merged_region_score));
    }

    // now just take the fraction
    double incremental_weight = log_numerator - log_denom + log_tau_ratio_term;
    // DONT KNOW WHY BUT VALIDATION ONLY WORKS IF YOU SUBTRACT LOG TAU
    // THEORY SAYS IT SHOULD BE ADD SO NEED TO INVESTIGATE

    // Check its not infinity
    if (!std::isfinite(incremental_weight)) {
        plan.Rprint(true);
        // plan_multigraph.Rprint();
        throw std::runtime_error(
            "One of the plan incremental weights is not finite!"
            "Try checking if constraint strength is too large and causing overflow errors.\n");
    }

    return incremental_weight;
}



// eventually need to modify to allow presaved options
// OLD DOCUMENTATION FROM GRAPH THING NEED TO UPDATE
// Compute the optimal log incremental weight of a plan
//
// Given a plan object this computes the minimum variance weights as derived in
// <PAPER NAME HERE>. This is equal to the inverse of a sum over all
// adjacent regions in a plan if using generalized region split and
// sum over all districts adajacent to the remainder if using one district
// split approach.
//
// @title Compute Optimal Incremental Weight of a plan
//
// @param g The underlying map graph
// @param plan A plan object
// @param split_district_only whether or not to compute the weights under
// the district only split scheme or not.
// @param target The target population for a single district
// @param pop_temper The population tempering parameter
//
// @details No modifications to inputs made
//
// @return the log of the incremental weight of the plan
//
double compute_log_optimal_incremental_weights(
    Plan const &plan, PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, USTSampler &ust_sampler,
    TreeSplitter &edge_splitter, SamplingSpace const sampling_space,
    ScoringFunction const &scoring_function, double const rho,
    double const whole_map_compactness_term, bool compute_log_splitting_prob,
    double const multidistrict_selection_alpha,
    bool is_final_split, bool const using_caching, WeightCache *weight_cache,
    GranularWeightTimes &granular_times) {
    // plan.Rprint();
    // bool for whether we'll need to compute spanning tree count
    bool const compute_log_tau = rho != 1;



    if constexpr (DEBUG_WEIGHTS_VERBOSE)
        Rprintf("Getting Pairs!\n");

    // build the schedule context once for this plan
    ScheduleContext schedule_ctx;
    splitting_schedule.build_context(plan, schedule_ctx);

    // get region pair to effective boundary length map
    auto time_get_pairs = maybe_now(); // optional timing 
    auto region_pair_log_eff_boundary_map =
        plan.get_valid_adj_regions_and_eff_log_boundary_lens(
            plan_multigraph, splitting_schedule, schedule_ctx, scoring_function, is_final_split,
            ust_sampler,
            edge_splitter);
    if constexpr (perf_config::track_granular_times){
        add_elapsed(granular_times.get_valid_pairs, time_get_pairs); // optional timing 
    }
    // vector to store log of each term in the weights before exponentiating and 
    // summing for numerical stability
    std::vector<double> pair_log_terms;
    pair_log_terms.reserve(region_pair_log_eff_boundary_map.size());

    // iterate over the pairs
    if constexpr (DEBUG_WEIGHTS_VERBOSE) {
        REprintf("There are %zu adjacent pairs!\n", region_pair_log_eff_boundary_map.size());
    }

    // compute plan score if needed
    double plan_score = 0.0;

    if (scoring_function.any_soft_plan_constraints) {
        auto time_plan_score = maybe_now(); // optional timing 
        plan_score += scoring_function.compute_plan_score(plan).second;
        if constexpr (perf_config::track_granular_times){
            add_elapsed(granular_times.plan_scores, time_plan_score); // optional timing 
        }

        if constexpr (DEBUG_WEIGHTS_VERBOSE) {
            REprintf("Entire Plan Score %f \n", plan_score);
        }
    }

    // Now iterate over adjacent region pairs and add splitting and pop temper
    for (const auto &pair_tuple : region_pair_log_eff_boundary_map) {
        double log_of_sum_term = 0.0;

        const int region1_id = std::get<0>(pair_tuple); // get the smaller region id
        const int region2_id = std::get<1>(pair_tuple); // get the bigger region id
        const double eff_log_boundary_len =
            std::get<2>(pair_tuple); // get the effective boundary length

        if constexpr (DEBUG_WEIGHTS_VERBOSE) {
            REprintf("Pair (%d,%d) - %f\n", region1_id, region2_id,
                     std::exp(eff_log_boundary_len));
        }

        // add to term
        log_of_sum_term += eff_log_boundary_len;

        // If not compute_log_splitting_prob then its just log(1) = 0
        if (compute_log_splitting_prob) {
            // in generalized region split find probability you would have
            // picked to split the union of the the two regions
            auto time_split_prob = maybe_now(); // optional timing 
            double log_splitting_prob = splitting_schedule.log_merged_region_selection_prob(
                plan, schedule_ctx, region1_id, region2_id, multidistrict_selection_alpha);
    
            if constexpr (perf_config::track_granular_times){
                add_elapsed(granular_times.splitting_prob, time_split_prob); // optional timing 
            }

            log_of_sum_term += log_splitting_prob;
            if constexpr (DEBUG_WEIGHTS_VERBOSE)
                REprintf("\n");
        }

        // compute score ratio if any constraints
        if (scoring_function.any_soft_region_constraints) {
            auto time_region_score = maybe_now(); // optional timing 
            // compute scoring functions
            const double region1_score =
                scoring_function.compute_region_soft_score(plan, region1_id);
            const double region2_score =
                scoring_function.compute_region_soft_score(plan, region2_id);
            const double merged_region_score =
                scoring_function
                    .compute_merged_region_full_score(plan, region1_id, region2_id
                                                      )
                    .second;
             
            if constexpr (perf_config::track_granular_times){
                add_elapsed(granular_times.region_scores, time_region_score); // optional timing
            }
            if constexpr (DEBUG_WEIGHTS_VERBOSE) {
                REprintf("Region (%u, %u) Scores (%f, %f) | Merged Score %f \n", region1_id,
                         region2_id, region1_score, region2_score, merged_region_score);
            }
            // log ratio is (log region1 + log region2) - log score merged
            double const score_ratio = region1_score + region2_score - merged_region_score;
            log_of_sum_term += score_ratio;
        }

        double merged_plan_score = 0.0;
        if (scoring_function.any_soft_plan_constraints) {
            auto time_plan_score = maybe_now(); // optional timing 
            merged_plan_score =
                scoring_function
                    .compute_merged_plan_score(plan, region1_id, region2_id)
                    .second; 
            if constexpr (perf_config::track_granular_times){
                add_elapsed(granular_times.region_scores, time_plan_score); // optional timing
            }
            if constexpr (DEBUG_WEIGHTS_VERBOSE) {
                REprintf("For Regions (%u, %u) Merged Entire Plan Score %f \n", region1_id,
                         region2_id, merged_plan_score);
            }
        }
        // term is e^- score(current plan) / e^- score(merged plan)
        // so we add the current plan score and subtract the log score
        // meaning add the difference
        log_of_sum_term += (plan_score - merged_plan_score);

        if constexpr (DEBUG_WEIGHTS_VERBOSE) {
            int region1_size = plan.region_sizes[region1_id];
            int region2_size = plan.region_sizes[region2_id];
            Rprintf("Adding (%d,%d) - sizes (%d, %d): len %f, split prob %f, ratio %f!\n",
                    region1_id, region2_id, region1_size, region2_size,
                    std::exp(eff_log_boundary_len), std::exp(log_of_sum_term), log_of_sum_term);
        }

        // Do taus if neccesary
        if (compute_log_tau) {
            auto time_log_tau = maybe_now(); // optional timing 
            if constexpr (DEBUG_WEIGHTS_VERBOSE) {
                REprintf("Computing log tau ratio!\n");
            }
            double log_tau_ratio = 0.0;
            // add merged region
            if (plan.num_regions == 2) {
                // if its two regions its just the tau for the whole map
                log_tau_ratio += whole_map_compactness_term;
            } else {
                log_tau_ratio +=
                    (rho - 1) * plan.compute_log_merged_region_spanning_trees(
                                    plan_multigraph.map_params, region1_id, region2_id);
            }

            if constexpr (DEBUG_WEIGHTS_VERBOSE) {
                REprintf("Computing log tau for each region pair!\n");
            }

            // subtract split regions
            // This function return (rho-1)*log compactness
            log_tau_ratio -= compute_or_fetch_log_region_compactness(
                plan, region1_id, plan_multigraph.map_params, rho, using_caching, weight_cache);
            log_tau_ratio -= compute_or_fetch_log_region_compactness(
                plan, region2_id, plan_multigraph.map_params, rho, using_caching, weight_cache);
            log_of_sum_term += log_tau_ratio;

            if constexpr (perf_config::track_granular_times){
                add_elapsed(granular_times.tau_terms, time_log_tau); // optional timing 
            }
        }

        pair_log_terms.push_back(log_of_sum_term);

    }

    // The incremental weight is the inverse of the sum of terms so 
    // negative log sum
    double const log_incremental_weight = -log_sum_exp(pair_log_terms);


    if constexpr (DEBUG_WEIGHTS_VERBOSE) {
        REprintf("Weight=%f, log weight = %f\n", std::exp(log_incremental_weight),
                 log_incremental_weight);
    }

    if (!std::isfinite(log_incremental_weight)) {
        plan.Rprint(true);
        if (log_incremental_weight == std::numeric_limits<double>::infinity()) {
            throw std::runtime_error(
                "One of the plan incremental weights is infinite! No adjacent region pair "
                "was a valid merge, so the plan has no viable ancestors. Check whether hard "
                "constraints are actually splittable.\n");
        } else if (log_incremental_weight == -std::numeric_limits<double>::infinity()) {
            throw std::runtime_error(
                "One of the plan incremental weights is 0! Some term in the sum was infinite. "
                "Check whether a constraint strength is large enough to make a score infinite.\n");
        }
        throw std::runtime_error(
            "One of the plan incremental weights is NaN! This usually means two terms in the "
            "sum were infinite with opposite signs.\n");
    }


    // now return the log of the inverse of the sum
    return log_incremental_weight;
}



std::pair<bool, double>
WeightCache::attempt_to_get_region_value(RegionID const region_id,
                                         int const current_region_order_added_num) {
    // If the current region order added is the same as the stored one then the region
    // has not changed and the old value is safe. If its not the same then the value is
    // stale and must be removed
    if (current_region_order_added_num == this_plan_order_added[region_id]) {
        return std::make_pair(true, region_cache_values[region_id]);
    } else {
        // the value is stale
        return std::make_pair(false, 0);
    }
};
