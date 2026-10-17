/********************************************************
 * Author: Cory McCartan
 * Institution: Harvard University
 * Date Created: 2021/02
 * Purpose: Merge-split MCMC redistricting sampler
 * (like Carter et al. 2019 but using the SMC proposal)
 ********************************************************/


#include <iomanip>
#include <locale>

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <Rcpp.h>

#include "base_plan_type.h"
#include "graph_plan_type.h"
#include "merging.h"
#include "redist_alg_helpers.h"
#include "scoring.h"
#include "tree_splitting.h"
#include "splitting_schedule_types.h"
#include "random.h"
#include "wilson.h"
#include "utils.h"
#include "threading_helpers.h"


#include <RcppThread.h>
#include <cli/progress.h>


constexpr bool DEBUG_PURE_MS_VERBOSE = false; // Compile-time constant

/*
 * Main entry point.
 *
 * USING MCMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
Rcpp::List ms_plans(
    int const nsims, std::int64_t const warmup, int const thin, int const ndists, int const total_seats,
    Rcpp::IntegerVector const &district_seat_sizes, Rcpp::List const &adj_list,
    const Rcpp::IntegerVector &counties, const Rcpp::IntegerVector &pop, double const target, double const lower,
    double const upper,
    double const rho, // compactness
    Rcpp::IntegerMatrix const &init_plan, Rcpp::IntegerMatrix const &init_seats,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &pair_rule, // method for setting probability of picking a pair to merge
    Rcpp::List const &control,     // control has pop temper, and k parameter value, and splitting schedule
    Rcpp::List const &constraints, // constraints
    int const verbosity = 3, bool const diagnostic_mode = false) {
    // whether or not to perform MH step
    bool do_mh = (bool) control["do_mh"];
    // whether or not to sample hierarchically but not enforce hierarchical validity
    bool enforce_hierarchical = Rcpp::as<bool>(control["enforce_hierarchical"]);

    if constexpr (DEBUG_PURE_MS_VERBOSE)
        Rprintf("Checkpoint 1!\n");
    Rcpp::List out; // return
    // re-seed MT so that `set.seed()` works in R
    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    RNGState rng_state(global_rng_seed);
    // Set the sampling space
    SamplingSpace sampling_space = get_sampling_space(sampling_space_str);
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // TODO: Legacy, in future remove
    global_seed_rng( (int) Rcpp::sample(INT_MAX, 1)[0]);

    // make sure both are 1 column matrix
    if (init_plan.ncol() > 1 || init_seats.ncol() > 1) {
        std::ostringstream oss;
        oss << "Mergesplit takes a single initial plan, so `init_plan` and `init_seats` must "
            << "each have exactly 1 column, but `init_plan` has " << init_plan.ncol()
            << " columns and `init_seats` has " << init_seats.ncol() << " columns.\n";
        throw Rcpp::exception(oss.str().c_str());
    }

    // Create map level graph and county level multigraph
    MapParams const map_params(list_to_graph(adj_list), 
        Rcpp::as<std::vector<unsigned int>>(counties), 
        Rcpp::as<std::vector<unsigned int>>(pop), 
        ndists, total_seats,
        Rcpp::as<std::vector<int>>(district_seat_sizes), lower, target, upper,
        sampling_space, !enforce_hierarchical);
    int V = map_params.g.size();

    int initial_num_regions = static_cast<int>(ndists);
    if (init_seats.nrow() != initial_num_regions) {
        REprintf("Inferred %d Initial Regions but Region Sizes Matrix has %u columns!\n",
                 initial_num_regions, init_seats.nrow());
        throw Rcpp::exception("Initial plan region sizes should match number of regions");
    }
    if constexpr (DEBUG_PURE_MS_VERBOSE)
        Rprintf("Checkpoint 2!\n");

    // get splitting type
    SplittingMethodType splitting_method =
        get_splitting_type(static_cast<std::string>(control["splitting_method"]));

    auto splitting_schedule_ptr = std::make_unique<PureMSSplittingSchedule>(
        ndists, total_seats, Rcpp::as<std::vector<int>>(district_seat_sizes));
    // splitting_schedule_ptr->print_current_step_splitting_info();
    bool const mmd_plans = map_params.is_mmd;

    // Add scoring function (constraints)
    // No population tempering
    ScoringFunction const scoring_function(map_params, constraints, 0, false);

    // Now create diagnostic information
    Rcpp::NumericVector log_mh_ratios(nsims); // stores log mh ratio
    Rcpp::IntegerMatrix saved_plans_mat(V, nsims);
    Rcpp::IntegerMatrix saved_district_pops_mat(ndists, nsims);
    Rcpp::IntegerMatrix saved_plan_sizes(mmd_plans ? ndists : 1, mmd_plans ? nsims : 1);
    int current_plan_mat_col = 0;
    std::vector<int> tree_sizes(total_seats, 0);
    std::vector<int> successful_tree_sizes(total_seats, 0);

    // Level 3
    // Saves proposal plans
    Rcpp::IntegerMatrix proposed_plans_mat(diagnostic_mode ? V : 1,
                                           diagnostic_mode ? nsims : 1);

    std::vector<std::vector<Graph>> all_steps_forests_adj_list;
    all_steps_forests_adj_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::ForestSpace) ? nsims : 0);

    if constexpr (DEBUG_PURE_MS_VERBOSE)
        Rprintf("Checkpoint 3!\n");

    // Create current plan and proposal

    int global_rng_seed2 = (int)Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;
    rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);

    Rcpp::IntegerVector mh_decisions(nsims);
    double mha;


    
    // Track the total number of successes during the warmup
    std::int64_t warmup_acceptances = 0;
    // Track total number of successes after warmup
    std::int64_t post_warump_acceptances = 0;

    std::int64_t const total_post_warmup_steps =
        static_cast<std::int64_t>(nsims) *
        static_cast<std::int64_t>(thin);

    std::int64_t const total_steps =
        total_post_warmup_steps + warmup;

    std::int64_t const start = 1 - warmup;

    // Track timing 
    double warmup_time;
    double nonwarmup_time;
    // granular timing if that's being tracked 
    GranularMCMCTimes granular_times;


    {
        USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
        PlanMultigraph current_plan_multigraph(map_params, sampling_space ==
                                                               SamplingSpace::LinkingEdgeSpace
                                                );
        PlanMultigraph proposed_plan_multigraph(
            map_params, sampling_space == SamplingSpace::LinkingEdgeSpace
            );

        // 1 means serial under the package threading convention
        RcppThread::ThreadPool pool = get_thread_pool(1);
        // underlying vector from plan
        PlanEnsemble plan_ensemble = get_plan_ensemble(
            map_params, *splitting_schedule_ptr, initial_num_regions, 1, sampling_space,
            init_plan, init_seats, rng_states, pool, verbosity);
        // now get for proposal plan
        PlanEnsemble proposal_plan_ensemble = get_plan_ensemble(
            map_params, *splitting_schedule_ptr, initial_num_regions, 1, sampling_space,
            init_plan, init_seats, rng_states, pool, verbosity);

        // splitter
        std::vector<std::unique_ptr<TreeSplitter>> tree_splitter_ptr_vec =
            get_tree_splitter_ptrs(map_params, splitting_method, sampling_space,
                control, nsims, 1);
        if constexpr (DEBUG_PURE_MS_VERBOSE)
            Rprintf("Checkpoint 4!\n");
        // sanity check make sure the plan is ok
        std::vector<bool> county_component_lookup(ndists * map_params.num_counties, false);
        bool hierarchically_valid = current_plan_multigraph.is_hierarchically_valid(
            *plan_ensemble.plan_ptr_vec[0], county_component_lookup);

        if (!hierarchically_valid) {
            plan_ensemble.plan_ptr_vec[0]->Rprint(true);
            throw Rcpp::exception("Initial Plan is not hierarchically valid. Either turn off "
                                  "counties or pass in a hierarchically valid plan\n");
        }
        // Build the schedule context for the current plan. It is swapped with the
        // proposed plan's context whenever a proposal is accepted, so it always
        // matches the current plan.
        ScheduleContext current_schedule_ctx, proposed_schedule_ctx;
        splitting_schedule_ptr->build_context(*plan_ensemble.plan_ptr_vec[0],
                                              current_schedule_ctx);

        // build multigraph on current plan and get pairs of adj districts to merge
        auto build_result =
            plan_ensemble.plan_ptr_vec[0]->attempt_to_get_valid_mergesplit_pairs(
                current_plan_multigraph, *splitting_schedule_ptr, current_schedule_ctx,
                scoring_function, true);
        // shouldn't be possible but just a sanity check
        if (!build_result.first) {
            throw Rcpp::exception("BIG ERROR: Plan registered as hierarchically valid but we "
                                  "failed to build hierarchical multigraph!\n");
        }
        auto current_plan_adj_region_pairs =
            plan_ensemble.plan_ptr_vec[0]
                ->attempt_to_get_valid_mergesplit_pairs(
                    current_plan_multigraph, *splitting_schedule_ptr, current_schedule_ctx,
                    scoring_function, true)
                .second;

        // get weights
        std::vector<double> current_plan_pair_unnoramalized_wgts = get_adj_pair_unnormalized_weights(
            *plan_ensemble.plan_ptr_vec[0], current_plan_adj_region_pairs, pair_rule);

        // Set or estimate k if doing graph space sampling
        if (sampling_space == SamplingSpace::GraphSpace) {
            int cut_k;
            bool try_to_estimate_cut_k = Rcpp::as<bool>(control["estimate_cut_k"]);
            if (try_to_estimate_cut_k) {
                double thresh = (double)control["adapt_k_thresh"];
                double tol = std::max(target - lower, upper - target) / target;

                cut_k = estimate_mergesplit_cut_k(map_params,
                    *plan_ensemble.plan_ptr_vec[0], current_plan_multigraph,
                    *splitting_schedule_ptr, thresh, tol, rng_state);

                if (verbosity >= 3) {
                    Rcpp::Rcout << " Using estimated k = " << cut_k << "\n";
                }
            } else {
                cut_k = Rcpp::as<int>(control["manual_k"]);
                if (verbosity >= 3) {
                    Rcpp::Rcout << "Using k = " << cut_k << "\n";
                }
            }
            // update the tree splitter
            tree_splitter_ptr_vec[0]->update_single_int_param(cut_k);
            out["cut_k"] = cut_k;
        }

        if constexpr (DEBUG_PURE_MS_VERBOSE)
            Rprintf("Checkpoint 5!\n");

        if constexpr (DEBUG_PURE_MS_VERBOSE) {
            Rprintf("Adj regions are:");
            for (auto const &a_pair : current_plan_adj_region_pairs) {
                Rprintf("(%d, %d), ", a_pair.first, a_pair.second);
            }
            Rprintf("\n");
        }

        if constexpr (DEBUG_PURE_MS_VERBOSE)
            Rprintf("Checkpoint 6!\n");
        // Loading Info
        if (verbosity >= 1) {
            Rcpp::Rcout.imbue(std::locale(""));
            Rcpp::Rcout << std::fixed << std::setprecision(0);
            Rcpp::Rcout << "MERGE SPLIT MONTE CARLO" << std::endl;
            Rcpp::Rcout << "Using " << sampling_space_to_str(sampling_space);
            Rcpp::Rcout << " Sampling space to sample " << nsims << " " << V << "-unit ";
            Rcpp::Rcout << "maps with " << ndists << " districts and population between " << lower
                  << " and " << upper;
            Rcpp::Rcout << " Using " << splitting_method_to_str(splitting_method) << "!" << std::endl;
            if (map_params.cg.size() > 1) {
                Rcpp::Rcout << "Sampling hierarchically with respect to the " << map_params.cg.size()
                      << " administrative units." << std::endl;
            }
            if (scoring_function.total_soft_constraints > 0) {
                Rcpp::Rcout << "Applying " << scoring_function.total_soft_constraints
                      << " soft constraints." << std::endl;
            }
            int const num_hard_constraints = scoring_function.num_hard_region_constraints +
                                             scoring_function.num_hard_plan_constraints;
            if (num_hard_constraints > 0) {
                Rcpp::Rcout << "Applying " << num_hard_constraints << " hard constraints."
                      << std::endl;
            }
        }

        // Declare the timer for the nonwarmup part 
        std::chrono::steady_clock::time_point nonwarmup_start_time;  // no call to now(), just declares variable

        // start timing the warmup 
        auto warmup_start_time = std::chrono::steady_clock::now();

        Rcpp::RObject bar = cli_progress_bar(total_steps, cli_config(false));
        for (std::int64_t i = start, step_num = 0; i <= total_post_warmup_steps; i++, step_num++) {
            if (i == 1){
                // stop timing warmup and switch to timing runs 
                // end timing 
                auto warmup_end_time = std::chrono::steady_clock::now();
                // compute the time taken
                std::chrono::duration<double, std::ratio<1>> warmup_diff = warmup_end_time - warmup_start_time;
                warmup_time = warmup_diff.count();
                // Now start timing the rest of the code 
                nonwarmup_start_time = std::chrono::steady_clock::now();
            }
            // Index 0 or less is warmup
            bool in_warmup = i <= 0;
            if constexpr (DEBUG_PURE_MS_VERBOSE) {
                Rprintf("Iter %lld and idx %d: Starting\n",
                        static_cast<long long>(i), current_plan_mat_col);
            }

            // attempt to mergesplit
            std::tuple<bool, bool, double, int> mergesplit_result = attempt_mergesplit_step(
                map_params, *splitting_schedule_ptr, scoring_function, rng_state,
                sampling_space, *plan_ensemble.plan_ptr_vec[0],
                *proposal_plan_ensemble.plan_ptr_vec[0], ust_sampler, *tree_splitter_ptr_vec[0],
                current_plan_multigraph, proposed_plan_multigraph, current_schedule_ctx,
                proposed_schedule_ctx, pair_rule,
                save_edge_selection_prob, current_plan_adj_region_pairs,
                current_plan_pair_unnoramalized_wgts, rho, true, do_mh, false, nullptr,
                granular_times);

            if constexpr (DEBUG_PURE_MS_VERBOSE) {
                Rprintf("Iter %lld and idx %d: Attempted Mergesplit \n",
                        static_cast<long long>(i), current_plan_mat_col);
            }
            // count size
            ++tree_sizes[std::get<3>(mergesplit_result) - 1];

            // copy if needed
            if (!in_warmup && i % thin == 0) {
                // Copy the plan into the matrix
                // since a Plan's region IDs are associated with
                std::copy(plan_ensemble.flattened_all_plans.begin(),
                          plan_ensemble.flattened_all_plans.end(),
                          saved_plans_mat.column(current_plan_mat_col)
                              .begin() // Start of column in Rcpp::IntegerMatrix
                );
                // copy region populations
                std::copy(plan_ensemble.flattened_all_region_pops.begin(),
                          plan_ensemble.flattened_all_region_pops.end(),
                          saved_district_pops_mat.column(current_plan_mat_col)
                              .begin() // Start of column in Rcpp::IntegerMatrix
                );
                // if mmd copy sizes
                if (mmd_plans) {
                    std::copy(plan_ensemble.flattened_all_region_sizes.begin(),
                              plan_ensemble.flattened_all_region_sizes.end(),
                              saved_plan_sizes.column(current_plan_mat_col)
                                  .begin() // Start of column in Rcpp::IntegerMatrix
                    );
                }
                if (diagnostic_mode) {
                    std::copy(proposal_plan_ensemble.plan_ptr_vec[0]->region_ids.begin(),
                              proposal_plan_ensemble.plan_ptr_vec[0]->region_ids.end(),
                              proposed_plans_mat.column(current_plan_mat_col)
                                  .begin() // Start of column in Rcpp::IntegerMatrix
                    );
                }

                // save ratio
                log_mh_ratios(current_plan_mat_col) = std::get<2>(mergesplit_result);
                mh_decisions(current_plan_mat_col) = std::get<1>(mergesplit_result);

                if constexpr (DEBUG_PURE_MS_VERBOSE)
                    Rprintf("log_mh_ratio = %f\n", std::get<2>(mergesplit_result));

                // increase the count of columne we're on
                current_plan_mat_col++;
            }

            // if successful then update acceptance count
            if (std::get<1>(mergesplit_result)) {
                ++successful_tree_sizes[std::get<3>(mergesplit_result) - 1];
                if (in_warmup) {
                    ++warmup_acceptances;
                } else {
                    ++post_warump_acceptances;
                }
            }

            if (verbosity >= 1 && CLI_SHOULD_TICK) {
                
                mha = static_cast<double>(warmup_acceptances + post_warump_acceptances) /
                      (step_num + 1);
                cli_progress_set_format(bar,
                                        "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} "
                                        "| MH Acceptance: %.2f",
                                        mha);
                cli_progress_set(bar, step_num + 1);
            }
            Rcpp::checkUserInterrupt();
        }
        cli_progress_done(bar);
        if constexpr (DEBUG_PURE_MS_VERBOSE)
            REprintf("Done with main MCMC loop\n");

        auto nonwarmup_end_time = std::chrono::steady_clock::now();
        // compute the time taken
        std::chrono::duration<double, std::ratio<1>> nonwarmup_diff = nonwarmup_end_time - nonwarmup_start_time;
        nonwarmup_time = nonwarmup_diff.count();
    }

    if (verbosity >= 1) {
        Rcpp::Rcout << "Acceptance rate: " << std::setprecision(2)
              << (100.0 * (warmup_acceptances + post_warump_acceptances)) / (total_steps) << "%"
              << std::endl;
    }

    /*
     * How many of the saved plans are distinct, ignoring how regions happen to
     * be numbered. Counted while the matrix is still 0-indexed, so region ids
     * run 0..ndists-1.
     *
     * Only the count is kept: the tally behind it is one entry per distinct
     * plan and is discarded here, so this adds a fraction of a percent to what
     * the saved plans matrix already occupies.
     */
    int const nunique_saved_plans = static_cast<int>(
        tally_unique_plans(
            [&](int const i) { return saved_plans_mat.column(i); },
            saved_plans_mat.ncol(), V, ndists
        ).size()
    );

    // now add 1 to plans
    std::transform(saved_plans_mat.begin(), saved_plans_mat.end(), saved_plans_mat.begin(),
                   [](int x) { return x + 1; });

    if constexpr (DEBUG_PURE_MS_VERBOSE)
        REprintf("Added one to plans, now creating diagnostic list.\n");

    out["plans"] = saved_plans_mat;
    out["region_pops"] = saved_district_pops_mat;
    out["seats"] = saved_plan_sizes;
    out["mhdecisions"] = mh_decisions;
    out["total_steps"] = total_steps;
    out["warmup_acceptances"] = warmup_acceptances;
    out["post_warump_acceptances"] = post_warump_acceptances;
    out["log_mh_ratio"] = log_mh_ratios;
    out["nunique_plans"] = nunique_saved_plans;
    out["tree_sizes"] = tree_sizes;
    out["successful_tree_sizes"] = successful_tree_sizes;
    out["warmup_time"] = warmup_time;
    out["nonwarmup_time"] = nonwarmup_time;

    if (diagnostic_mode) {
        // now add 1 to plans
        std::transform(proposed_plans_mat.begin(), proposed_plans_mat.end(),
                       proposed_plans_mat.begin(), [](int x) { return x + 1; });
        out["proposed_plans"] = proposed_plans_mat;
    }


    // add granular time info 
    Rcpp::List granular_timing;
    if constexpr (perf_config::track_granular_times){
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = true,
            Rcpp::_["wilson_call_times"] = granular_times.wilson_time,
            Rcpp::_["wilson_backfill_call_times"] = granular_times.backfill_wilson_time,
            Rcpp::_["selecting_merge_pair"] = granular_times.selecting_merge_pair,
            Rcpp::_["hard_constraint_split_times"] = granular_times.hard_constraint_time,
            Rcpp::_["eff_boundary_times"] = granular_times.eff_boundary_length,
            Rcpp::_["computing_plan_scores_times"] = granular_times.plan_scores,
            Rcpp::_["computing_region_scores_times"] = granular_times.region_scores,
            Rcpp::_["get_valid_mergepairs_times"] = granular_times.get_valid_mergepairs,
            Rcpp::_["computing_spanning_tree_count_times"] = granular_times.tau_terms,   
            Rcpp::_["plan_updating_times"] = granular_times.plan_copying
        );

    }else{
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = false
        );
    }

    out["granular_times"] = granular_timing;

    if constexpr (DEBUG_PURE_MS_VERBOSE)
        REprintf("Done now returning!\n");

    return out;
}

