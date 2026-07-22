/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2024/10
 * Purpose: Run SMC with MCMC merge-split steps mixed in
 ********************************************************/

constexpr bool DEBUG_GSMC_PLANS_VERBOSE = false; // Compile-time constant


#include <atomic>
#include <chrono>

#include <cmath>
#include <functional>
#include <iostream>
#include <string>

#include "redist_alg_helpers.h"
#include "splitting_schedule_types.h"
#include "graph_plan_type.h"
#include "merging.h"
#include "scoring.h"
#include "wilson.h"
#include "weight_caching.h"
#include "weights.h"
#include "threading_helpers.h"
#include "tree_splitting.h"
#include "utils.h"
#include "random.h"
#include <RcppThread.h>
#include <cli/progress.h>

/*
 *  Use SMC Sampler method to split a multidistrict in all of the plans
 *
 *  Using the procedure outlined in <PAPER HERE> this function attempts to split
 *  a multidistrict in a previous step's plan until `M` successful splits have been made.
 *  This is based on the `split_maps` function in smc.cpp
 *
 *  @param map_params Map parameters (adj graph, population, etc.)
 *  @param old_plans_ptr_vec A vector of smart pointers (unique) to plans from the
 *  previous step
 *  @param new_plans_ptr_vec A vector of smart pointers (unique) to plans which
 *  will be filled with plans that had a multidistrict split to make them
 *  @param tree_splitters_ptr_vec A vector of smart pointers (unique) to TreeSplitter
 *  objects which will be used to split a plan from `old_plans_ptr_vec` to make a
 *  plan in `new_plans_ptr_vec`.
 *  @param parent_index_vec A vector used to track the index of the previous plan
 *  sampled that was successfully split. The value of `parent_index_vec[i]` is the
 *  index of the old plan from which the new plan `new_plans_ptr_vec[i]` was
 *  successfully split from. In other words `new_plans_ptr_vec[i]` is equal to
 *  `attempt_region_split(old_plans_ptr_vec[parent_index_vec[i]], ...)`
 *  @param normalized_cumulative_weights A vector of weights used to sample indices
 *  of the `old_plans_ptr_vec`. The value of `normalized_cumulative_weights[i]` is
 *  the normalized cumulative probability of the weights up to index i
 *  i.e. the probability that index `i` is selected is
 *  normalized_cumulative_weights[i]-normalized_cumulative_weights[i-1]
 *  @param draw_tries_vec A vector used to keep track of how many plan split
 *  attempts were made for index i. The value `draw_tries_vec[i]` represents how
 *  many split attempts were made for the i-th new plan (including the successful
 *  split). For example, `draw_tries_vec[i] = 1` means that the first split
 *  attempt was successful.
 *  @param parent_unsuccessful_tries_vec A vector used to keep track of how many times the
 *  previous rounds plans were sampled and unsuccessfully split. The value
 *  `parent_unsuccessful_tries_vec[i]` represents how many times `old_plans_ptr_vec[i]` was
 * sampled and then unsuccessfully split while creating all `M` of the new plans. THIS MAY NOT
 * BE THREAD SAFE
 *  @param accept_rate The number of accepted splits over the total number of
 *  attempted splits. This is equal to `sum(draw_tries_vec)/M`
 *  @param n_unique_parent_indices The number of unique parent indices, ie the
 *  number of previous plans that had at least one descendant amongst the new
 *  plans. This is equal to `unique(parent_index_vec)`
 *  @param ancestors Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
 *  WHAT IT IS DOING
 *  @param lags Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
 *  WHAT IT IS DOING
 *  @param split_district_only Whether or not to only allow for single district
 *  splits. If set to `true` will only attempt to split off one district at a
 *  time
 *  @param pool A threadpool for multithreading
 *  @param verbosity A parameter controlling the amount of detail printed out
 *  during the algorithms running
 *  @param diagnostic_level What type of diagnostics to track. Not implemented
 *  yet.
 *
 *  @details Modifications
 *     - The `new_plans_ptr_vec` is updated with all the newly split plans
 *     - The `old_plans_ptr_vec` is updated with all the newly split plans as well.
 *     Note that the reason both this and `new_plans_ptr_vec` are updated is because
 *     of the nature of the code you need both vectors and so both are passed by
 *     reference to save memory.
 *     - The `original_ancestor_vec` is updated to contain the indices of the
 *     original ancestors of the new plans
 *     - The `parent_index_vec` is updated to contain the indices of the parents of the
 *     new plans
 *     - If two new valid regions are split then the new_region_ids is updated so the
 *     first entry is the first new region and the second entry is the second new region
 *     - The `draw_tries_vec` is updated to contain the number of tries for each
 *     of the new plans
 *     - The `parent_unsuccessful_tries_vec` is updated to contain the number of unsuccessful
 *     samples of the old plans
 *     - The `accept_rate` is updated to contain the average acceptance rate for
 *     this iteration
 *     - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
 *     with the unique number of parents and original ancestors for all the new
 *     plans respectively
 *     - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
 *     I DO NOT KNOW WHAT IT MEANS
 *
 */
void run_smc_step(const MapParams &map_params, SplittingSchedule const &splitting_schedule,
                  std::vector<ScoringFunction> const &scoring_functions,
                  std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
                  std::unique_ptr<PlanEnsemble> &old_plan_ensemble,
                  std::unique_ptr<PlanEnsemble> &new_plan_ensemble,
                  std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters,
                  const arma::vec &normalized_cumulative_weights,
                  SMCDiagnostics &smc_diagnostics, int const smc_step_num, int const step_num,
                  bool const is_final_split, arma::umat &ancestors, const std::vector<int> &lags,
                  bool const estimated_unbiased_normalizing_constant,
                  RcppThread::ThreadPool &pool, int verbosity, int diagnostic_level,
                  int const max_split_tries) {
    // important constants
    int const num_threads = get_num_threads(pool);
    const int M = old_plan_ensemble->nsims;
    bool const smd_split_district_only =
        splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;

    // PREVIOUS SMC CODE I DONT KNOW WHAT IT DOES
    const int dist_ctr = old_plan_ensemble->plan_ptr_vec.at(0)->num_regions;
    const int n_lags = lags.size();
    arma::umat ancestors_new(M, n_lags); // lags/ancestor thing

    // Because of multithreading we have to add specific checks for if the user
    // wants to quit the program
    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50;         // check for interrupts every _ iterations
    const bool check_max_split_tries = max_split_tries > 0; // only check if greater than 0

    // The new region in the split plans is the number of regions in a split plan minus
    // one so the number of regions in a presplit plan
    int new_region_id = old_plan_ensemble->plan_ptr_vec.at(0)->num_regions;

    // we only save for linking edge
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;

    // temporary buffers to hold info 
    std::vector<int> parent_index_buffer(M, 0);
    std::vector<int> draw_tries_buffer(M, 0);
    std::vector<int> parent_unsuccessful_by_thread(num_threads * M, 0);
    // std::vector<std::atomic<int>> parent_unsuccessful_tries_atomic(M);

    // for (int i = 0; i < M; ++i) {
    //     parent_unsuccessful_tries_atomic[i].store(0, std::memory_order_relaxed);
    // }


    // count the sizes we draw trees on
    std::vector<std::vector<int>> thread_tree_sizes(
        num_threads, std::vector<int>(map_params.total_seats, 0));
    // count the sizes of regions successful trees drawn on
    std::vector<std::vector<int>> thread_successful_tree_sizes(
        num_threads, std::vector<int>(map_params.total_seats, 0));

    // granular time tracking 
    std::vector<double> wilson_call_times (
        perf_config::track_granular_times ? num_threads : 0
    );
    std::vector<double> md_selection_times (
        perf_config::track_granular_times ? num_threads : 0
    );
    std::vector<double> plan_updating_times (
        perf_config::track_granular_times ? num_threads : 0
    );
    std::vector<double> hard_constraint_split_times (
        perf_config::track_granular_times ? num_threads : 0
    );
    

    
    // Trick to give each thread a unique id
    // We need the extra steps to avoid the problem where
    // only some threads persist from previous calls but the counter resets
    // resulting in multiple threads with the same id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};
    // optional extra checking 
    std::vector<std::atomic<int>> active_users(
        perf_config::check_threadpool_integrity ? num_threads : 0);

    if constexpr (perf_config::check_threadpool_integrity){
        for (int t = 0; t < num_threads; ++t) {
            active_users[t].store(0, std::memory_order_relaxed);
        }
    }

    // now make the vectors of important variables to be used by threads
    std::vector<USTSampler> ust_samplers_vec;
    ust_samplers_vec.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++) {
        ust_samplers_vec.emplace_back(map_params, splitting_schedule);
    }

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("About to start SMC Step for %d plans\n", M);
    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id = -1;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }

        if (thread_id < 0 || thread_id >= num_threads) {
            std::ostringstream oss;
            oss << "In `run_smc_step` Thread id broke, thread id is " << thread_id
                              << " but num threads is  " << num_threads << std::endl;
            throw std::runtime_error(oss.str());
        }
        // UNCOMMENT FOR THREADPOOL CHECKING
        // std::unique_ptr<ActiveUserGuard> active_guard;
        // if constexpr (perf_config::check_threadpool_integrity) {
        //     active_guard = std::make_unique<ActiveUserGuard>(active_users[thread_id]);
        // }

        // optional time tracker for granular time tracking 
        // no call to now() if constexpr is false, since it just declares the variable it should be optimized out
        auto total_plan_start_time = maybe_now(); // optional timing 

        bool ok = false;
        int idx;
        int reject_ct = 0;
        RcppThread::checkUserInterrupt(i % check_int == 0);

        while (!ok) {
            if (check_max_split_tries && reject_ct >= max_split_tries) {
                throw Rcpp::exception(
                    "Failed to split a single plan after `max_split_tries` attempts!\n");
            }
            // increase the number of tries for particle i by 1
            draw_tries_buffer[i]++;
            // sample previous plan
            idx = rng_states[thread_id].r_int_wgt(normalized_cumulative_weights);

            // Get region id the split
            int region_id_to_split;
            auto md_selection_time = maybe_now();
            if (smd_split_district_only) {
                // if just doing district splits just use remainder region
                // which is always the highest id
                region_id_to_split = old_plan_ensemble->plan_ptr_vec[idx]->num_regions - 1;
            } else {
                // if generalized split pick a region to try to split
                region_id_to_split =
                    old_plan_ensemble->plan_ptr_vec[idx]->choose_multidistrict_to_split(
                        splitting_schedule.valid_region_sizes_to_split, rng_states[thread_id]);
            }
            if constexpr (perf_config::track_granular_times){
                add_elapsed(md_selection_times[thread_id], md_selection_time); // optional timing
            }
            

            int region_to_split_size =
                old_plan_ensemble->plan_ptr_vec[idx]->region_sizes[region_id_to_split];
            if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                REprintf("Sampled idx %d - Region %d, Size %d!\n", idx, region_id_to_split,
                         region_to_split_size);
            }

            if constexpr(perf_config::supposedly_safe_input_checks){
                // count tree size
                std::ostringstream oss;
                oss << "Calling on thread_tree_sizes, thread id " << thread_id <<" in `run smc step`, ";
                oss << "smc index i=" << i << "\n";
                oss << old_plan_ensemble->plan_ptr_vec[idx]->debug_string(true);
                tree_size_check(
                    map_params, 
                    region_to_split_size - 1, 
                    thread_tree_sizes[thread_id],
                    oss.str()
                );
            }
            // increase the count
            ++thread_tree_sizes[thread_id][region_to_split_size - 1];

            // Try to split the region
            auto wilson_time = maybe_now();
            std::pair<bool, EdgeCut> edge_search_result =
                ust_samplers_vec[thread_id].attempt_to_find_valid_tree_split(
                    rng_states[thread_id], scoring_functions[thread_id],
                    *tree_splitters[thread_id], *old_plan_ensemble->plan_ptr_vec[idx],
                    region_id_to_split, new_region_id, save_edge_selection_prob);
            if constexpr (perf_config::track_granular_times){
                add_elapsed(wilson_call_times[thread_id], wilson_time); // optional timing
            }


            if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                REprintf("idx %d - Splti %s\n", idx,
                         (std::get<0>(edge_search_result) ? "SUCCESS" : "FAILURE"));
            }

            // if successful update the new plan and check if satisfies any other hard
            // constraints
            if (std::get<0>(edge_search_result)) {
                if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                    Rprintf("Tree on Plan %d Successfully split\n", i);
                }
                // make the new plan a copy of the old one
                auto plan_updating_time = maybe_now();
                new_plan_ensemble->plan_ptr_vec[i]->shallow_copy(
                    *old_plan_ensemble->plan_ptr_vec[idx]);
                // now split that region we found on the old one
                new_plan_ensemble->plan_ptr_vec[i]->update_from_successful_split(
                    *tree_splitters[thread_id], ust_samplers_vec[thread_id],
                    std::get<1>(edge_search_result), region_id_to_split, new_region_id, true);
                if constexpr (perf_config::track_granular_times){
                    add_elapsed(plan_updating_times[thread_id], plan_updating_time); // optional timing 
                }


                // check if there are any additional hard constraints
                if (!scoring_functions[thread_id].any_hard_constraints) {
                    ok = true;
                } else {
                    // If custom hard constraints are used then
                    // the thread pool can only have a single thread or else everything will
                    // break
                    // Now we check if the new plan satisfies all hard constraints 
                    auto split_hard_constraint_time = maybe_now();
                    ok = scoring_functions[thread_id].satisfies_hard_constraints(
                        *new_plan_ensemble->plan_ptr_vec[i], region_id_to_split, new_region_id); 
                    if constexpr (perf_config::track_granular_times){
                        add_elapsed(hard_constraint_split_times[thread_id], split_hard_constraint_time); // optional timing
                    }
                    if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                        Rprintf("Plan %d - New split has %s probability\n", i,
                                (ok ? "POSITIVE" : "ZERO"));
                    }
                }
            }

            if (ok) {
                if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                    Rprintf("Success, updating Plan %d\n", i);
                }
                RcppThread::checkUserInterrupt(i % check_int == 0);
                // means idx was ok
                // record index of new plan's parent
                parent_index_buffer[i] = idx;
                if constexpr(perf_config::supposedly_safe_input_checks){
                    // count tree size
                    std::ostringstream oss;
                    oss << "Calling on thread_successful_tree_sizes, thread id " << thread_id <<" in `run smc step`, ";
                    oss << "smc index i=" << i << "\n";
                    oss << old_plan_ensemble->plan_ptr_vec[idx]->debug_string(true);
                    tree_size_check(
                        map_params, 
                        region_to_split_size - 1, 
                        thread_successful_tree_sizes[thread_id],
                        oss.str()
                    );
                }
                // add as successful tree size
                ++thread_successful_tree_sizes[thread_id][region_to_split_size - 1];
            } else { // else bad sample so try again
                     // check for user interrupt
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);
                // if diagnostic level 2 or higher get unsuccessful count
                if (diagnostic_level >= 0) {
                    // not atomic so technically not thread safe but doesn't seem to differ in
                    // practice
                    // parent_unsuccessful_tries_atomic[idx].fetch_add(
                    //     1,
                    //     std::memory_order_relaxed
                    // );
                    ++parent_unsuccessful_by_thread[thread_id * M + idx];
                }
            }
        }

        // ORIGINAL SMC CODE I DONT KNOW WHAT THIS DOES
        // save ancestors/lags
        for (int j = 0; j < n_lags; j++) {
            if (dist_ctr <= lags[j]) {
                ancestors_new(i, j) = i;
            } else {
                ancestors_new(i, j) = ancestors(idx, j);
            }
        }
        if constexpr (perf_config::track_granular_times){
            // set the time spent successfully sampling a plan
            add_elapsed(
                smc_diagnostics.total_plan_smc_split_times(i, smc_step_num),
                total_plan_start_time
            ); // optional timing
        }

        if (verbosity >= 3) {
            ++bar;
        }
        
    });



    // Wait for all the threads to finish
    pool.wait();

    Rcpp::checkUserInterrupt();

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
        REprintf("Done splitting!\n");
    }

    if(estimated_unbiased_normalizing_constant){
        if constexpr(DEBUG_GSMC_PLANS_VERBOSE){
            REprintf("Starting estimation: ");
        }

    // Now we sample one more plan and discard it to allow for unbiased 
    // normalization constant estimation 
    bool extra_plan_sampled = false;
    int extra_particle_reject_ct = 0;

    while(!extra_plan_sampled) {
        if constexpr(DEBUG_GSMC_PLANS_VERBOSE) REprintf("%d ", extra_particle_reject_ct);
        if (check_max_split_tries && extra_particle_reject_ct >= max_split_tries) {
            throw Rcpp::exception(
                "Failed to split a single plan after `max_split_tries` attempts!\n");
        }
        // increase the number of tries for particle i by 1
        extra_particle_reject_ct++;
        // sample previous plan
        int idx = rng_states[0].r_int_wgt(normalized_cumulative_weights);
        // Get region id the split
        int region_id_to_split;
        if (smd_split_district_only) {
            // if just doing district splits just use remainder region
            // which is always the highest id
            region_id_to_split = old_plan_ensemble->plan_ptr_vec[idx]->num_regions - 1;
        } else {
            // if generalized split pick a region to try to split
            region_id_to_split =
                old_plan_ensemble->plan_ptr_vec[idx]->choose_multidistrict_to_split(
                    splitting_schedule.valid_region_sizes_to_split, rng_states[0]);
        }

        // Try to split the region
        std::pair<bool, EdgeCut> edge_search_result =
            ust_samplers_vec[0].attempt_to_find_valid_tree_split(
                rng_states[0], scoring_functions[0],
                *tree_splitters[0], *old_plan_ensemble->plan_ptr_vec[idx],
                region_id_to_split, new_region_id, save_edge_selection_prob);



        // if successful update the new plan and check if satisfies any other hard
        // constraints
        if (std::get<0>(edge_search_result)) {
            // check if there are any additional hard constraints
            if (!scoring_functions[0].any_hard_constraints) {
                // if not we can stop trying 
                extra_plan_sampled = true;
            } else {
                // If custom hard constraints are used then
                // the thread pool can only have a single thread or else everything will
                // break
                // TODO: Make it possible to check this new plan without actually copying anything
                // since this is an extra plan we discard
                extra_plan_sampled = true;
                // make the new plan a copy of the old one
                // new_plan_ensemble->plan_ptr_vec[i]->shallow_copy(
                //     *old_plan_ensemble->plan_ptr_vec[idx]);
                // // now split that region we found on the old one
                // new_plan_ensemble->plan_ptr_vec[i]->update_from_successful_split(
                //     *tree_splitters[thread_id], ust_samplers_vec[thread_id],
                //     std::get<1>(edge_search_result), region_id_to_split, new_region_id, true);

                // ok = scoring_functions[thread_id].satisfies_hard_constraints(
                //     *new_plan_ensemble->plan_ptr_vec[i], region_id_to_split, new_region_id,
                //     is_final_split);
            }
        }
    }
    if constexpr(DEBUG_GSMC_PLANS_VERBOSE){
        REprintf("Done!\n");
    }
    

    // now save the number of failed attempts before sampling the extra plan
    smc_diagnostics.tries_before_extra_particle[smc_step_num] = extra_particle_reject_ct;

    }

    // now swap the old plans with the new ones. This avoids needing to actually copy
    std::swap(old_plan_ensemble, new_plan_ensemble);

    // Add the buffer info in
    for (int i = 0; i < M; ++i) {
        smc_diagnostics.parent_index_mat(i, smc_step_num) = parent_index_buffer[i]; 
        smc_diagnostics.draw_tries_mat(i, step_num) = draw_tries_buffer[i];

        // smc_diagnostics.parent_unsuccessful_tries_mat(i, smc_step_num) =
        //     parent_unsuccessful_tries_atomic[i].load(std::memory_order_relaxed);
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
        {
            smc_diagnostics.parent_unsuccessful_tries_mat(i, smc_step_num) +=
            parent_unsuccessful_by_thread[thread_id * M + i];
        }
    }


    // update tree sizes counts
    for (size_t region_size = 0; region_size < map_params.total_seats; region_size++) {
        for (size_t a_thread_id = 0; a_thread_id < num_threads; a_thread_id++) {
            smc_diagnostics.tree_sizes_mat(region_size, step_num) +=
                thread_tree_sizes[a_thread_id][region_size];
            smc_diagnostics.successful_tree_sizes_mat(region_size, step_num) +=
                thread_successful_tree_sizes[a_thread_id][region_size];
        }
    }

    // add the wilson call times 
    if constexpr (perf_config::track_granular_times){
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
        {
            smc_diagnostics.wilson_call_times[step_num] += wilson_call_times[thread_id];
            smc_diagnostics.md_selection_times[smc_step_num] += md_selection_times[thread_id];
            smc_diagnostics.hard_constraint_split_times[step_num] += hard_constraint_split_times[thread_id];
            smc_diagnostics.plan_updating_times[step_num] += plan_updating_times[thread_id];
        }
    }

    

    // now compute acceptance rate and unique parents and original ancestors
    double accept_rate = M / static_cast<double>(std::accumulate(draw_tries_buffer.begin(), draw_tries_buffer.end(), 0));
    smc_diagnostics.acceptance_rates.at(step_num) = accept_rate;

    // Get number of unique parents
    std::set<int> unique_parents(parent_index_buffer.begin(), parent_index_buffer.end());
    smc_diagnostics.nunique_parents.at(smc_step_num) = unique_parents.size();
    // Get the number of unique plans
    smc_diagnostics.nunique_plans[step_num] = old_plan_ensemble->count_unique_plans(pool);

    if (verbosity >= 3) {
        Rcpp::Rcout << "  " << std::setprecision(2) << 100.0 * accept_rate << "% acceptance rate. "
              << 100.0 * smc_diagnostics.nunique_parents.at(smc_step_num) / M
              << "% of previous step's plans survived," << " and there are now "
              << smc_diagnostics.nunique_plans[step_num] << " unique plans." << std::endl;
    }

    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;
}

void run_merge_split_step_on_all_plans(
    RcppThread::ThreadPool &pool, MapParams const &map_params,
    const SplittingSchedule &splitting_schedule,
    std::vector<ScoringFunction> const &scoring_functions, std::vector<RNGState> &rng_states,
    SamplingSpace const sampling_space, std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec,
    std::vector<std::unique_ptr<TreeSplitter>> &tree_splitter_ptrs_vec,
    std::string const merge_prob_type, double const rho, bool const is_final,
    int const nsteps_to_run, int const merge_split_step_num, int const step_num,
    SMCDiagnostics &smc_diagnostics, WeightCacheEnsemble &cache_ensemble, int verbosity) {

    const int check_int = 15; // check for interrupts every _ iterations
    int nsims = (int)plan_ptrs_vec.size();
    int const num_threads = get_num_threads(pool);
    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Going to run %d steps!\n", nsteps_to_run);

    // Diagnostics
    std::vector<int> success_count_buffer(nsims, 0);
    std::vector<double> total_plan_mcmc_time_buffer(
        perf_config::track_granular_times ? nsims : 0,
        0.0
    );

    // count the sizes we draw trees on
    std::vector<std::vector<int>> thread_tree_sizes(
        num_threads, std::vector<int>(map_params.total_seats, 0));
    std::vector<std::vector<int>> thread_successful_tree_sizes(
        num_threads, std::vector<int>(map_params.total_seats, 0));

    // thread safe id counter
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};
    // debug thing
    std::vector<std::atomic<int>> active_users(
        perf_config::check_threadpool_integrity ? num_threads : 0);

    // now make the vectors of important variables to be used by threads
    std::vector<USTSampler> ust_samplers_vec;
    ust_samplers_vec.reserve(num_threads);
    std::vector<PlanMultigraph> current_plan_multigraphs_vec;
    current_plan_multigraphs_vec.reserve(num_threads);
    std::vector<PlanMultigraph> proposed_plan_multigraphs_vec;
    proposed_plan_multigraphs_vec.reserve(num_threads);

    for (size_t i = 0; i < num_threads; i++) {
        ust_samplers_vec.emplace_back(map_params, splitting_schedule);
        current_plan_multigraphs_vec.emplace_back(
            map_params, sampling_space == SamplingSpace::LinkingEdgeSpace);
        proposed_plan_multigraphs_vec.emplace_back(
            map_params, sampling_space == SamplingSpace::LinkingEdgeSpace);
    }
    std::vector<GranularMCMCTimes> granular_weight_times(num_threads);

    // create a progress bar
    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id = -1;
        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }
        if (thread_id < 0 || thread_id >= num_threads) {
            std::ostringstream oss;
            oss << "In `run_merge_split_step_on_all_plans` Thread id broke, thread id is " << thread_id
                              << " but num threads is  " << num_threads << std::endl;
            throw std::runtime_error(oss.str());
        }
        // UNCOMMENT FOR THREADPOOL CHECKING
        // std::unique_ptr<ActiveUserGuard> active_guard;
        // if constexpr (perf_config::check_threadpool_integrity) {
        //     active_guard = std::make_unique<ActiveUserGuard>(active_users[thread_id]);
        // }

        auto total_plan_start_time = maybe_now();

        // store the number of succesful runs
        if (cache_ensemble.using_caching) {
            success_count_buffer[i] = run_merge_split_steps(
                map_params, splitting_schedule, scoring_functions[thread_id],
                rng_states[thread_id], sampling_space, *plan_ptrs_vec[i], *new_plan_ptrs_vec[i],
                ust_samplers_vec[thread_id], *tree_splitter_ptrs_vec[thread_id],
                current_plan_multigraphs_vec[thread_id],
                proposed_plan_multigraphs_vec[thread_id], merge_prob_type, rho, is_final,
                nsteps_to_run, thread_tree_sizes[thread_id],
                thread_successful_tree_sizes[thread_id], true,
                cache_ensemble.weight_cache_ptr_vec[i].get(),
                granular_weight_times[thread_id]);
        } else {
            success_count_buffer[i] = run_merge_split_steps(
                map_params, splitting_schedule, scoring_functions[thread_id],
                rng_states[thread_id], sampling_space, *plan_ptrs_vec[i], *new_plan_ptrs_vec[i],
                ust_samplers_vec[thread_id], *tree_splitter_ptrs_vec[thread_id],
                current_plan_multigraphs_vec[thread_id],
                proposed_plan_multigraphs_vec[thread_id], merge_prob_type, rho, is_final,
                nsteps_to_run, thread_tree_sizes[thread_id],
                thread_successful_tree_sizes[thread_id], false, nullptr,
                granular_weight_times[thread_id]);
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);

        if constexpr (perf_config::track_granular_times){
            add_elapsed(
                total_plan_mcmc_time_buffer[i],
                total_plan_start_time
            );
        }

        if (verbosity >= 3) {
            ++bar;
        }
    });

    // Wait for all the threads to finish
    pool.wait();

    Rcpp::checkUserInterrupt();

    // update
    for (size_t i = 0; i < nsims; i++)
    {
        smc_diagnostics.merge_split_successes_mat(i, merge_split_step_num) = success_count_buffer[i];
        if constexpr (perf_config::track_granular_times){
            smc_diagnostics.total_plan_mcmc_times(i, merge_split_step_num) = total_plan_mcmc_time_buffer[i];
        }
    }
    

    // update tree sizes counts
    for (size_t region_size = 0; region_size < map_params.total_seats; region_size++) {
        for (size_t a_thread_id = 0; a_thread_id < num_threads; a_thread_id++) {
            smc_diagnostics.tree_sizes_mat(region_size, step_num) +=
                thread_tree_sizes[a_thread_id][region_size];
            smc_diagnostics.successful_tree_sizes_mat(region_size, step_num) +=
                thread_successful_tree_sizes[a_thread_id][region_size];
        }
    }

    // add granular time if that's being tracked 
    if constexpr (perf_config::track_granular_times){
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
        {
            smc_diagnostics.wilson_call_times[step_num] += granular_weight_times[thread_id].wilson_time;
            smc_diagnostics.get_valid_pairs_times[step_num] += granular_weight_times[thread_id].get_valid_pairs;
            smc_diagnostics.selecting_merge_pair_times[merge_split_step_num] += granular_weight_times[thread_id].selecting_merge_pair;
            smc_diagnostics.hard_constraint_split_times[step_num] += granular_weight_times[thread_id].hard_constraint_time;
            smc_diagnostics.eff_boundary_times[merge_split_step_num] += granular_weight_times[thread_id].eff_boundary_length;
            smc_diagnostics.plan_scores_times[step_num] += granular_weight_times[thread_id].plan_scores;
            smc_diagnostics.region_scores_times[step_num] += granular_weight_times[thread_id].region_scores;
            smc_diagnostics.log_tau_times[step_num] += granular_weight_times[thread_id].tau_terms;
            smc_diagnostics.plan_updating_times[step_num] += granular_weight_times[thread_id].plan_copying;
        }
    }

    return;
}

// Different diagnostic levels
//      - level 0 - Does not capture any ancestry information or retain intermediate weights
//      - level 1 - Saves ancestry information, intermediate weights and the number of tries
//      - level 2 - Captures the parent tries mat
//      - level 3 - Saves intermediate region dvals and plan ids

// Run SMC (optionally with Merge Split steps too)
//
// Uses smc method with optimal weights and merge split steps to generate a sample of `nsims`
// plans in `c++`
//
//
// Using the procedure outlined in <PAPER HERE> this function uses Sequential
// Monte Carlo (SMC) methods to generate a sample of `M` plans
//
//
// @param ndists The number of districts the final plans will have
// @param adj_list A 0-indexed adjacency list representing the undirected graph
// which represents the underlying map the plans are to be drawn on
// @param counties Vector of county labels of each vertex in `g`
// @param pop A vector of the population associated with each vertex in `g`
// @param target Ideal population of a valid district. This is what deviance is calculated
// relative to
// @param lower Acceptable lower bounds on a valid district's population
// @param upper Acceptable upper bounds on a valid district's population
// @param nsims The number of plans (samples) to draw
// @param k_param The k parameter from the SMC algorithm, you choose among the top k_param
// edges ' @param control Named list of additional parameters. ' @param num_threads The number
// of threads the threadpool should use ' @param verbosity What level of detail to print out
// while the algorithm is ' running <ADD OPTIONS> ' @keywords internal ' @noRd
// [[Rcpp::export]]
Rcpp::List run_redist_smc(
    int const nsims, int const total_seats, int const ndists,
    Rcpp::IntegerVector const district_seat_sizes, int const initial_num_regions,
    Rcpp::List const &adj_list, Rcpp::IntegerVector const &counties, const Rcpp::IntegerVector &pop,
    Rcpp::CharacterVector const &step_types, double const target, double const lower,
    double const upper,
    double const rho,                      // compactness
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    Rcpp::List const &control, // control has pop temper, and k parameter value, and splitting method are allowed
    Rcpp::List const &constraints, // constraints
    int const verbosity, int const diagnostic_level, Rcpp::IntegerMatrix const &region_id_mat,
    Rcpp::IntegerMatrix const &region_sizes_mat, arma::vec &log_weights) {
    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        REprintf("Inside c++ code!\n");
    bool diagnostic_mode = diagnostic_level == 1;
    // set the number of threads
    int num_threads = (int)control["num_threads"];
    if (num_threads <= 0)
        num_threads = std::thread::hardware_concurrency();


    // Set the sampling space
    SamplingSpace sampling_space = get_sampling_space(sampling_space_str);

    
    // Create map level graph and county level multigraph
    MapParams const map_params(list_to_graph(adj_list), 
        Rcpp::as<std::vector<unsigned int>>(counties), 
        Rcpp::as<std::vector<unsigned int>>(pop),
         ndists, total_seats,
                               Rcpp::as<std::vector<int>>(district_seat_sizes), lower, target, upper,
                            sampling_space);
    int V = map_params.g.size();

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
        REprintf("District Seat Sizes: ");
        for (int i = 1; i <= total_seats; i++) {
            REprintf("%d - %s,", i,
                     (map_params.is_district[i] ? "YES DISTRICT" : "NOT A DISTRICT"));
        }
        REprintf("\n");
    }

    // Add scoring functions (constraints)
    // one per thread
    std::vector<ScoringFunction> scoring_functions;
    scoring_functions.reserve(num_threads);
    for (size_t thread_id = 0; thread_id < num_threads; thread_id++) {
        scoring_functions.emplace_back(map_params, constraints,
                                       Rcpp::as<double>(control["pop_temper"]), true, thread_id);
    }

    // if we have any custom hard constraints then we have to single thread everything
    if (scoring_functions[0].any_hard_custom_constraints) {
        num_threads = 1;
    };

    // get seq_alpha
    double weights_alpha = Rcpp::as<double>(control["seq_alpha"]);
    bool const apply_weights_alpha = weights_alpha != 1;

    // re-seed MT so that `set.seed()` works in R
    int global_rng_seed = (int)Rcpp::sample(INT_MAX, 1)[0];
    int num_rng_states = num_threads;
    std::vector<RNGState> rng_states;
    rng_states.reserve(num_rng_states);
    for (size_t i = 1; i <= num_rng_states; i++) {
        // same seed with i*3 long_jumps for state
        rng_states.emplace_back(global_rng_seed, i * 3);
    }

    // Do not support hard plan constraints with linking edge
    if (sampling_space == SamplingSpace::LinkingEdgeSpace &&
        scoring_functions[0].any_hard_plan_constraints) {
        // The issue right now is for a merged plan we need to know what pairs in the merged
        // plan are valid For region based constraints merging two regions doens't affect the
        // others but theoretically for the entire plan a merge in the original plan could be ok
        // and then cease to be ok after two other regions are merged. Checks need to be added
        // for that
        throw Rcpp::exception("Whole Plan constraints with thresholding is not supported for "
                              "linking edge space sampling yet!\n");
    }

    // Legacy, in future remove
    RNGState rng_state((int)Rcpp::sample(INT_MAX, 1)[0]);
    global_seed_rng((int)Rcpp::sample(INT_MAX, 1)[0]);
    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        REprintf("RNG States created!\n");

    // unpack control params
    // lags thing (copied from original smc code, don't understand what its doing)
    std::vector<int> lags = Rcpp::as<std::vector<int>>(control["lags"]);
    arma::umat ancestors(nsims, lags.size(), arma::fill::zeros);
    // weight type
    std::string wgt_type = Rcpp::as<std::string>(control["weight_type"]);
    // whether or not to cache the weights
    bool const using_caching = Rcpp::as<bool>(control["cache_weights"]);
    // max tries value
    int const max_split_tries = Rcpp::as<int>(control["max_split_tries"]);
    // unbiased normalizing estimate
    bool const estimated_unbiased_normalizing_constant = Rcpp::as<bool>(control["est_norm_unbiased"]); 

    if(estimated_unbiased_normalizing_constant && scoring_functions[0].any_hard_constraints){
        Rcpp::warning("Unbiased normalizing constant estimation si not support right now for hard constraints!");
    }

    // total number of steps to run
    int total_steps = static_cast<int>(step_types.size());
    int total_ms_steps = 0;
    int total_smc_steps = 0;
    std::vector<bool> merge_split_step_vec(step_types.size());
    for (size_t i = 0; i < step_types.size(); i++) {
        if (static_cast<std::string>(step_types.at(i)) == "smc") {
            merge_split_step_vec.at(i) = false;
            total_smc_steps++;
        } else if (static_cast<std::string>(step_types.at(i)) == "ms") {
            merge_split_step_vec.at(i) = true;
            total_ms_steps++;
        } else {
            REprintf("Invalid step type: %s\n",
                     static_cast<std::string>(step_types.at(i)).c_str());
            throw Rcpp::exception("Invalid step type passed!");
        }
    }
    // sanity check we're not splitting more than ndists districts
    if (initial_num_regions + total_smc_steps > ndists) {
        REprintf("Trying to do %d splits with %d initial regions will "
                 "create more than ndists=%d districts!\n",
                 total_smc_steps, initial_num_regions, ndists);
        throw Rcpp::exception(
            "Desired number of splits will produce more than ndist districts!");
    }
    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        REprintf("Step types vec created!\n");

    // see if we are splitting plans all the way or just creating partial plans
    bool splitting_all_the_way = ndists == initial_num_regions + total_smc_steps;

    // multipler for number of merge split steps
    double ms_steps_multiplier;
    std::string merge_prob_type;
    if (total_ms_steps > 0) {
        ms_steps_multiplier = Rcpp::as<double>(control["mh_accept_per_smc"]);
        merge_prob_type = Rcpp::as<std::string>(control["pair_rule"]);
    }

    double tol = std::max(target - lower, upper - target) / target;

    // get splitting type
    SplittingMethodType splitting_method =
        get_splitting_type(static_cast<std::string>(control["splitting_method"]));

    // get the splitting size regime
    SplittingSizeScheduleType splitting_size_regime =
        get_splitting_size_regime(static_cast<std::string>(control["splitting_size_regime"]));
    auto splitting_schedule_ptr = get_splitting_schedule(
        total_smc_steps, ndists, total_seats, Rcpp::as<std::vector<int>>(district_seat_sizes),
        splitting_size_regime, control);
    // it wants presplit number of regions so make initial regions - 1
    // Needed for initializing linking edge plans
    splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(0, initial_num_regions -
                                                                               1);

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        REprintf("Splitting Schedule Obj created!\n");

    // Whether or not to only do district splits only
    bool const multi_member_districting =
        (splitting_size_regime == SplittingSizeScheduleType::DistrictOnlyMMD ||
         splitting_size_regime == SplittingSizeScheduleType::AnyValidSizeMMD);
    bool split_district_only =
        splitting_size_regime == SplittingSizeScheduleType::DistrictOnlySMD;

    // Do some input checking
    // Make sure first merge split argument isn't true
    if (merge_split_step_vec.at(0)) {
        throw Rcpp::exception("The first entry of merge_split_step_vec cannot be true.");
    };

    // Now create diagnostic information
    SMCDiagnostics smc_diagnostics(
        sampling_space, splitting_method, splitting_size_regime, merge_split_step_vec, V, nsims,
        ndists, total_seats, initial_num_regions, total_smc_steps, total_ms_steps,
        estimated_unbiased_normalizing_constant, 
        diagnostic_level, splitting_all_the_way, split_district_only);

    // Create a threadpool
    RcppThread::ThreadPool pool(num_threads > 1 ? num_threads : 0);

    // If hard custom then switch to no threading
    if (scoring_functions[0].any_hard_custom_constraints) {
        pool.setNumThreads(0);
    }

    // Now we add everything here to a scope since it won't be needed for the end
    // create the ensemble
    std::unique_ptr<PlanEnsemble> plan_ensemble_ptr = get_plan_ensemble_ptr(
        map_params, *splitting_schedule_ptr, initial_num_regions, nsims, sampling_space,
        region_id_mat, region_sizes_mat, rng_states, pool, verbosity);

    // compute the log of the unnormalized density of the entire map 
    // which is log spanning tree count  - score 
    double log_blank_map_target_density;

    {
        // Ensemble of dummy plans for copying
        std::unique_ptr<PlanEnsemble> dummy_plan_ensemble_ptr = get_plan_ensemble_ptr(
            map_params, *splitting_schedule_ptr, initial_num_regions, nsims, sampling_space,
            region_id_mat, region_sizes_mat, rng_states, pool, verbosity);

        // Get the tree splitter
        std::vector<std::unique_ptr<TreeSplitter>> tree_splitter_ptrs_vec =
            get_tree_splitter_ptrs(map_params, splitting_method, sampling_space,
                 control, nsims, num_threads);

        bool use_naive_k_splitter = splitting_method == SplittingMethodType::NaiveTopK;
        // adaptive k estimation threshold
        bool try_to_estimate_cut_k;
        double thresh;
        // k param values to potentially use. If set to 0 or lower then estimate
        std::vector<int> k_params;
        if (use_naive_k_splitter) {
            try_to_estimate_cut_k = Rcpp::as<bool>(control["estimate_cut_k"]);
            if (try_to_estimate_cut_k) {
                thresh = (double)control["adapt_k_thresh"];
                k_params.resize(total_smc_steps);
            } else {
                k_params = Rcpp::as<std::vector<int>>(control["manual_k_params"]);
            }
        }

        // Start off all the unnormalized weights at at exp of log weights
        arma::vec unnormalized_sampling_weights = arma::exp(log_weights);
        // now get initial normalized weights
        arma::vec normalized_cumulative_weights = arma::cumsum(unnormalized_sampling_weights);
        normalized_cumulative_weights =
            normalized_cumulative_weights / normalized_cumulative_weights[nsims - 1];

        // Create the weight cache's if needed
        std::unique_ptr<WeightCacheEnsemble> cache_ensemble_ptr =
            std::make_unique<WeightCacheEnsemble>(using_caching, map_params, nsims, rho,
                                                  sampling_space);

        std::unique_ptr<WeightCacheEnsemble> dummy_cache_ensemble_ptr =
            std::make_unique<WeightCacheEnsemble>(using_caching, map_params, nsims, rho,
                                                  sampling_space);

        double entire_map_compactness = 0.0;
        // compute the whole map compactness if we're starting from scratch
        if (initial_num_regions == 1) {
            double entire_map_log_st_count = plan_ensemble_ptr->plan_ptr_vec[0]->compute_log_region_spanning_trees(
                    map_params, 0);
            entire_map_compactness = (rho - 1) * entire_map_log_st_count;
            log_blank_map_target_density = rho * entire_map_log_st_count;

            // now subtract the score 
            auto region_score_result = scoring_functions[0].compute_region_full_score(
                *plan_ensemble_ptr->plan_ptr_vec[0], 0);
            auto plan_only_score_result =
                scoring_functions[0].compute_plan_score(*plan_ensemble_ptr->plan_ptr_vec[0]);
            log_blank_map_target_density -= region_score_result.second;
            log_blank_map_target_density -= plan_only_score_result.second;
        }

        // Loading Info
        if (verbosity >= 1) {
            Rcpp::Rcout.imbue(std::locale(""));
            Rcpp::Rcout << std::fixed << std::setprecision(0);
            if (!split_district_only) {
                Rcpp::Rcout << "GENERALIZED SEQUENTIAL MONTE CARLO";
            } else {
                Rcpp::Rcout << "SEQUENTIAL MONTE CARLO";
            }
            if (total_ms_steps > 0) {
                Rcpp::Rcout << " WITH MERGE SPLIT";
            }
            Rcpp::Rcout << std::endl;
            Rcpp::Rcout << "Using " << sampling_space_to_str(sampling_space);
            Rcpp::Rcout << " Sampling space to sample " << nsims << " " << V << "-unit ";
            Rcpp::Rcout << "maps with " << ndists << " districts and population between " << lower
                  << " and " << upper << "." << std::endl;
            if (verbosity >= 3) {
                Rcpp::Rcout << "Using " << (pool.getNumThreads() == 0 ? 1 : pool.getNumThreads())
                      << " threads, " << total_ms_steps << " merge split steps, ";
                if (splitting_size_regime == SplittingSizeScheduleType::DistrictOnlySMD ||
                    splitting_size_regime == SplittingSizeScheduleType::DistrictOnlyMMD) {
                    Rcpp::Rcout << "and only performing 1-district splits.";
                } else if (splitting_size_regime ==
                           SplittingSizeScheduleType::AnyValidSizeSMD) {
                    Rcpp::Rcout << "and generalized region splits.";
                } else if (splitting_size_regime == SplittingSizeScheduleType::OneCustomSize) {
                    Rcpp::Rcout << "and custom size region splits.";
                }
                Rcpp::Rcout << " Using " << splitting_method_to_str(splitting_method) << " with "
                      << (wgt_type == "optimal" ? "Optimal" : "Simple") << " Weights!\n";
                if (map_params.cg.size() > 1) {
                    Rcpp::Rcout << "Ensuring no more than " << ndists - 1 << " splits of the "
                          << map_params.cg.size() << " administrative units.\n";
                }
                if (scoring_functions[0].total_soft_constraints > 0) {
                    Rcpp::Rcout << "Applying " << scoring_functions[0].total_soft_constraints
                          << " soft constraints.\n";
                }
                if (scoring_functions[0].num_hard_plan_constraints > 0) {
                    Rcpp::Rcout << "Applying " << scoring_functions[0].num_hard_plan_constraints
                          << " hard constraints.\n";
                }
            }
        }

        // keep track of if we need to swap at the end.
        // counts the number of smc steps
        int smc_step_num = 0;
        int merge_split_step_num = 0;

        std::string bar_fmt =
            "Split [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} | ETA{cli::pb_eta}";
        Rcpp::RObject bar = cli_progress_bar(total_steps, cli_config(false, bar_fmt.c_str()));
        // Now for each run through split the map
        try {
            for (int step_num = 0; step_num < total_steps; step_num++) {
                if (verbosity > 1) {
                    if (merge_split_step_vec[step_num]) {
                        Rprintf("Iteration %d: Merge Split Step %d \n", step_num + 1,
                                merge_split_step_num + 1);
                    } else {
                        Rprintf("Iteration %d: SMC Step %d of %d \n", step_num + 1,
                                smc_step_num + 1, total_smc_steps);
                    }
                }
                // its the final splitting step if step_num + 1 == total_smc steps
                bool const is_final_splitting_step = step_num + 1 == total_smc_steps;
                // If we have any custom hard constraints then must switch to single threading
                // for everything

                // Check what step type
                if (!merge_split_step_vec[step_num]) {

                    // set the splitting schedule
                    splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
                        smc_step_num, plan_ensemble_ptr->plan_ptr_vec[0]->num_regions);

                    // splitting_schedule_ptr->print_current_step_splitting_info();

                    // Print if needed
                    if (verbosity >= 3 &&
                        splitting_size_regime == SplittingSizeScheduleType::OneCustomSize) {
                        splitting_schedule_ptr->print_current_step_splitting_info();
                    }

                    // check if k is passed in or estimate
                    if (use_naive_k_splitter) {
                        if (try_to_estimate_cut_k) {
                            // Start timing 
                            auto smc_param_estimation_start_time = std::chrono::steady_clock::now();
                            // est k
                            int est_cut_k;
                            int last_k = smc_step_num == 0 ? std::max(1, V - 5)
                                                           : k_params.at(smc_step_num - 1);
                            if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                                Rprintf("About to try to estimate cut k!\n");
                            estimate_cut_k(map_params, *splitting_schedule_ptr, rng_state,
                                           est_cut_k, last_k, unnormalized_sampling_weights,
                                           thresh, tol, plan_ensemble_ptr->plan_ptr_vec,
                                           split_district_only, verbosity);
                            // end timing 
                            auto smc_param_estimation_time = std::chrono::steady_clock::now();
                            // add the time 
                            std::chrono::duration<double, std::ratio<1>> smc_param_diff = smc_param_estimation_time - smc_param_estimation_start_time;
                            smc_diagnostics.smc_step_parameter_estimation_times[smc_step_num] = smc_param_diff.count();
                            if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                                Rprintf("Estimated cut k!\n");
                            k_params.at(smc_step_num) = est_cut_k;

                            if (verbosity >= 3) {
                                Rcpp::Rcout << " (using estimated k = " << k_params.at(smc_step_num)
                                      << ")" << std::endl;
                                Rcpp::Rcout << " Estimation took " 
                                      << smc_diagnostics.smc_step_parameter_estimation_times[smc_step_num]
                                      << " seconds."
                                      << std::endl;
                            }
                        } else {
                            if (verbosity >= 3) {
                                Rcpp::Rcout << " (using input k = " << k_params.at(smc_step_num)
                                      << ")\n";
                            }
                        }
                        for (auto &tree_splitter_ptr : tree_splitter_ptrs_vec) {
                            tree_splitter_ptr->update_single_int_param(
                                k_params.at(smc_step_num));
                        }

                        smc_diagnostics.cut_k_values.at(step_num) =
                            tree_splitter_ptrs_vec[0]->get_single_int_param();
                    }


                    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                        Rprintf("About to run smc step %d!\n", smc_step_num);
                    // optional integrity check 
                    if constexpr (perf_config::object_integrity_checking){
                        std::ostringstream oss;
                        oss << "Plan integrity check failed BEFORE calling `run_smc_step` ";
                        oss << "on smc_step =" << smc_step_num;
                        oss << " total step " << step_num << "\n";
                        plan_ensemble_ptr->check_all_plans_valid(
                            map_params,
                            "oss.str()"
                        );                        
                    }
                    // start timing the smc split
                    auto smc_splitting_start_time = std::chrono::steady_clock::now();
                    // split the map
                    run_smc_step(map_params, *splitting_schedule_ptr, scoring_functions,
                                 rng_states, sampling_space, plan_ensemble_ptr,
                                 dummy_plan_ensemble_ptr, tree_splitter_ptrs_vec,
                                 normalized_cumulative_weights, smc_diagnostics, smc_step_num,
                                 step_num, is_final_splitting_step, ancestors, lags, 
                                 estimated_unbiased_normalizing_constant, pool,
                                 verbosity, diagnostic_mode ? 3 : 0, max_split_tries);
                    // end timing 
                    auto smc_splitting_end_time = std::chrono::steady_clock::now();
                    // add the time 
                    std::chrono::duration<double, std::ratio<1>> smc_split_diff = smc_splitting_end_time - smc_splitting_start_time;
                    smc_diagnostics.smc_split_times[smc_step_num] = smc_split_diff.count();

                    if (verbosity >= 3) {
                        Rcpp::Rcout << "  Performing SMC splits took " 
                              << smc_diagnostics.smc_split_times[smc_step_num]
                              << " seconds."
                              << std::endl;
                    }

                    if constexpr (perf_config::object_integrity_checking){
                        std::ostringstream oss;
                        oss << "Plan integrity check failed AFTER calling `run_smc_step` ";
                        oss << "on smc_step =" << smc_step_num;
                        oss << " total step " << step_num << "\n";
                        plan_ensemble_ptr->check_all_plans_valid(
                            map_params,
                            "oss.str()"
                        );                        
                    }

                    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                        Rprintf("Ran smc step %d!\n", smc_step_num);
                    if (smc_step_num == 0 && initial_num_regions == 1) {
                        // For the first ancestor one make every ancestor themselves
                        std::iota(smc_diagnostics.parent_index_mat.column(0).begin(),
                                  smc_diagnostics.parent_index_mat.column(0).end(), 0);
                    }

                    // update the weights if caching
                    if (using_caching) {
                        // use the dummy one as a temp
                        for (size_t i = 0; i < nsims; i++) {
                            // for plan i it was split from plan
                            // smc_diagnostics.parent_index_mat(i, smc_step_num) so we want to
                            // make plan i in the dummy one to be the cache from the split plan
                            dummy_cache_ensemble_ptr->weight_cache_ptr_vec[i]->copy_from(
                                *cache_ensemble_ptr
                                     ->weight_cache_ptr_vec[smc_diagnostics.parent_index_mat(
                                         i, smc_step_num)]);
                        }
                        // now we swap the two pointers so the cache ensemble one is aligned
                        std::swap(cache_ensemble_ptr, dummy_cache_ensemble_ptr);
                    }

                    // compute splitting probability if MMD or if Anysplits SMD and num regions
                    // isn't number of districts
                    bool compute_log_splitting_prob =
                        (multi_member_districting ||
                         (splitting_schedule_ptr->schedule_type ==
                              SplittingSizeScheduleType::AnyValidSizeSMD &&
                          plan_ensemble_ptr->plan_ptr_vec[0]->num_regions != ndists));

                    // if soft custom constraints and no hard custom then temporarily switch to
                    // 1 thread
                    if (scoring_functions[0].any_soft_custom_constraints &&
                        !scoring_functions[0].any_hard_custom_constraints) {
                        pool.setNumThreads(0);
                    }
                    // start timing the smc split
                    auto smc_weight_start_time = std::chrono::steady_clock::now();
                    if (wgt_type == "optimal") {
                        // TODO make more princicpal in the future
                        // for now its just if not district only and not final round
                        if (verbosity >= 3)
                            Rprintf("Computing Optimal Weights:\n");
                        compute_all_plans_log_optimal_incremental_weights(
                            pool, map_params, *splitting_schedule_ptr, sampling_space,
                            scoring_functions, rho, entire_map_compactness,
                            plan_ensemble_ptr->plan_ptr_vec, tree_splitter_ptrs_vec,
                            compute_log_splitting_prob, is_final_splitting_step,
                            smc_diagnostics.log_incremental_weights_mat.col(smc_step_num),
                            *cache_ensemble_ptr, 
                            smc_diagnostics, smc_step_num, step_num,
                            verbosity);
                    } else if (wgt_type == "simple") {
                        if (verbosity >= 3)
                            Rprintf("Computing Simple Backwards Kernel Weights:\n");
                        compute_all_plans_log_simple_incremental_weights(
                            pool, map_params, *splitting_schedule_ptr, sampling_space,
                            scoring_functions, rho, plan_ensemble_ptr->plan_ptr_vec,
                            tree_splitter_ptrs_vec, compute_log_splitting_prob,
                            is_final_splitting_step,
                            smc_diagnostics.log_incremental_weights_mat.col(smc_step_num),
                            verbosity);
                    } else {
                        throw Rcpp::exception("invalid weight type!");
                    }
                    // end timing 
                    auto smc_weight_end_time = std::chrono::steady_clock::now();
                    // add the time 
                    std::chrono::duration<double, std::ratio<1>> smc_weight_diff = smc_weight_end_time - smc_weight_start_time;
                    smc_diagnostics.smc_weight_times[smc_step_num] = smc_weight_diff.count();

                    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                        Rprintf("Done computing weights!\n");

                    // now swap back
                    if (scoring_functions[0].any_soft_custom_constraints &&
                        !scoring_functions[0].any_hard_custom_constraints) {
                        pool.setNumThreads(num_threads);
                    }

                    // do seq_alpha if needed
                    if (apply_weights_alpha) {
                        // reindex first by making the log weight at index i the index of the
                        // parent we use unnormalized_sampling_weights to hold new values then
                        // swap later
                        for (size_t i = 0; i < nsims; i++) {
                            unnormalized_sampling_weights[i] =
                                log_weights[smc_diagnostics.parent_index_mat(i, smc_step_num)];
                        }
                        std::swap(unnormalized_sampling_weights, log_weights);
                        // add incremental weights to the current log weights
                        log_weights =
                            log_weights +
                            smc_diagnostics.log_incremental_weights_mat.col(smc_step_num);

                        // if using seq_alpha then our sampling weights for next round are
                        // proportional to exp(alpha* (prev_log_weights + incremental_weights))
                        unnormalized_sampling_weights = arma::exp(weights_alpha * log_weights);
                        if (!is_final_splitting_step) {
                            // if not the end then multiply by 1-alpha
                            log_weights = (1 - weights_alpha) * log_weights;
                        }
                    } else {
                        // if no seq alpha then log weights are just the incremental weights
                        // and sampling weights are just exp of exponential weights
                        log_weights =
                            smc_diagnostics.log_incremental_weights_mat.col(smc_step_num);
                        unnormalized_sampling_weights = arma::exp(log_weights);
                    }
                    normalized_cumulative_weights = arma::cumsum(unnormalized_sampling_weights);

                    // compute log weight sd
                    smc_diagnostics.log_wgt_stddevs.at(smc_step_num) =
                        arma::stddev(log_weights);
                    // compute effective sample size
                    smc_diagnostics.n_eff.at(smc_step_num) =
                        normalized_cumulative_weights[nsims - 1] *
                        normalized_cumulative_weights[nsims - 1] /
                        arma::sum(arma::square(unnormalized_sampling_weights));
                    // Now normalize the weights
                    normalized_cumulative_weights = normalized_cumulative_weights /
                                                    normalized_cumulative_weights[nsims - 1];

                    if (verbosity >= 3) {
                        Rcpp::Rcout << "  " << std::setprecision(2)
                              << 100 * smc_diagnostics.n_eff.at(smc_step_num) / nsims
                              << "% efficiency." << std::setprecision(4)
                              << " Log Weight Standard Deviation: "
                              << smc_diagnostics.log_wgt_stddevs.at(smc_step_num) << std::endl;
                        Rcpp::Rcout << "  Calculating log weights took " 
                              << smc_diagnostics.smc_weight_times[smc_step_num]
                              << " seconds."
                              << std::endl;
                    }

                    // only increase if we have smc steps left else it will cause index issues
                    // with merge split
                    if (smc_step_num < total_smc_steps - 1) {
                        smc_step_num++;
                    }
                } else if (merge_split_step_vec[step_num]) { // check if its a merge split step
                    // run merge split
                    // Set the number of steps to run at 1 over previous stage acceptance rate
                    // if not 0
                    int prev_acceptance_index =
                        merge_split_step_num == 0 ? step_num - 1 : step_num - 2;
                    double prev_acceptance_rate =
                        smc_diagnostics.acceptance_rates.at(prev_acceptance_index);
                    // if the acceptance is zero just default to 5
                    prev_acceptance_rate = prev_acceptance_rate > 0 ? prev_acceptance_rate : .1;

                    int nsteps_to_run =
                        std::ceil(ms_steps_multiplier * std::ceil((1 / prev_acceptance_rate)));
                    smc_diagnostics.num_merge_split_attempts_vec.at(merge_split_step_num) =
                        nsteps_to_run;

                    if (verbosity >= 3) {
                        Rprintf("  Running %d Merge Split Steps per plan, %d in total!\n",
                                nsteps_to_run, nsteps_to_run * nsims);
                    }

                    splitting_schedule_ptr->update_cut_sizes_for_mergesplit_step(
                        smc_step_num, plan_ensemble_ptr->plan_ptr_vec[0]->num_regions);

                    // If only soft custom constraints then need to temporarily switch to single
                    // threading
                    if (scoring_functions[0].any_soft_custom_constraints &&
                        !scoring_functions[0].any_hard_custom_constraints) {
                        pool.setNumThreads(0);
                    }
                    
                    // optional integrity check 
                    if constexpr (perf_config::object_integrity_checking){
                        std::ostringstream oss;
                        oss << "Plan integrity check failed BEFORE calling `run_merge_split_step_on_all_plans` ";
                        oss << "on ms_step =" << merge_split_step_num;
                        oss << " total step " << step_num << "\n";
                        plan_ensemble_ptr->check_all_plans_valid(
                            map_params,
                            "oss.str()"
                        );                        
                    }

                    // start timing
                    auto ms_round_start_time = std::chrono::steady_clock::now();
                    run_merge_split_step_on_all_plans(
                        pool, map_params, *splitting_schedule_ptr, scoring_functions,
                        rng_states, sampling_space, plan_ensemble_ptr->plan_ptr_vec,
                        dummy_plan_ensemble_ptr->plan_ptr_vec, tree_splitter_ptrs_vec,
                        merge_prob_type, rho, is_final_splitting_step, nsteps_to_run,
                        merge_split_step_num, step_num, smc_diagnostics, *cache_ensemble_ptr,
                        verbosity);

                    // end timing 
                    auto ms_round_end_time = std::chrono::steady_clock::now();
                    // add the time 
                    std::chrono::duration<double, std::ratio<1>> ms_round_diff = ms_round_end_time - ms_round_start_time;
                    smc_diagnostics.ms_step_times[merge_split_step_num] = ms_round_diff.count();

                    // optional integrity check 
                    if constexpr (perf_config::object_integrity_checking){
                        std::ostringstream oss;
                        oss << "Plan integrity check failed AFTER calling `run_merge_split_step_on_all_plans` ";
                        oss << "on ms_step =" << merge_split_step_num;
                        oss << " total step " << step_num << "\n";
                        plan_ensemble_ptr->check_all_plans_valid(
                            map_params,
                            "oss.str()"
                        );                        
                    }


                    // now switch back
                    if (scoring_functions[0].any_soft_custom_constraints &&
                        !scoring_functions[0].any_hard_custom_constraints) {
                        pool.setNumThreads(num_threads);
                    }

                    if (use_naive_k_splitter) {
                        smc_diagnostics.cut_k_values.at(step_num) =
                            tree_splitter_ptrs_vec[0]->get_single_int_param();
                    }


                    // set the acceptance rate
                    int total_ms_successes = Rcpp::sum(
                        smc_diagnostics.merge_split_successes_mat.column(merge_split_step_num));
                    int total_ms_attempts = nsims * nsteps_to_run;

                    smc_diagnostics.acceptance_rates.at(step_num) =
                        total_ms_successes / static_cast<double>(total_ms_attempts);

                    // add number of unique plans
                    smc_diagnostics.nunique_plans[step_num] =
                        plan_ensemble_ptr->count_unique_plans(pool);

                    if (verbosity >= 3) {
                        Rcpp::Rcout << "  " << std::setprecision(2)
                              << 100.0 * smc_diagnostics.acceptance_rates.at(step_num)
                              << "% acceptance rate. " << "There are now "
                              << smc_diagnostics.nunique_plans[step_num] << " unique plans."
                              << std::endl;
                        Rcpp::Rcout << "  Performing MCMC round took " 
                              << smc_diagnostics.ms_step_times[merge_split_step_num] 
                              << " seconds."
                              << std::endl;
                    }

                    // Access the column
                    Rcpp::IntegerMatrix::Column col = smc_diagnostics.draw_tries_mat(Rcpp::_, step_num);
                    // Set all elements in the column to the value nsteps_to_run
                    std::fill(col.begin(), col.end(), nsteps_to_run);

                    merge_split_step_num++;
                }

                // Add details diagnostics if needed
                // Now update the diagnostic info if needed, region labels, dval column of the
                // matrix
                if (diagnostic_mode) {
                    if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
                        REprintf("Adding full diagnostics!\n");
                    }
                    smc_diagnostics.add_full_step_diagnostics(
                        total_steps, splitting_all_the_way, step_num, merge_split_step_num,
                        smc_step_num, !merge_split_step_vec[step_num], sampling_space, pool,
                        *plan_ensemble_ptr, *dummy_plan_ensemble_ptr, *splitting_schedule_ptr);
                }

                if (verbosity == 1 && CLI_SHOULD_TICK) {
                    cli_progress_set(bar, step_num);
                }
                Rcpp::checkUserInterrupt();
            }
        }  catch (const RcppThread::UserInterruptException &e) {
                cli_progress_done(bar);
                Rcpp::stop("redist_smc was interrupted by the user.");
            } catch (const Rcpp::internal::InterruptedException &e) {
                cli_progress_done(bar);
                Rcpp::stop("redist_smc was interrupted by the user.");
            } catch (const std::exception &e) {
                cli_progress_done(bar);
                Rcpp::stop("redist_smc failed with error message:\n%s", e.what());
            } catch (...) {
                cli_progress_done(bar);
                Rcpp::stop("redist_smc failed with an unknown C++ exception.");
            }
        cli_progress_done(bar);
        if constexpr (DEBUG_GSMC_PLANS_VERBOSE) {
            REprintf("Done with main for loop!\n");
        }

        // if not all diagnostic mode and split district only reorder the plans
        if (!diagnostic_mode &&
            splitting_size_regime != SplittingSizeScheduleType::DistrictOnlySMD) {
            reorder_all_plans(pool, plan_ensemble_ptr->plan_ptr_vec,
                              dummy_plan_ensemble_ptr->plan_ptr_vec);
        }
        // end of scope
    }

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Exiting main loop and going to do diagnostics!\n");

    bool plan_sizes_saved;

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Plans saved!\n");

    // if only sampling partial plans or MMD plans then return the size matrix
    if (!splitting_all_the_way ||
        splitting_schedule_ptr->schedule_type == SplittingSizeScheduleType::AnyValidSizeMMD ||
        splitting_schedule_ptr->schedule_type == SplittingSizeScheduleType::DistrictOnlyMMD) {
        if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
            Rprintf("Getting ready to save region sizes!\n");
        plan_sizes_saved = true;
    } else {
        plan_sizes_saved = false;
    }

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Plan matrix (and sizes potentially) saved!\n");

    // Return results
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["plans_mat"] =
            plan_ensemble_ptr->get_R_plans_matrix(), // integer matrix to store final plans
        Rcpp::_["seats"] = plan_sizes_saved
                         ? plan_ensemble_ptr->get_R_sizes_matrix(pool)
                         : Rcpp::IntegerMatrix(1, 1), // saves sizes matrix if needed
        Rcpp::_["region_pops"] = plan_ensemble_ptr->get_region_pops_matrix(pool),
        Rcpp::_["plan_seats_saved"] = plan_sizes_saved, Rcpp::_["log_weights"] = log_weights,
        Rcpp::_["ancestors"] = ancestors, Rcpp::_["step_types"] = step_types,
        Rcpp::_["merge_split_steps"] = merge_split_step_vec,
        Rcpp::_["log_blank_map_target_density"] = log_blank_map_target_density,
        Rcpp::_["multidistrict_selection_alpha"] = SELECTION_ALPHA 
    );

    // to try to save memory kill the plan vector
    plan_ensemble_ptr->flattened_all_plans.clear();
    plan_ensemble_ptr->flattened_all_plans.shrink_to_fit();
    // add all the diagnostics
    smc_diagnostics.add_diagnostics_to_out_list(out);

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Returning to R!\n");
    return out;
}
