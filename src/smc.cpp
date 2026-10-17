/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2024/10
 * Purpose: Run SMC with MCMC merge-split steps mixed in
 ********************************************************/

constexpr bool DEBUG_GSMC_PLANS_VERBOSE = false; // Compile-time constant


#include <exception>
#include <iomanip>
#include <locale>
#include <thread>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

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
 *  Fill `wgts` with `exp(scale * log_wgts)`, shifted so the largest entry is
 *  exactly 1.
 *
 *  Everything these weights feed is a ratio in which a shared constant factor
 *  cancels (the normalized cumulative weights, and n_eff = (sum w)^2 / sum w^2),
 *  so the shift is mathematically a no-op. It matters numerically: without it
 *  `exp` overflows once a log weight exceeds ~709, and the sum of squares
 *  overflows at only ~355, which would silently turn n_eff into 0 or NaN. After
 *  shifting, the largest weight is 1 so the sum lies in [1, n] and cannot
 *  overflow or be zero.
 *
 *  `log_wgts` is carried across SMC steps and returned to R, so it is left
 *  untouched; only the temporary sampling weights are shifted.
 */
static void fill_shifted_exp_weights(std::vector<double> const &log_wgts,
                                     double const scale,
                                     std::vector<double> &wgts) {
    std::size_t const n = log_wgts.size();

    // computed on the exponent itself so a negative `scale` is handled too
    double max_exponent = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < n; i++) {
        double const exponent = scale * log_wgts[i];
        if (exponent > max_exponent) max_exponent = exponent;
    }

    for (std::size_t i = 0; i < n; i++) {
        wgts[i] = std::exp(scale * log_wgts[i] - max_exponent);
    }
}

/*
 *  Sample standard deviation, with denominator n - 1.
 *
 *  This is the corrected two-pass formula `arma::stddev` used,
 *      var = (sum (m - x)^2 - (sum (m - x))^2 / n) / (n - 1),
 *  where the second term cancels the rounding error left in the mean.
 */
static double compute_stddev(std::vector<double> const &x) {
    std::size_t const n = x.size();
    if (n < 2) return 0.0;

    double sum = 0.0;
    for (double const v : x) sum += v;
    double const mean = sum / static_cast<double>(n);

    double sum_sq_diff = 0.0;
    double sum_diff = 0.0;
    for (double const v : x) {
        double const diff = mean - v;
        sum_sq_diff += diff * diff;
        sum_diff += diff;
    }

    return std::sqrt((sum_sq_diff - sum_diff * sum_diff / static_cast<double>(n)) /
                     static_cast<double>(n - 1));
}

// Wrapper object for all non-essential SMC diagnostics
class SMCDiagnostics {

  public:
    SMCDiagnostics(SamplingSpace const sampling_space,
                   SplittingMethodType const splitting_method_type,
                   SplittingSizeScheduleType const splitting_schedule_type,
                   std::vector<bool> const &merge_split_step_vec, int const V, int const nsims,
                   int const ndists, int const total_seats, int const initial_num_regions,
                   int const total_smc_steps, int const total_ms_steps,
                   bool const estimated_unbiased_normalizing_constant,
                   int const diagnostic_level, bool const splitting_all_the_way,
                   bool const split_district_only, bool const final_resample);

    int const diagnostic_level;
    int const total_steps;
    // Whether a final resampling step will be run. When it is, the per-step
    // diagnostics below carry one extra slot for it.
    bool const final_resample;
    // Level 0
    // Essential Diagnostics that are always created
    std::vector<double> log_wgt_stddevs;  // log weight std devs
    std::vector<double> acceptance_rates; // Tracks the acceptance rate - total number of tries
                                          // over nsims - for each round
    std::vector<int> nunique_parents;     // number of unique parents
    std::vector<int> nunique_plans;       // number of unique plans after each step
    std::vector<double> n_eff; // Tracks the effective sample size for the weights of each round
    // For each merge split step this counts the number of attempts that were made
    std::vector<int> num_merge_split_attempts_vec;
    // Only required for graph sampling
    std::vector<int> cut_k_values; // k value used at each step

    // Level 1
    // These are all nsims by number of smc steps
    Rcpp::NumericMatrix log_incremental_weights_mat; // entry [i][s] is the log unnormalized
                                                     // weight of particle i AFTER split s
    Rcpp::IntegerMatrix draw_tries_mat; // Entry [i][s] is the number of tries it took to form
                                        // particle i on split s
    Rcpp::IntegerMatrix
        parent_index_mat; // Entry [i][s] is the index of the parent of particle i at split s
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of
    // successful merge splits performed for plan i on merge split round s
    Rcpp::IntegerMatrix merge_split_successes_mat;
    // counts the size of the trees
    Rcpp::IntegerMatrix tree_sizes_mat;            // ndists by total_steps matrix
    Rcpp::IntegerMatrix successful_tree_sizes_mat; // ndists by total_steps matrix
    std::vector<int> tries_before_extra_particle; // length number of smc stesp 
    // These store time info
    std::vector<double> smc_step_parameter_estimation_times; // length number of smc steps
    std::vector<double> smc_split_times; // length number of smc steps
    std::vector<double> smc_weight_times; // length number of smc steps
    std::vector<double> ms_step_parameter_estimation_times; // length number of ms rounds
    std::vector<double> ms_step_times; // length number of ms rounds

    // Level 2
    Rcpp::IntegerMatrix parent_unsuccessful_tries_mat;

    // level 3
    std::vector<Rcpp::IntegerMatrix> all_steps_plan_region_ids_list;
    // num_edges by nsims logical matrices, see PlanEnsemble::get_R_forest_matrix
    std::vector<Rcpp::LogicalMatrix> all_steps_forests_list;
    // long format data frames, see PlanEnsemble::get_R_linking_edges
    std::vector<Rcpp::DataFrame> all_steps_linking_edge_list;
    std::vector<Rcpp::IntegerMatrix> region_sizes_mat_list;

    // These are granular time stuff that is only tracked when 
    // TRACK_GRANULAR_PERFORMANCE_TIMES (in redist_constants.h) is set to true 
    std::vector<double> wilson_call_times; // time spent drawing spanning trees with wilson
    std::vector<double> wilson_backfill_call_times; // time spent drawing spanning trees with wilson to replace skipped deterministic ones
    std::vector<double> md_selection_times; // time spent picking a multidistrict to split 
    std::vector<double> plan_updating_times; // Times spent updating a plan object after a split 
    std::vector<double> hard_constraint_split_times; // Time spent checking hard constraints in splitting 
    Rcpp::NumericMatrix total_plan_smc_split_times; // Time spent in smc splitting step per plan
    std::vector<double> get_valid_smc_pairs_times; // Time spent finding adjacent pairs and effective boundary lengths for smc. Just time getting pairs for mcmc
    std::vector<double> get_valid_mergepairs_times; // Time spent constructing the list of valid adjacent pairs to merge for mcmc
    std::vector<double> plan_scores_times; // Time spent computing plan based scores in smc weights and mergesplit
    std::vector<double> region_scores_times; // Time spent computing region based scores in smc weights and mergesplit
    std::vector<double> log_tau_times; // Time spent computing log spanning trees in smc weights and mergesplit 
    std::vector<double> retro_splitting_prob_times; // Time spent computing retroactive splitting prob in smc weights
    Rcpp::NumericMatrix total_plan_smc_weight_times; // total time spent on smc weights 
    std::vector<double> selecting_merge_pair_times; // Time spent related to picking pair to merge 
    std::vector<double> eff_boundary_times; // time spent computing effective boundary lengths 
    Rcpp::NumericMatrix total_plan_mcmc_times; // total time spent running mergesplit on a plan


    // Records the final resampling step in the extra diagnostic slots
    void record_final_resample(std::vector<int> const &resample_index,
                               int const nunique_resampled_plans);

    // adds full diagnostics (takes a lot of memory)
    void add_full_step_diagnostics(int const total_steps, bool const splitting_all_the_way,
                                   int const step_num, int const merge_split_step_num,
                                   SamplingSpace const sampling_space, MapParams const &map_params,
                                   RcppThread::ThreadPool &pool, PlanEnsemble &plan_ensemble,
                                   PlanEnsemble &new_plans_ensemble,
                                   SplittingSchedule const &splitting_schedule);

    // adds full diagnostics for the plans after the final resample as an extra
    // entry after the last step's
    void add_final_resample_full_diagnostics(bool const splitting_all_the_way,
                                             SamplingSpace const sampling_space,
                                             MapParams const &map_params,
                                             RcppThread::ThreadPool &pool,
                                             PlanEnsemble &plan_ensemble,
                                             SplittingSchedule const &splitting_schedule);

    // Updates the out list with all the diagnostics
    void add_diagnostics_to_out_list(Rcpp::List &out);
};



SMCDiagnostics::SMCDiagnostics(SamplingSpace const sampling_space,
                               SplittingMethodType const splitting_method_type,
                               SplittingSizeScheduleType const splitting_schedule_type,
                               std::vector<bool> const &merge_split_step_vec, int const V,
                               int const nsims, int const ndists, int const total_seats,
                               int const initial_num_regions, int const total_smc_steps,
                               int const total_ms_steps,
                               bool const estimated_unbiased_normalizing_constant,
                               int const diagnostic_level,
                               bool const splitting_all_the_way, bool const split_district_only,
                               bool const final_resample)
    : diagnostic_level(diagnostic_level), total_steps(total_smc_steps + total_ms_steps),
      final_resample(final_resample),
      log_wgt_stddevs(total_smc_steps), acceptance_rates(total_steps),
      // the final resampling step gets a slot of its own in anything that
      // tracks how much of the ensemble survived each step
      nunique_parents(total_smc_steps + (final_resample ? 1 : 0)),
      nunique_plans(total_steps + (final_resample ? 1 : 0)), n_eff(total_smc_steps),
      num_merge_split_attempts_vec(total_ms_steps),
      cut_k_values(sampling_space == SamplingSpace::GraphSpace ? total_steps : 0),
      tries_before_extra_particle(estimated_unbiased_normalizing_constant ? total_smc_steps : 0),
      smc_step_parameter_estimation_times(total_smc_steps),
      smc_split_times(total_smc_steps),
      smc_weight_times(total_smc_steps),
      ms_step_parameter_estimation_times(total_ms_steps),
      ms_step_times(total_ms_steps),
      wilson_call_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      wilson_backfill_call_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      md_selection_times(perf_config::track_granular_times ? total_smc_steps : 0),
      plan_updating_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      hard_constraint_split_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      total_plan_smc_split_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_smc_steps : 0
      ),
      get_valid_smc_pairs_times(perf_config::track_granular_times ? total_smc_steps : 0),
      get_valid_mergepairs_times(perf_config::track_granular_times ? total_ms_steps : 0),
      plan_scores_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      region_scores_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      log_tau_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      retro_splitting_prob_times(perf_config::track_granular_times ? total_smc_steps : 0),
      total_plan_smc_weight_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_smc_steps : 0
      ),
      selecting_merge_pair_times(perf_config::track_granular_times ? total_ms_steps : 0),
      eff_boundary_times(perf_config::track_granular_times ? total_ms_steps : 0),
      total_plan_mcmc_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_ms_steps : 0
      )
       {
    // Level 1 Diagnostics. Not too big relative to plan size
    log_incremental_weights_mat = Rcpp::NumericMatrix(
        nsims, total_smc_steps); // entry [i][s] is the log unnormalized
                                 // weight of particle i AFTER split s
    draw_tries_mat =
        Rcpp::IntegerMatrix(nsims, total_steps); // Entry [i][s] is the number of tries it took
                                                 // to form particle i on split s
    parent_index_mat = Rcpp::IntegerMatrix(
        nsims,
        // Entry [i][s] is the index of the parent of particle i at split s. When
        // the run ends with a resampling step it gets one extra column holding
        // the resampling index, which is the parent of particle i after it.
        total_smc_steps + (final_resample ? 1 : 0));
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of
    // successful merge splits performed for plan i on merge split round s
    merge_split_successes_mat =
        Rcpp::IntegerMatrix(total_ms_steps > 0 ? nsims : 1, total_ms_steps);
    // counts the size of the trees
    tree_sizes_mat = Rcpp::IntegerMatrix(total_seats, total_steps);
    successful_tree_sizes_mat = Rcpp::IntegerMatrix(total_seats, total_steps);

    // Level 2
    parent_unsuccessful_tries_mat = Rcpp::IntegerMatrix(nsims, total_smc_steps);

    bool diagnostic_mode = diagnostic_level == 1;
    // level 3
    all_steps_plan_region_ids_list.reserve(diagnostic_mode ? total_steps : 0);
    all_steps_forests_list.resize(
        (diagnostic_mode && sampling_space != SamplingSpace::GraphSpace) ? total_steps : 0);
    all_steps_linking_edge_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::LinkingEdgeSpace) ? total_steps
                                                                               : 0);

    // Store size at every step but last one if needed
    int plan_dval_list_size = (diagnostic_mode & !split_district_only) ? total_steps - 1 : 0;
    if (!splitting_all_the_way)
        plan_dval_list_size++;

    region_sizes_mat_list.reserve(plan_dval_list_size);

    // If diagnostic mode track vertex region ids from every round
    if (diagnostic_mode) {
        // The number of regions starts at 1
        int curr_num_regions = initial_num_regions;
        for (size_t i = 0; i < total_steps; i++) {
            all_steps_plan_region_ids_list.emplace_back(V, nsims);
            // Create V by nsims matrix for the plan
            // This is a vector where every entry is a V by nsims Rcpp::IntegerMatrix

            // increase number of regions by 1 if that step is an smc one
            if (!merge_split_step_vec.at(i))
                curr_num_regions++;

            // If not doing district only splits, and its not the final one or
            // we're only doing partial plans then make size matrix
            if (!split_district_only && (i < total_steps - 1 || !splitting_all_the_way)) {
                // This is number of regions by nsims
                region_sizes_mat_list.emplace_back(curr_num_regions, nsims);
            }
        }
    }
}

void SMCDiagnostics::record_final_resample(std::vector<int> const &resample_index,
                                           int const nunique_resampled_plans) {
    if (!final_resample) {
        throw Rcpp::exception(
            "SMCDiagnostics was not sized for a final resampling step!\n");
    }
    // `resample_index[i]` is the plan that became plan i, which is exactly the
    // parent relationship the rest of this matrix records, so it goes in the
    // extra final column. It is stored 0-indexed and shifted along with every
    // other entry in `add_diagnostics_to_out_list`.
    std::copy(resample_index.begin(), resample_index.end(),
              parent_index_mat.column(parent_index_mat.ncol() - 1).begin());

    std::set<int> const surviving_parents(resample_index.begin(), resample_index.end());
    nunique_parents.back() = static_cast<int>(surviving_parents.size());
    nunique_plans.back() = nunique_resampled_plans;
}

void SMCDiagnostics::add_full_step_diagnostics(
    int const total_steps, bool const splitting_all_the_way, int const step_num,
    int const merge_split_step_num,
    SamplingSpace const sampling_space, MapParams const &map_params,
    RcppThread::ThreadPool &pool, PlanEnsemble &plan_ensemble, PlanEnsemble &new_plans_ensemble,
    SplittingSchedule const &splitting_schedule) {
    // if(diagnostic_mode){ // record if in diagnostic mode and generalized splits
    //  reorder the plans by oldest split if either we'vxe done any merge split or
    //  its generalized region splits

    bool const split_district_only =
        splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;
    if (merge_split_step_num > 0 || !split_district_only) {
        reorder_all_plans(pool, plan_ensemble.plan_ptr_vec, new_plans_ensemble.plan_ptr_vec);
    }

    // Copy the vertex plan matrix
    all_steps_plan_region_ids_list.at(step_num) = plan_ensemble.get_R_plans_matrix();

    // store the spanning forests and linking edges if the sampling space has them
    if (sampling_space != SamplingSpace::GraphSpace) {
        all_steps_forests_list.at(step_num) =
            plan_ensemble.get_R_forest_matrix(map_params.num_edges, pool);
    }
    if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
        all_steps_linking_edge_list.at(step_num) =
            plan_ensemble.get_R_linking_edges(map_params.graph_edge_index);
    }

    // Copy the sizes if neccesary
    if (!split_district_only && (step_num < total_steps - 1 || !splitting_all_the_way)) {
        region_sizes_mat_list.at(step_num) = plan_ensemble.get_R_sizes_matrix(pool);
    }

    return;
}

void SMCDiagnostics::add_final_resample_full_diagnostics(
    bool const splitting_all_the_way, SamplingSpace const sampling_space,
    MapParams const &map_params, RcppThread::ThreadPool &pool, PlanEnsemble &plan_ensemble,
    SplittingSchedule const &splitting_schedule) {
    if (!final_resample) {
        throw Rcpp::exception(
            "SMCDiagnostics was not sized for a final resampling step!\n");
    }
    bool const split_district_only =
        splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;

    all_steps_plan_region_ids_list.push_back(plan_ensemble.get_R_plans_matrix());
    if (sampling_space != SamplingSpace::GraphSpace) {
        all_steps_forests_list.push_back(
            plan_ensemble.get_R_forest_matrix(map_params.num_edges, pool));
    }
    if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
        all_steps_linking_edge_list.push_back(
            plan_ensemble.get_R_linking_edges(map_params.graph_edge_index));
    }
    // same rule as the last step, sizes are only needed for partial plans
    if (!split_district_only && !splitting_all_the_way) {
        region_sizes_mat_list.push_back(plan_ensemble.get_R_sizes_matrix(pool));
    }
}

void SMCDiagnostics::add_diagnostics_to_out_list(Rcpp::List &out) {
    // make parent index 1 indexed in place
    std::transform(parent_index_mat.begin(), parent_index_mat.end(), parent_index_mat.begin(),
                   [](int x) { return x + 1; });

    // add granular time info 
    Rcpp::List granular_timing;
    if constexpr (perf_config::track_granular_times){
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = true,
            Rcpp::_["wilson_call_times"] = wilson_call_times,
            Rcpp::_["wilson_backfill_call_times"] = wilson_backfill_call_times,
            Rcpp::_["multidistrict_selection_times"] = md_selection_times,
            Rcpp::_["plan_updating_times"] = plan_updating_times,
            Rcpp::_["hard_constraint_split_times"] = hard_constraint_split_times,
            Rcpp::_["total_plan_smc_split_times"] = total_plan_smc_split_times,
            Rcpp::_["get_valid_smc_pairs_times"] = get_valid_smc_pairs_times,
            Rcpp::_["getting_valid_mergepairs_times"] = get_valid_mergepairs_times,
            Rcpp::_["computing_plan_scores_times"] = plan_scores_times,
            Rcpp::_["computing_region_scores_times"] = region_scores_times,
            Rcpp::_["computing_spanning_tree_count_times"] = log_tau_times,   
            Rcpp::_["computing_retro_splitting_prob_times"] = retro_splitting_prob_times,     
            Rcpp::_["total_plan_smc_weight_times"] = total_plan_smc_weight_times,
            Rcpp::_["selecting_merge_pair"] = selecting_merge_pair_times,
            Rcpp::_["eff_boundary_times"] = eff_boundary_times,
            Rcpp::_["total_plan_mcmc_times"] = total_plan_mcmc_times
        );

    }else{
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = false
        );
    }

    // optional add special timing 
    if constexpr (perf_config::special_timing){
        Rcpp::List special_timing_list = Rcpp::List::create(
        );

        granular_timing["special_timing_list"] = special_timing_list;
    }
    
    out["acceptance_rates"] = acceptance_rates;
    out["draw_tries_mat"] = draw_tries_mat;
    out["parent_index"] = parent_index_mat;
    out["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat;
    out["step_n_eff"] = n_eff;
    out["nunique_parent_indices"] = nunique_parents;
    out["nunique_plans"] = nunique_plans;
    out["tree_sizes"] = tree_sizes_mat;
    out["successful_tree_sizes"] = successful_tree_sizes_mat;
    out["log_weight_stddev"] = log_wgt_stddevs;
    out["cut_k_vals"] = cut_k_values;
    out["log_incremental_weights_mat"] = log_incremental_weights_mat;
    out["ms_step_counts"] = num_merge_split_attempts_vec;
    out["merge_split_success_mat"] = merge_split_successes_mat;
    out["tries_before_extra_particle"] = tries_before_extra_particle;
    out["smc_step_parameter_estimation_times"] = smc_step_parameter_estimation_times;
    out["smc_split_times"] = smc_split_times;
    out["smc_weight_times"] = smc_weight_times;
    out["ms_step_parameter_estimation_times"] = ms_step_parameter_estimation_times;
    out["ms_step_times"] = ms_step_times;
    out["granular_times"] = granular_timing;
    out["region_ids_mat_list"] = all_steps_plan_region_ids_list;
    out["region_seats_mat_list"] = region_sizes_mat_list;
    out["forest_adjs_list"] = all_steps_forests_list;
    out["linking_edges_list"] = all_steps_linking_edge_list;


    return;
}



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
 *
 */
void run_smc_step(const MapParams &map_params, SplittingSchedule const &splitting_schedule,
                  std::vector<ScoringFunction> const &scoring_functions,
                  std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
                  std::unique_ptr<PlanEnsemble> &old_plan_ensemble,
                  std::unique_ptr<PlanEnsemble> &new_plan_ensemble,
                  std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters,
                  const std::vector<double> &normalized_cumulative_weights,
                  SMCDiagnostics &smc_diagnostics, int const smc_step_num, int const step_num,
                  bool const is_final_split,
                  bool const estimated_unbiased_normalizing_constant,
                  RcppThread::ThreadPool &pool, int verbosity, int diagnostic_level,
                  int const max_split_tries, double const multidistrict_selection_alpha) {
    // important constants
    int const num_threads = get_num_threads(pool);
    const int M = old_plan_ensemble->nsims;
    bool const smd_split_district_only =
        splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;

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
    std::vector<double> backfill_wilson_times (
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
    // one schedule context per thread, rebuilt for each plan a thread tries to split
    std::vector<ScheduleContext> schedule_ctxs(num_threads);


    // sample extra plan if needed
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
                "Failed to split a plan after `max_split_tries` attempts!\n");
        }
        // increase the number of tries for particle i by 1
        Rcpp::checkUserInterrupt();
        extra_particle_reject_ct++;
        // sample previous plan
        int idx = rng_states[0].r_int_wgt(normalized_cumulative_weights);
        Plan const &old_plan = *old_plan_ensemble->plan_ptr_vec[idx];
        splitting_schedule.build_context(old_plan, schedule_ctxs[0]);
        // Get region id the split
        int region_id_to_split;
        if (smd_split_district_only) {
            // if just doing district splits just use remainder region
            // which is always the highest id
            region_id_to_split = old_plan.num_regions - 1;
        } else {
            // if generalized split pick a region to try to split
            region_id_to_split = old_plan.choose_multidistrict_to_split(
                splitting_schedule, schedule_ctxs[0], rng_states[0],
                multidistrict_selection_alpha);
        }

        // Try to split the region
        CutSizeOptions const cut_options =
            splitting_schedule.region_cut_sizes(old_plan, schedule_ctxs[0], region_id_to_split);
        std::pair<bool, EdgeCut> edge_search_result =
            ust_samplers_vec[0].attempt_to_find_valid_tree_split(
                rng_states[0], scoring_functions[0],
                *tree_splitters[0], old_plan,
                region_id_to_split, new_region_id, cut_options, save_edge_selection_prob);



        // if successful update the new plan and check if satisfies any other hard
        // constraints
        if (std::get<0>(edge_search_result)) {
            // check if there are any additional hard constraints
            if (!scoring_functions[0].any_hard_constraints) {
                // if not we can stop trying 
                extra_plan_sampled = true;
            } else {
                // make the new plan a copy of the old one
                // we can just use 0 since it will be overwritten later 
                new_plan_ensemble->plan_ptr_vec[0]->shallow_copy(
                    *old_plan_ensemble->plan_ptr_vec[idx]);
                // // now split that region we found on the old one
                new_plan_ensemble->plan_ptr_vec[0]->update_from_successful_split(
                    *tree_splitters[0], ust_samplers_vec[0],
                    std::get<1>(edge_search_result), region_id_to_split, new_region_id, true);

                extra_plan_sampled = scoring_functions[0].satisfies_hard_constraints(
                    *new_plan_ensemble->plan_ptr_vec[0], region_id_to_split, new_region_id
                );
            }
        }
        
    }
    if constexpr(DEBUG_GSMC_PLANS_VERBOSE){
        REprintf("Done!\n");
    }
    

    // now save the number of failed attempts before sampling the extra plan
    smc_diagnostics.tries_before_extra_particle[smc_step_num] = extra_particle_reject_ct;

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
            Plan const &old_plan = *old_plan_ensemble->plan_ptr_vec[idx];
            ScheduleContext &schedule_ctx = schedule_ctxs[thread_id];
            splitting_schedule.build_context(old_plan, schedule_ctx);

            // Get region id the split
            int region_id_to_split;
            auto md_selection_time = maybe_now();
            if (smd_split_district_only) {
                // if just doing district splits just use remainder region
                // which is always the highest id
                region_id_to_split = old_plan.num_regions - 1;
            } else {
                // if generalized split pick a region to try to split
                region_id_to_split = old_plan.choose_multidistrict_to_split(
                    splitting_schedule, schedule_ctx, rng_states[thread_id],
                    multidistrict_selection_alpha);
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
            CutSizeOptions const cut_options =
                splitting_schedule.region_cut_sizes(old_plan, schedule_ctx, region_id_to_split);
            std::pair<bool, EdgeCut> edge_search_result =
                ust_samplers_vec[thread_id].attempt_to_find_valid_tree_split(
                    rng_states[thread_id], scoring_functions[thread_id],
                    *tree_splitters[thread_id], old_plan,
                    region_id_to_split, new_region_id, cut_options, save_edge_selection_prob);
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
                auto backfill_wilson_time = maybe_now();
                // Now fill in the deterministic subtrees
                new_plan_ensemble->plan_ptr_vec[i]->fill_in_skipped_subtrees(
                    ust_samplers_vec[thread_id], rng_states[thread_id]
                );
                if constexpr (perf_config::track_granular_times){
                    add_elapsed(backfill_wilson_times[thread_id], backfill_wilson_time); // optional timing
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
            smc_diagnostics.wilson_backfill_call_times[step_num] += backfill_wilson_times[thread_id];
            smc_diagnostics.md_selection_times[smc_step_num] += md_selection_times[thread_id];
            smc_diagnostics.hard_constraint_split_times[step_num] += hard_constraint_split_times[thread_id];
            smc_diagnostics.plan_updating_times[step_num] += plan_updating_times[thread_id];
        }
    }

    

    // now compute acceptance rate and unique parents and original ancestors
    // int64 to avoid overflow
    std::int64_t total_split_attempts =
        std::accumulate(
            draw_tries_buffer.begin(),
            draw_tries_buffer.end(),
            std::int64_t{0}
        );

    double accept_rate =
        static_cast<double>(M) /
        static_cast<double>(total_split_attempts);

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
            smc_diagnostics.wilson_backfill_call_times[step_num] += granular_weight_times[thread_id].backfill_wilson_time;
            smc_diagnostics.get_valid_mergepairs_times[merge_split_step_num] += granular_weight_times[thread_id].get_valid_mergepairs;
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




// NEED TO UPDATE THIS IS OLD DOCUMENTATION FOR GRAPH STUFF
// Computes log unnormalized weights for vector of plans
//
// Using the procedure outlined in <PAPER HERE> this function computes the log
// incremental weights and the unnormalized weights for a vector of plans (which
// may or may not be the same depending on the parameters).
//
// @title Compute Log Unnormalized Weights
//
// @param pool A threadpool for multithreading
// @param g A graph (adjacency list) passed by reference
// @param plans_ptr_vec A vector of plans to compute the log unnormalized weights
// of
// @param split_district_only whether or not to compute the weights under
// the district only split scheme or not. If `split_district_only` is true
// then uses optimal weights from one-district split scheme.
// @param log_incremental_weights A vector of the log incremental weights
// computed for the plans. The value of `log_incremental_weights[i]` is
// the log incremental weight for `plans_ptr_vec[i]`
// @param unnormalized_sampling_weights A vector of the unnormalized sampling
// weights to be used with sampling the `plans_ptr_vec` in the next iteration of the
// algorithm. Depending on the other hyperparameters this may or may not be the
// same as `exp(log_incremental_weights)`
// @param target Target population of a single district
// @param pop_temper <DETAILS NEEDED>
//
// @details Modifications
//    - The `log_incremental_weights` is updated to contain the incremental
//    weights of the plans
//    - The `unnormalized_sampling_weights` is updated to contain the unnormalized
//    sampling weights of the plans for the next round
void compute_all_plans_log_optimal_incremental_weights(
    RcppThread::ThreadPool &pool, const MapParams &map_params,
    const SplittingSchedule &splitting_schedule, SamplingSpace const sampling_space,
    std::vector<ScoringFunction> const &scoring_functions, double rho,
    double const whole_map_compactness_term, std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    std::vector<std::unique_ptr<TreeSplitter>> &tree_splitter_ptrs_vec,
    bool compute_log_splitting_prob, 
    double const multidistrict_selection_alpha,
    bool is_final_plans, 
    std::vector<double> &log_incremental_weights, WeightCacheEnsemble &cache_ensemble,
    SMCDiagnostics &smc_diagnostics, int const smc_step_num, int const step_num,
    int verbosity) {
    const int nsims = static_cast<int>(plans_ptr_vec.size());
    const int check_int = 50; // check for interrupts every _ iterations

    int const num_threads = get_num_threads(pool);
    // thread safe id counter
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};
    std::vector<std::atomic<int>> active_users(
        perf_config::check_threadpool_integrity ? num_threads : 0);

    // now make the vectors of important variables to be used by threads
    std::vector<USTSampler> ust_samplers_vec;
    ust_samplers_vec.reserve(num_threads);
    std::vector<PlanMultigraph> plan_multigraphs_vec;
    plan_multigraphs_vec.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++) {
        ust_samplers_vec.emplace_back(map_params, splitting_schedule);
        plan_multigraphs_vec.emplace_back(map_params,
                                          sampling_space == SamplingSpace::LinkingEdgeSpace);
    }
    std::vector<GranularWeightTimes> granular_weight_times(num_threads);

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

        auto total_weight_time = maybe_now(); // optional timing 


        if (cache_ensemble.using_caching) {
            log_incremental_weights[i] = compute_log_optimal_incremental_weights(
                *plans_ptr_vec[i], plan_multigraphs_vec[thread_id], splitting_schedule,
                ust_samplers_vec[thread_id], *tree_splitter_ptrs_vec[thread_id], sampling_space,
                scoring_functions[thread_id], rho, whole_map_compactness_term,
                compute_log_splitting_prob, multidistrict_selection_alpha,
                is_final_plans, cache_ensemble.using_caching,
                cache_ensemble.weight_cache_ptr_vec[i].get(), 
                granular_weight_times[thread_id]);
        } else {
            log_incremental_weights[i] = compute_log_optimal_incremental_weights(
                *plans_ptr_vec[i], plan_multigraphs_vec[thread_id], splitting_schedule,
                ust_samplers_vec[thread_id], *tree_splitter_ptrs_vec[thread_id], sampling_space,
                scoring_functions[thread_id], rho, whole_map_compactness_term,
                compute_log_splitting_prob, multidistrict_selection_alpha,
                is_final_plans, cache_ensemble.using_caching,
                nullptr, granular_weight_times[thread_id]);
        }

        

        if constexpr (perf_config::track_granular_times){
            add_elapsed(
                smc_diagnostics.total_plan_smc_weight_times(i, smc_step_num), 
                total_weight_time
            ); // optional timing 
        }
        

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);

    });

    // Wait for all the threads to finish
    pool.wait();

    // add granular time if that's being tracked 
    if constexpr (perf_config::track_granular_times){
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
        {
            smc_diagnostics.get_valid_smc_pairs_times[smc_step_num] += granular_weight_times[thread_id].get_valid_pairs;
            smc_diagnostics.retro_splitting_prob_times[smc_step_num] += granular_weight_times[thread_id].splitting_prob;
            smc_diagnostics.plan_scores_times[step_num] += granular_weight_times[thread_id].plan_scores;
            smc_diagnostics.region_scores_times[step_num] += granular_weight_times[thread_id].region_scores;
            smc_diagnostics.log_tau_times[step_num] += granular_weight_times[thread_id].tau_terms;
        }
    }

    return;
}


void compute_all_plans_log_simple_incremental_weights(
    RcppThread::ThreadPool &pool, const MapParams &map_params,
    const SplittingSchedule &splitting_schedule, SamplingSpace const sampling_space,
    std::vector<ScoringFunction> const &scoring_functions, double rho,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec,
    std::vector<std::unique_ptr<TreeSplitter>> &tree_splitter_ptrs_vec,
    bool compute_log_splitting_prob, double const multidistrict_selection_alpha,
    bool is_final_plans,
    std::vector<double> &log_incremental_weights, WeightCacheEnsemble &cache_ensemble,
    SMCDiagnostics &smc_diagnostics, int const smc_step_num, int const step_num,
    int verbosity) {
    int const nsims = (int)plans_ptr_vec.size();
    const int check_int = 50; // check for interrupts every _ iterations

    int const num_threads = get_num_threads(pool);
    // thread safe id counter
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};
    std::vector<std::atomic<int>> active_users(
        perf_config::check_threadpool_integrity ? num_threads : 0);
    // now make the vectors of important variables to be used by threads
    std::vector<USTSampler> ust_samplers_vec;
    ust_samplers_vec.reserve(num_threads);
    std::vector<PlanMultigraph> plan_multigraphs_vec;
    plan_multigraphs_vec.reserve(num_threads);
    for (size_t i = 0; i < num_threads; i++) {
        ust_samplers_vec.emplace_back(map_params, splitting_schedule);
        plan_multigraphs_vec.emplace_back(map_params,
                                          sampling_space == SamplingSpace::LinkingEdgeSpace);
    }
    std::vector<GranularWeightTimes> granular_weight_times(num_threads);

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

        auto total_weight_time = maybe_now(); // optional timing

        log_incremental_weights[i] = compute_simple_log_incremental_weight(
            *plans_ptr_vec[i], plan_multigraphs_vec[thread_id], splitting_schedule,
            ust_samplers_vec[thread_id], *tree_splitter_ptrs_vec[thread_id], sampling_space,
            scoring_functions[thread_id], rho, compute_log_splitting_prob, multidistrict_selection_alpha,
            is_final_plans, cache_ensemble.using_caching,
            cache_ensemble.using_caching ? cache_ensemble.weight_cache_ptr_vec[i].get() : nullptr,
            granular_weight_times[thread_id]);

        if constexpr (perf_config::track_granular_times) {
            add_elapsed(
                smc_diagnostics.total_plan_smc_weight_times(i, smc_step_num),
                total_weight_time
            ); // optional timing
        }

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();

    // add granular time if that's being tracked
    if constexpr (perf_config::track_granular_times){
        for (size_t thread_id = 0; thread_id < num_threads; thread_id++)
        {
            smc_diagnostics.get_valid_smc_pairs_times[smc_step_num] += granular_weight_times[thread_id].get_valid_pairs;
            smc_diagnostics.retro_splitting_prob_times[smc_step_num] += granular_weight_times[thread_id].splitting_prob;
            smc_diagnostics.plan_scores_times[step_num] += granular_weight_times[thread_id].plan_scores;
            smc_diagnostics.region_scores_times[step_num] += granular_weight_times[thread_id].region_scores;
            smc_diagnostics.log_tau_times[step_num] += granular_weight_times[thread_id].tau_terms;
        }
    }

    return;
}

/*
 *  Resample the ensemble with a systematic (low-variance) resampler.
 *
 *  This runs at the end of the sampler rather than in R because everything it
 *  needs is already allocated here. `dummy_plan_ensemble` is the scratch
 *  ensemble the SMC steps have been splitting into, so gathering into it and
 *  swapping costs no extra memory, and the plans are exported to R exactly
 *  once, already in their final order. Resampling in R instead means holding
 *  the exported plan matrix plus a permutation buffer, and forces a `cbind`
 *  that reallocates the whole parent index matrix.
 *
 *  A gather cannot be done in place because `resample_index` repeats entries,
 *  which is why the scratch ensemble is needed at all.
 *
 *  @param log_wgts The unnormalized log weights of the plans.
 *
 *  @details Modifications
 *     - `plan_ensemble` holds the resampled plans on return, and
 *     `dummy_plan_ensemble` holds the pre-resampling plans
 *
 *  @returns The 0-indexed resampling index, where entry `i` is the plan that
 *  became plan `i`.
 */
static std::vector<int> resample_plan_ensemble(
    RcppThread::ThreadPool &pool, RNGState &rng_state, std::vector<double> const &log_wgts,
    std::unique_ptr<PlanEnsemble> &plan_ensemble,
    std::unique_ptr<PlanEnsemble> &dummy_plan_ensemble) {
    int const nsims = plan_ensemble->nsims;

    // Exponentiate shifted so the largest weight is 1, which keeps the sum
    // finite, then turn the weights into a normalized cumulative sum. The shift
    // cancels in the normalization.
    std::vector<double> cumulative_wgts(nsims);
    fill_shifted_exp_weights(log_wgts, 1.0, cumulative_wgts);
    double weight_total = 0.0;
    for (int i = 0; i < nsims; i++) {
        weight_total += cumulative_wgts[i];
        cumulative_wgts[i] = weight_total;
    }
    for (double &w : cumulative_wgts) w /= weight_total;

    // One uniform draw placed into `nsims` evenly spaced strata
    std::vector<int> resample_index(nsims);
    double const r = rng_state.r_unif() / nsims;
    int idx = 0;
    for (int n = 0; n < nsims; n++) {
        double const u = r + n / static_cast<double>(nsims);
        // the bound on `idx` stops rounding in the cumulative sum from walking
        // off the end when `u` lands just past the final weight
        while (u > cumulative_wgts[idx] && idx + 1 < nsims) ++idx;
        resample_index[n] = idx;
    }

    pool.parallelFor(0, nsims, [&](int i) {
        dummy_plan_ensemble->plan_ptr_vec[i]->shallow_copy(
            *plan_ensemble->plan_ptr_vec[resample_index[i]]);
    });
    pool.wait();

    // the gathered plans are in the dummy ensemble, so swap it in
    std::swap(plan_ensemble, dummy_plan_ensemble);

    return resample_index;
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
    Rcpp::IntegerMatrix const &region_sizes_mat) {
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
    // weight type
    std::string wgt_type = Rcpp::as<std::string>(control["weight_type"]);
    // whether or not to cache the weights
    bool const using_caching = Rcpp::as<bool>(control["cache_weights"]);
    // max tries value
    int const max_split_tries = Rcpp::as<int>(control["max_split_tries"]);
    // unbiased normalizing estimate
    bool const estimated_unbiased_normalizing_constant = Rcpp::as<bool>(control["est_norm_unbiased"]);
    // whether to resample the plans once the run is over
    bool const final_resample = Rcpp::as<bool>(control["resample"]);
    // Values carried over from an earlier run that this run continues
    // (`keep_init_plan_info`) so it behaves like a single run through every step.
    // - continue_from_init_plans: the initial weights are what a single run would
    //   carry into this step, so they are split between sampling and carrying
    //   with seq_alpha like any other step
    // - initial_last_cut_k: the earlier run's last estimated cut k, which the
    //   first cut k estimate starts its search from
    // - initial_ms_acceptance_rate: the earlier run's last mergesplit acceptance
    //   rate, which sets the number of steps in this run's first mergesplit round
    bool const continue_from_init_plans =
        control.containsElementNamed("continue_from_init_plans") &&
        Rcpp::as<bool>(control["continue_from_init_plans"]);
    // where the first cut k estimate starts its search
    int const initial_last_cut_k = control.containsElementNamed("initial_last_cut_k")
                                       ? Rcpp::as<int>(control["initial_last_cut_k"])
                                       : std::max(1, V - 5);
    bool const has_initial_ms_acceptance_rate =
        control.containsElementNamed("initial_ms_acceptance_rate");
    double const initial_ms_acceptance_rate =
        has_initial_ms_acceptance_rate
            ? Rcpp::as<double>(control["initial_ms_acceptance_rate"])
            : 0.0;
    // The acceptance rate of the most recent mergesplit round, which sets the
    // number of steps in the next round. Before any mergesplit round has been
    // run (here or in an earlier run this continues) the previous step is used.
    bool have_last_ms_acceptance_rate = has_initial_ms_acceptance_rate;
    double last_ms_acceptance_rate = initial_ms_acceptance_rate;
    // whether to return the final spanning forests and linking edges
    bool const return_augmented_samples =
        control.containsElementNamed("return_augmented_samples") &&
        Rcpp::as<bool>(control["return_augmented_samples"]);
    if (return_augmented_samples && sampling_space == SamplingSpace::GraphSpace) {
        throw Rcpp::exception("Augmented samples are not available in graph space!\n");
    }
    double tmp_multidistrict_selection_alpha;
    // custom multidistrict selection alpha
    if(control.containsElementNamed("md_alpha")){
        tmp_multidistrict_selection_alpha = Rcpp::as<double>(control["md_alpha"]); 
    }else{
        tmp_multidistrict_selection_alpha = SELECTION_ALPHA;
    }

    double const multidistrict_selection_alpha = tmp_multidistrict_selection_alpha;

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
            throw std::runtime_error("Invalid step type \"" +
                                     static_cast<std::string>(step_types.at(i)) +
                                     "\" at step " + std::to_string(i + 1) +
                                     ". Step types must be \"smc\" or \"ms\".\n");
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
        splitting_size_regime);
    // it wants presplit number of regions so make initial regions - 1
    // Needed for initializing linking edge plans. A blank map has no previous
    // split so there is nothing to set (and a presplit count of 0 is invalid)
    if (initial_num_regions > 1) {
        splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
            0, initial_num_regions - 1);
    }

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
        diagnostic_level, splitting_all_the_way, split_district_only, final_resample);

    // Create a threadpool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);

    // If hard custom then switch to no threading
    if (scoring_functions[0].any_hard_custom_constraints) {
        pool.setNumThreads(0);
    }

    // The log weights of the plans carried between steps. Any initial weights
    // are only used to sample the parents of the first step so the particles
    // after that are equally weighted and these start at 0, unless this run
    // continues an earlier one with seq_alpha (see the initial sampling weights).
    std::vector<double> log_weights(nsims, 0.0);

    // Now we add everything here to a scope since it won't be needed for the end
    // create the ensemble, using any saved forests and linking edges passed in
    std::unique_ptr<PlanEnsemble> plan_ensemble_ptr = get_plan_ensemble_ptr(
        map_params, *splitting_schedule_ptr, initial_num_regions, nsims, sampling_space,
        region_id_mat, region_sizes_mat, rng_states, pool, verbosity,
        get_initial_augmented_samples(control, map_params, nsims));

    // make sure any partial plans passed in can actually be completed
    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (initial_num_regions > 1) {
            for (auto const &plan_ptr : plan_ensemble_ptr->plan_ptr_vec) {
                splitting_schedule_ptr->check_initial_plan(*plan_ptr);
            }
        }
    }

    // compute the log of the unnormalized density of the entire map 
    // which is log spanning tree count  - score 
    double log_blank_map_target_density = NA_REAL;

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

        // The first step samples parents with the initial weights if they are
        // passed in through control (assumed to already be normalized) else
        // all plans are weighted equally
        std::vector<double> unnormalized_sampling_weights(nsims, 1.0);
        if (control.containsElementNamed("initial_weights")) {
            unnormalized_sampling_weights =
                Rcpp::as<std::vector<double>>(control["initial_weights"]);
            // When continuing an earlier run the initial weights are the full
            // weights a single run would have after this point, so like any other
            // step the parents are sampled proportional to w^alpha and
            // (1 - alpha) * log(w) is carried forward
            if (continue_from_init_plans && apply_weights_alpha) {
                for (int i = 0; i < nsims; i++) {
                    double const log_initial_weight =
                        std::log(unnormalized_sampling_weights[i]);
                    log_weights[i] = (1 - weights_alpha) * log_initial_weight;
                    unnormalized_sampling_weights[i] =
                        std::exp(weights_alpha * log_initial_weight);
                }
            }
        }
        // Reused each SMC step to collect that step's incremental weights before
        // they are copied into the diagnostics matrix.
        std::vector<double> log_incremental_weights(nsims);
        // now get initial normalized weights
        std::vector<double> normalized_cumulative_weights(nsims);
        double weight_total = 0.0;
        for (size_t i = 0; i < nsims; i++) {
            weight_total += unnormalized_sampling_weights[i];
            normalized_cumulative_weights[i] = weight_total;
        }
        for (double &w : normalized_cumulative_weights) w /= weight_total;

        // Create the weight cache's if needed
        std::unique_ptr<WeightCacheEnsemble> cache_ensemble_ptr =
            std::make_unique<WeightCacheEnsemble>(using_caching, map_params, nsims, rho,
                                                  sampling_space);

        std::unique_ptr<WeightCacheEnsemble> dummy_cache_ensemble_ptr =
            std::make_unique<WeightCacheEnsemble>(using_caching, map_params, nsims, rho,
                                                  sampling_space);

        // If the run ends in mergesplit rounds the final resampling is done
        // after the last split instead of at the very end, so those rounds run
        // on the resampled plans and can rebuild the diversity resampling loses.
        // The first step is always a split so there is always a last split.
        int last_smc_step_num = total_steps - 1;
        while (merge_split_step_vec[last_smc_step_num]) --last_smc_step_num;
        bool const resample_before_final_ms =
            final_resample && last_smc_step_num < total_steps - 1;

        // Resamples the plans with the current weights and records it. The
        // weight caches are reordered the same way so they stay with their plans.
        auto do_final_resample = [&]() {
            std::vector<int> const resample_index = resample_plan_ensemble(
                pool, rng_states[0], log_weights, plan_ensemble_ptr,
                dummy_plan_ensemble_ptr);
            if (using_caching) {
                for (int i = 0; i < nsims; i++) {
                    dummy_cache_ensemble_ptr->weight_cache_ptr_vec[i]->copy_from(
                        *cache_ensemble_ptr->weight_cache_ptr_vec[resample_index[i]]);
                }
                std::swap(cache_ensemble_ptr, dummy_cache_ensemble_ptr);
            }
            smc_diagnostics.record_final_resample(
                resample_index, plan_ensemble_ptr->count_unique_plans(pool));
            if (verbosity >= 3) {
                Rcpp::Rcout << "Resampling: " << std::setprecision(2)
                            << 100.0 * smc_diagnostics.nunique_parents.back() / nsims
                            << "% of plans were preserved in the resample, there are now "
                            << smc_diagnostics.nunique_plans.back() << " unique plans."
                            << std::endl;
            }
            // the per step lists are sized to the steps so this is still last
            if (diagnostic_mode) {
                smc_diagnostics.add_final_resample_full_diagnostics(
                    splitting_all_the_way, sampling_space, map_params, pool,
                    *plan_ensemble_ptr, *splitting_schedule_ptr);
            }
        };

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
                int const num_hard_constraints =
                    scoring_functions[0].num_hard_region_constraints +
                    scoring_functions[0].num_hard_plan_constraints;
                if (num_hard_constraints > 0) {
                    Rcpp::Rcout << "Applying " << num_hard_constraints << " hard constraints.\n";
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
                bool const is_final_splitting_step = smc_step_num + 1 == total_smc_steps;
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
                            int last_k = smc_step_num == 0 ? initial_last_cut_k
                                                           : k_params.at(smc_step_num - 1);
                            if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
                                Rprintf("About to try to estimate cut k!\n");
                            estimate_cut_k(map_params, *splitting_schedule_ptr, rng_state,
                                           est_cut_k, last_k,
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
                                 step_num, is_final_splitting_step,
                                 estimated_unbiased_normalizing_constant, pool,
                                 verbosity, diagnostic_mode ? 3 : 0, max_split_tries,
                                 multidistrict_selection_alpha);
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
                            compute_log_splitting_prob, 
                            multidistrict_selection_alpha,
                            is_final_splitting_step,
                            log_incremental_weights,
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
                            multidistrict_selection_alpha,
                            is_final_splitting_step,
                            log_incremental_weights,
                            *cache_ensemble_ptr,
                            smc_diagnostics, smc_step_num, step_num,
                            verbosity);
                    } else {
                        throw std::runtime_error("Invalid weight type \"" + wgt_type +
                                                 "\". Weight type must be \"optimal\" or "
                                                 "\"simple\".\n");
                    }
                    // save this step's incremental weights into the diagnostics
                    std::copy(
                        log_incremental_weights.begin(), log_incremental_weights.end(),
                        smc_diagnostics.log_incremental_weights_mat.column(smc_step_num).begin());
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
                        for (size_t i = 0; i < nsims; i++) {
                            log_weights[i] += log_incremental_weights[i];
                        }

                        // if using seq_alpha then our sampling weights for next round are
                        // proportional to exp(alpha* (prev_log_weights + incremental_weights))
                        fill_shifted_exp_weights(log_weights, weights_alpha,
                                                 unnormalized_sampling_weights);
                        if (!is_final_splitting_step) {
                            // if not the end then multiply by 1-alpha
                            for (double &w : log_weights) w *= (1 - weights_alpha);
                        }
                    } else {
                        // if no seq alpha then log weights are just the incremental weights
                        // and sampling weights are just exp of exponential weights
                        log_weights = log_incremental_weights;
                        fill_shifted_exp_weights(log_weights, 1.0,
                                                 unnormalized_sampling_weights);
                    }
                    weight_total = 0.0;
                    for (size_t i = 0; i < nsims; i++) {
                        weight_total += unnormalized_sampling_weights[i];
                        normalized_cumulative_weights[i] = weight_total;
                    }

                    // compute log weight sd
                    smc_diagnostics.log_wgt_stddevs.at(smc_step_num) =
                        compute_stddev(log_weights);
                    // compute effective sample size. The sampling weights are
                    // shifted so the largest is 1, which keeps this sum finite.
                    double sum_sq_wgts = 0.0;
                    for (double const w : unnormalized_sampling_weights) {
                        sum_sq_wgts += w * w;
                    }
                    smc_diagnostics.n_eff.at(smc_step_num) =
                        weight_total * weight_total / sum_sq_wgts;
                    // Now normalize the weights
                    for (double &w : normalized_cumulative_weights) w /= weight_total;

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
                    // The first mergesplit round of a run continuing an earlier one
                    // uses the earlier run's last mergesplit acceptance rate if it had
                    // any mergesplit rounds
                    double prev_acceptance_rate =
                        have_last_ms_acceptance_rate
                            ? last_ms_acceptance_rate
                            : smc_diagnostics.acceptance_rates.at(step_num - 1);
                    // if the acceptance is zero just default to 5
                    prev_acceptance_rate = prev_acceptance_rate > 0 ? prev_acceptance_rate : .1;

                    int nsteps_to_run =
                        std::ceil(ms_steps_multiplier * std::ceil((1 / prev_acceptance_rate)));
                    smc_diagnostics.num_merge_split_attempts_vec.at(merge_split_step_num) =
                        nsteps_to_run;

                    if (verbosity >= 3) {
                        std::int64_t total_ms_steps =
                            static_cast<std::int64_t>(nsteps_to_run) * nsims;

                        Rcpp::Rcout
                            << "  Running " << nsteps_to_run
                            << " Merge Split Steps per plan, "
                            << total_ms_steps
                            << " in total!" << std::endl;
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
                    std::int64_t total_ms_successes = 0;

                    for (int i = 0; i < nsims; ++i) {
                        total_ms_successes +=
                            smc_diagnostics.merge_split_successes_mat(
                                i,
                                merge_split_step_num
                            );
                    }
                    // make int64 to avoid overflow 
                    std::int64_t total_ms_attempts =
                        static_cast<std::int64_t>(nsims) * nsteps_to_run;

                    smc_diagnostics.acceptance_rates.at(step_num) =
                        static_cast<double>(total_ms_successes) /
                        static_cast<double>(total_ms_attempts);
                    // the next mergesplit round is sized off this one
                    last_ms_acceptance_rate = smc_diagnostics.acceptance_rates.at(step_num);
                    have_last_ms_acceptance_rate = true;

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
                        sampling_space, map_params, pool,
                        *plan_ensemble_ptr, *dummy_plan_ensemble_ptr, *splitting_schedule_ptr);
                }

                // resample after the last split if mergesplit rounds follow it
                if (resample_before_final_ms && step_num == last_smc_step_num) {
                    do_final_resample();
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

        // Every region of a finished plan must be an allowed district size. The
        // splitting schedule should guarantee this, so this only catches bugs.
        if constexpr (perf_config::redundancy_checks) {
            if (splitting_all_the_way) {
                for (auto const &plan_ptr : plan_ensemble_ptr->plan_ptr_vec) {
                    for (int region_id = 0; region_id < plan_ptr->num_regions; region_id++) {
                        int const region_size = plan_ptr->region_sizes[region_id];
                        if (!map_params.is_district[region_size]) {
                            Rcpp::stop("redist_smc produced a plan with a %d seat district, "
                                       "which is not an allowed district size.\n%s",
                                       region_size, plan_ptr->debug_string(true));
                        }
                    }
                }
            }
        }

        // if not all diagnostic mode and split district only reorder the plans
        if (!diagnostic_mode &&
            splitting_size_regime != SplittingSizeScheduleType::DistrictOnlySMD) {
            reorder_all_plans(pool, plan_ensemble_ptr->plan_ptr_vec,
                              dummy_plan_ensemble_ptr->plan_ptr_vec);
        }

        // Resample while the dummy ensemble is still in scope to gather into,
        // unless that was already done before the final mergesplit rounds.
        // After this the plans are in their final order, so everything below
        // exports them to R exactly once.
        if (final_resample && !resample_before_final_ms) {
            do_final_resample();
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

    // Export the per-region matrices first, then drop everything in the
    // ensemble except the plan region ids. R stores the plans 4 bytes per unit
    // where the sampler stores them in 1, so the plan matrix is the single
    // largest allocation of the whole run and it is worth making it with as
    // little else resident as possible.
    Rcpp::IntegerMatrix seats_mat = plan_sizes_saved
                                        ? plan_ensemble_ptr->get_R_sizes_matrix(pool)
                                        : Rcpp::IntegerMatrix(1, 1);
    Rcpp::IntegerMatrix region_pops_mat = plan_ensemble_ptr->get_region_pops_matrix(pool);
    // The forests and linking edges have to be exported before the ensemble
    // releases them. This is after any final resampling so they match the plans.
    // RObject keeps these protected from R's GC while the plans matrix is made
    Rcpp::RObject spanning_forests;
    Rcpp::RObject linking_edges;
    if (return_augmented_samples) {
        spanning_forests = plan_ensemble_ptr->get_R_forest_matrix(map_params.num_edges, pool);
        if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
            linking_edges = plan_ensemble_ptr->get_R_linking_edges(map_params.graph_edge_index);
        }
    }
    plan_ensemble_ptr->release_all_but_plan_ids();

    Rcpp::IntegerMatrix plans_mat = plan_ensemble_ptr->get_R_plans_matrix();
    // the ids now live in R's memory, so hand the sampler's copy back
    std::vector<RegionID>().swap(plan_ensemble_ptr->flattened_all_plans);

    // Return results
    Rcpp::List out = Rcpp::List::create(
        Rcpp::_["plans_mat"] = plans_mat, // integer matrix to store final plans
        Rcpp::_["seats"] = seats_mat,     // saves sizes matrix if needed
        Rcpp::_["region_pops"] = region_pops_mat,
        Rcpp::_["plan_seats_saved"] = plan_sizes_saved, Rcpp::_["log_weights"] = log_weights,
        Rcpp::_["step_types"] = step_types,
        Rcpp::_["merge_split_steps"] = merge_split_step_vec,
        Rcpp::_["log_blank_map_target_density"] = log_blank_map_target_density,
        Rcpp::_["multidistrict_selection_alpha"] = multidistrict_selection_alpha
    );
    out["spanning_forests"] = spanning_forests;
    out["linking_edges"] = linking_edges;

    // add all the diagnostics
    smc_diagnostics.add_diagnostics_to_out_list(out);

    if constexpr (DEBUG_GSMC_PLANS_VERBOSE)
        Rprintf("Returning to R!\n");
    return out;
}
