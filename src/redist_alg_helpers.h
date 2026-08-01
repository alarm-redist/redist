#pragma once
#ifndef SMC_ALG_HELPERS_H
#define SMC_ALG_HELPERS_H


#include <memory>
#include <string_view>
#include <vector>
#include <RcppArmadillo.h>

#include "advanced_types.h"
#include "base_plan_type.h"



class TreeSplitter;
class SplittingSchedule;
class RNGState;

namespace RcppThread {
    class ThreadPool;
}


Rcpp::List maximum_input_sizes();

/*
 * Convert zero-indxed R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const Rcpp::List &l);

/*
 * Creates a reindexing vector for the plan with two regions merged
 * does this by making region2_id map to region1_id. Nothing else
 * changes
 */
void set_merged_region_reindex_vec(int const num_regions, std::vector<int> &region_reindex_vec,
                                   int const region1_id, int const region2_id);

/*
 *  Reorders all the plans in the vector by order a region was split
 *
 *  Takes a vector of plans and uses the vector of dummy plans to reorder
 *  each of the plans by the order a region was split.
 *
 *
 *  @title Reorders all the plans in the vector by order a region was split
 *
 *  @param pool A threadpool for multithreading
 *  @param plan_ptrs_vec A vector of pointers to plans
 *  @param dummy_plans_vec A vector of pointers to dummy plans
 *
 *  @details Modifications
 *     - Each plan in the `plans_vec` object is reordered by when the region was split
 *     - Each plan is a shallow copy of the plans in `plans_vec`
 *
 */
void reorder_all_plans(RcppThread::ThreadPool &pool,
                       std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
                       std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec);

std::vector<std::unique_ptr<TreeSplitter>>
get_tree_splitter_ptrs(MapParams const &map_params, SplittingMethodType const splitting_method,
                       SamplingSpace const sampling_space,
                       Rcpp::List const &control, int const nsims, int const num_threads);

// lightweight container for plans
class PlanEnsemble {

  public:
    // constructor for empty plans
    PlanEnsemble(MapParams const &map_params, int const total_pop, int const nsims,
                 SamplingSpace const sampling_space, RcppThread::ThreadPool &pool,
                 int const verbosity = 3);
    // constructor for non-empty starting plans
    PlanEnsemble(MapParams const &map_params, SplittingSchedule const &splitting_schedule,
                 int const num_regions, int const nsims, SamplingSpace const sampling_space,
                 Rcpp::IntegerMatrix const &plans_mat,
                 Rcpp::IntegerMatrix const &region_sizes_mat, std::vector<RNGState> &rng_states,
                 RcppThread::ThreadPool &pool, int const verbosity = 3);

    int const nsims;
    int const V;
    int const ndists;
    int const total_seats;
    SamplingSpace const sampling_space;
    std::vector<RegionID> flattened_all_plans;
    std::vector<RegionID> flattened_all_region_sizes;
    std::vector<int> flattened_all_region_pops;
    std::vector<int> flattened_all_region_order_added;
    // Empty unless sampling space is ForestSpace or LinkingEdgeSpace.
    int const num_forest_edge_bit_words_per_plan;
    std::vector<EdgeBitWord> flattened_all_forest_edge_bits;
    std::vector<std::unique_ptr<Plan>> plan_ptr_vec;

    // exports current plans to 1-indexed Rcpp matrix
    Rcpp::IntegerMatrix get_R_plans_matrix();
    // export current region sizes to Rcpp matrix
    Rcpp::IntegerMatrix get_R_sizes_matrix(RcppThread::ThreadPool &pool);
    // export current region populations to Rcpp matrix
    Rcpp::IntegerMatrix get_region_pops_matrix(RcppThread::ThreadPool &pool);
    // counts the number of unique plans in the ensemble
    int count_unique_plans(RcppThread::ThreadPool &pool) const;

    // debugging methods
    // checks all plans are valid. 
    void check_all_plans_valid(
        MapParams const &map_params,
        std::string_view where
    );
};

PlanEnsemble get_plan_ensemble(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity);

std::unique_ptr<PlanEnsemble> get_plan_ensemble_ptr(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity);


// converts trees to compact edge list form 
// This only supports undirected trees 
std::vector<bool> vector_tree_to_edge_vector(
    GraphEdgeIndex const &edge_index,
    Tree const &tree
);

// Converts a graph edge index to an R deciperable list
// where its a list of length edge_index.num_edges
// and each element is the pair (u,v) of the vertices in the 
// edge it represents 
Rcpp::List graph_edge_index_to_list(
    GraphEdgeIndex const &edge_index
);

// TODO: This should be moved to smc.cpp
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
                   bool const split_district_only);

    int const diagnostic_level;
    int const total_steps;
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
    arma::dmat log_incremental_weights_mat; // entry [i][s] is the log unnormalized weight of
                                            // particle i AFTER split s
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
    std::vector<std::vector<Tree>> all_steps_forests_adj_list;
    std::vector<std::vector<std::vector<std::array<double, 3>>>> all_steps_linking_edge_list;
    std::vector<std::vector<int>> all_steps_valid_region_sizes_to_split;
    std::vector<std::vector<int>> all_steps_valid_split_region_sizes;
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


    // adds full diagnostics (takes a lot of memory)
    void add_full_step_diagnostics(int const total_steps, bool const splitting_all_the_way,
                                   int const step_num, int const merge_split_step_num,
                                   int const smc_step_num, bool const is_smc_step,
                                   SamplingSpace const sampling_space,
                                   RcppThread::ThreadPool &pool, PlanEnsemble &plan_ensemble,
                                   PlanEnsemble &new_plans_ensemble,
                                   SplittingSchedule const &splitting_schedule);

    // Updates the out list with all the diagnostics
    void add_diagnostics_to_out_list(Rcpp::List &out);
};


#endif
