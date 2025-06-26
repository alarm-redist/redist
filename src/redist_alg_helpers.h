#pragma once
#ifndef SMC_ALG_HELPERS_H
#define SMC_ALG_HELPERS_H

// [[Rcpp::depends(redistmetrics)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <limits>
#include <RcppThread.h>

#include "base_plan_type.h"
#include "graph_plan_type.h"
#include "forest_plan_type.h"
#include "linking_edge_plan_type.h"
#include "redist_types.h"
#include "tree_splitter_types.h"


// [[Rcpp::export]]
Rcpp::List maximum_input_sizes();



//' Get canonically relabeled plans matrix
//'
//' Given a matrix of 1-indexed plans (or partial plans) this function 
//' returns a new plans matrix with all the plans labeled canonically. 
//' The canonical labelling of a plan is the one where the region of the 
//' first vertex gets mapped to 1, the region of the next smallest vertex
//' in a different region than the first gets mapped to 2, and so on. This
//' is guaranteed to result in the same labelling for any plan where the 
//' region ids have been permuted. 
//'
//'
//' @param plans_mat A matrix of 1-indexed plans
//' @param num_regions The number of regions in the plan
//' @param num_threads The number of threads to use. Defaults to number of machine threads.
//'
//' @details Modifications
//'    - None
//'
//' @returns A matrix of canonically labelled plans
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_canonical_plan_labelling(
    Rcpp::IntegerMatrix const &plans_mat,
    int const num_regions,
    int const num_threads
);


/*
 * Creates a Rcpp Threadpool
 *
 * 
 */
RcppThread::ThreadPool get_thread_pool(int const num_threads);


//' Reorders all the plans in the vector by order a region was split
//'
//' Takes a vector of plans and uses the vector of dummy plans to reorder
//' each of the plans by the order a region was split.
//'
//'
//' @title Reorders all the plans in the vector by order a region was split
//'
//' @param pool A threadpool for multithreading
//' @param plan_ptrs_vec A vector of pointers to plans 
//' @param dummy_plans_vec A vector of pointers to dummy plans 
//'
//' @details Modifications
//'    - Each plan in the `plans_vec` object is reordered by when the region was split
//'    - Each plan is a shallow copy of the plans in `plans_vec`
//'
//' @noRd
//' @keywords internal
void reorder_all_plans(
    RcppThread::ThreadPool &pool,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec);


std::unique_ptr<TreeSplitter> get_tree_splitters(
    MapParams const &map_params,
    SplittingMethodType const splitting_method,
    Rcpp::List const &control,
    int const nsims
);

// lightweight container for plans 
class PlanEnsemble {

    public:
        // constructor for empty plans
        PlanEnsemble(
            MapParams const &map_params,
            int const total_pop, int const nsims, 
            SamplingSpace const sampling_space,
            RcppThread::ThreadPool &pool,
            int const verbosity = 3
        );
        // constructor for non-empty starting plans 
        PlanEnsemble(
            MapParams const &map_params,
            int const num_regions, int const nsims, 
            SamplingSpace const sampling_space,
            Rcpp::IntegerMatrix const &plans_mat, 
            Rcpp::IntegerMatrix const &region_sizes_mat,
            std::vector<RNGState> &rng_states,
            RcppThread::ThreadPool &pool,
            int const verbosity = 3 
        );

    
        int const nsims;
        int const V;
        int const ndists;
        int const total_seats;
        std::vector<RegionID> flattened_all_plans;
        std::vector<RegionID> flattened_all_region_sizes;
        std::vector<int> flattened_all_region_pops;
        std::vector<int> flattened_all_region_order_added;
        std::vector<std::unique_ptr<Plan>> plan_ptr_vec;
        
        
        // exports current plans to 1-indexed Rcpp matrix 
        Rcpp::IntegerMatrix get_R_plans_matrix();
        // export current region sizes to Rcpp matrix 
        Rcpp::IntegerMatrix get_R_sizes_matrix(
            RcppThread::ThreadPool &pool
        );
    
};

PlanEnsemble get_plan_ensemble(
    MapParams const &map_params,
    int const num_regions, int const nsims,
    SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, 
    Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states,
    RcppThread::ThreadPool &pool,
    int const verbosity
);


std::unique_ptr<PlanEnsemble> get_plan_ensemble_ptr(
    MapParams const &map_params,
    int const num_regions, int const nsims, 
    SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, 
    Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states,
    RcppThread::ThreadPool &pool,
    int const verbosity
);

// swaps the contents of two plan ensembles 
void swap_plan_ensembles(
    PlanEnsemble &plan_ensemble1,
    PlanEnsemble &plan_ensemble2
);

// Wrapper object for all non-essential SMC diagnostics 
class SMCDiagnostics{

    public:
        SMCDiagnostics(
            SamplingSpace const sampling_space, SplittingMethodType const splitting_method_type,
            SplittingSizeScheduleType const splitting_schedule_type, 
            std::vector<bool> const &merge_split_step_vec,
            int const V, int const nsims,
            int const ndists, int const total_seats, int const initial_num_regions,
            int const total_smc_steps, int const total_ms_steps,
            int const diagnostic_level,
            bool const splitting_all_the_way, bool const split_district_only
        );

    
    int const diagnostic_level;
    int const total_steps;
    // Level 0
    // Essential Diagnostics that are always created 
    std::vector<double> log_wgt_stddevs; // log weight std devs
    std::vector<double> acceptance_rates; // Tracks the acceptance rate - total number of tries over nsims - for each round
    std::vector<int> nunique_parents; // number of unique parents
    std::vector<double> n_eff; // Tracks the effective sample size for the weights of each round
    // For each merge split step this counts the number of attempts that were made
    std::vector<int> num_merge_split_attempts_vec;
    // Only required for graph sampling 
    std::vector<int> cut_k_values; // k value used at each step

    
    // Level 1
    // These are all nsims by number of smc steps 
    arma::dmat log_incremental_weights_mat; // entry [i][s] is the log unnormalized weight of particle i AFTER split s
    Rcpp::IntegerMatrix draw_tries_mat; // Entry [i][s] is the number of tries it took to form particle i on split s
    Rcpp::IntegerMatrix parent_index_mat; // Entry [i][s] is the index of the parent of particle i at split s
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of 
    // successful merge splits performed for plan i on merge split round s
    Rcpp::IntegerMatrix merge_split_successes_mat;
    // counts the size of the trees
    Rcpp::IntegerMatrix tree_sizes_mat; // ndists by total_steps matrix
    Rcpp::IntegerMatrix successful_tree_sizes_mat; // ndists by total_steps matrix


    // Level 2
    Rcpp::IntegerMatrix parent_unsuccessful_tries_mat;

    // level 3
    std::vector<Rcpp::IntegerMatrix> all_steps_plan_region_ids_list;
    std::vector<std::vector<VertexGraph>> all_steps_forests_adj_list;
    std::vector<std::vector<std::vector<std::array<double, 3>>>> all_steps_linking_edge_list;
    std::vector<std::vector<int>> all_steps_valid_region_sizes_to_split;
    std::vector<std::vector<int>> all_steps_valid_split_region_sizes;
    std::vector<Rcpp::IntegerMatrix> region_sizes_mat_list;

    // adds full diagnostics (takes a lot of memory)
    void add_full_step_diagnostics(
        int const total_steps, bool const splitting_all_the_way,
        int const step_num, int const merge_split_step_num, int const smc_step_num,
        bool const is_smc_step,
        SamplingSpace const sampling_space,
        RcppThread::ThreadPool &pool,
        PlanEnsemble &plan_ensemble,
        PlanEnsemble &new_plans_ensemble,
        SplittingSchedule const &splitting_schedule
    );

    // Updates the out list with all the diagnostics 
    void add_diagnostics_to_out_list(Rcpp::List &out);
};


// Resamples plans in place using weights 
// [[Rcpp::export]]
Rcpp::IntegerVector resample_plans_lowvar(
    Rcpp::NumericVector const &normalized_weights,
    Rcpp::IntegerMatrix &plans_mat,
    Rcpp::IntegerMatrix &region_sizes_mat,
    bool const reorder_sizes_mat
);

#endif
