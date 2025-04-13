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
#include "gredist_types.h"
#include "tree_splitter_types.h"


//' Copies data from a vector of `Plan` objects into an Rcpp Matrix
//'
//' Takes a vector of plans and copies all the data into an RcppMatrix
//' of the same size using the Rcpp Threadpool to copy in parallel. 
//'
//'
//' @title Copies data from an arma Matrix into an Rcpp Matrix
//'
//' @param pool A threadpool for multithreading
//' @param arma_mat Subview of an arma unsigned integer matrix 
//' @param rcpp_mat A matrix of integers with the same size as the arma_mat
//'
//' @details Modifications
//'    - The `rcpp_mat` is filled in with the data om the arma matrix subview
//'
//' @noRd
//' @keywords internal
void copy_plans_to_rcpp_mat(
    RcppThread::ThreadPool &pool,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
    Rcpp::IntegerMatrix &rcpp_mat,
    bool const copy_sizes_not_ids
);



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


// Wrapper object for all non-essential diagnostics 
class SMCDiagnostics{

    public:
        SMCDiagnostics(
            SamplingSpace sampling_space, SplittingMethodType splitting_method_type,
            SplittingSizeScheduleType splitting_schedule_type, 
            std::vector<bool> &merge_split_step_vec,
            int V, int nsims,
            int ndists, int initial_num_regions,
            int total_smc_steps, int total_ms_steps,
            int diagnostic_level,
            bool splitting_all_the_way, bool split_district_only
        );

    
    int diagnostic_level;
    int total_steps;
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
    std::vector<std::vector<Graph>> all_steps_forests_adj_list;
    std::vector<std::vector<std::vector<std::array<double, 3>>>> all_steps_linking_edge_list;
    std::vector<std::vector<int>> all_steps_valid_region_sizes_to_split;
    std::vector<std::vector<int>> all_steps_valid_split_region_sizes;
    std::vector<Rcpp::IntegerMatrix> region_sizes_mat_list;


    // 
    void add_full_step_diagnostics(
        int const total_steps, bool const splitting_all_the_way,
        int const step_num, int const merge_split_step_num, int const smc_step_num,
        bool const is_smc_step,
        SamplingSpace const sampling_space,
        RcppThread::ThreadPool &pool,
        std::vector<std::unique_ptr<Plan>>  &plans_ptr_vec, 
        std::vector<std::unique_ptr<Plan>>  &new_plans_ptr_vec,
        SplittingSchedule const &splitting_schedule
    );

    // Updates the out list with all the diagnostics 
    void add_diagnostics_to_out_list(Rcpp::List &out);
};



#endif
