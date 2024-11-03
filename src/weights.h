#pragma once
#ifndef WEIGHTS_H
#define WEIGHTS_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"


#include <string>
#include <random>
#include <unordered_map>
#include <set>
#include <utility> // for std::pair
#include <cmath>
#include "redist_types.h"
#include "tree_op.h"


//' Computes the effective sample size from log incremental weights
//'
//' Takes a vector of log incremental weights and computes the effective sample
//' size which is the sum of the weights squared divided by the sum of squared
//' weights
//'
//'
//' @title Compute Effective Sample Size
//'
//' @param log_wgt vector of log incremental weights
//'
//' @details No modifications to inputs made
//'
//' @return sum of weights squared over sum of squared weights (sum(wgt)^2 / sum(wgt^2))
//'
double compute_n_eff(const std::vector<double> &log_wgt);


void get_all_adj_pairs(
    Graph const &g, std::vector<std::pair<int, int>> &adj_pairs_vec,
    std::vector<int> const &vertex_region_ids,
    std::vector<bool> const valid_regions
);


double get_log_mh_ratio(
    const Graph &g, 
    const int region1_id, const int region2_id,
    const std::vector<int> &old_vertex_region_ids,
    const std::vector<int> &new_vertex_region_ids,
    const int num_old_adj_regions, const int num_new_adj_regions
);

//' Computes log unnormalized weights for vector of plans
//'
//' Using the procedure outlined in <PAPER HERE> this function computes the log
//' incremental weights and the unnormalized weights for a vector of plans (which
//' may or may not be the same depending on the parameters).
//'
//' @title Compute Log Unnormalized Weights
//'
//' @param pool A threadpool for multithreading
//' @param g A graph (adjacency list) passed by reference
//' @param plans_vec A vector of plans to compute the log unnormalized weights
//' of
//' @param split_district_only whether or not to compute the weights under 
//' the district only split scheme or not. If `split_district_only` is true
//' then uses optimal weights from one-district split scheme.
//' @param log_incremental_weights A vector of the log incremental weights
//' computed for the plans. The value of `log_incremental_weights[i]` is
//' the log incremental weight for `plans_vec[i]`
//' @param unnormalized_sampling_weights A vector of the unnormalized sampling
//' weights to be used with sampling the `plans_vec` in the next iteration of the
//' algorithm. Depending on the other hyperparameters this may or may not be the
//' same as `exp(log_incremental_weights)`
//' @param target Target population of a single district
//' @param pop_temper <DETAILS NEEDED>
//'
//' @details Modifications
//'    - The `log_incremental_weights` is updated to contain the incremental
//'    weights of the plans
//'    - The `unnormalized_sampling_weights` is updated to contain the unnormalized
//'    sampling weights of the plans for the next round
void get_all_plans_log_gsmc_weights(
        RcppThread::ThreadPool &pool,
        const Graph &g, std::vector<Plan> &plans_vec,
        bool split_district_only,
        std::vector<double> &log_incremental_weights,
        std::vector<double> &unnormalized_sampling_weights,
        double target, double pop_temper
);


#endif
