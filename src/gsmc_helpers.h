#pragma once
#ifndef GENERALIZED_SMC_HELPERS_H
#define GENERALIZED_SMC_HELPERS_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"


#include <string>
#include <random>
#include <unordered_map>
#include <set>
#include <utility> // for std::pair
#include <cmath>
#include "redist_types.h"

double compute_n_eff(const std::vector<double> &log_wgt);

double choose_multidistrict_to_split(
        Plan const&plan, int &region_id_to_split);

//' Compute the log incremental weight of a plan
//'
//' Given a plan object this computes the minimum variance weights as derived in
//' <PAPER NAME HERE>. This is equivalent to the inverse of a sum over all
//' adjacent regions in a plan.
//'
//' @title Compute Incremental Weight of a plan
//'
//' @param plan A plan object
//' @param g The underlying map graph
//' @param target The target population for a single district
//' @param pop_temper The population tempering parameter
//'
//' @details No modifications to inputs made
//'
//' @return the log of the incremental weight of the plan
//'
double compute_log_incremental_weight(
     const Graph &g, const Plan &plan,
     const double target, const double pop_temper);

Graph get_region_graph(const Graph &g, const Plan &plan);

#endif
