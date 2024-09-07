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

double compute_log_incremental_weight(const Graph &g, const Plan &plan);

Graph get_region_graph(const Graph &g, const Plan &plan);

#endif
