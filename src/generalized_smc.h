#pragma once
#ifndef GENERALIZED_SMC_H
#define GENERALIZED_SMC_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <format>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "redist_types.h"
#include "generalized_smc_helpers.h"


//' @export
// [[Rcpp::export]]
List generalized_smc_plans(
        int N, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        double target, double lower, double upper,
        int M, int k_param,// Number of particles aka number of different plans
        List control,
        int ncores = -1, int verbosity = 3);




bool attempt_region_split(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                           Plan &plan, const std::string region_to_split,
                           std::vector<std::string> &new_region_labels,
                           std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                           double &lower, double upper, double target, int k_param);

#endif
