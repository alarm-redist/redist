#pragma once
#ifndef SPLITTING_INSPECTION_H
#define SPLITTING_INSPECTION_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "active_dev.h"
#include "smc_and_mcmc.h"


// [[Rcpp::export]]
List perform_a_valid_region_split(
        List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int k_param, int region_id_to_split,
        double target, double lower, double upper,
        int N, int num_regions, int num_districts,
        std::vector<int> region_ids, std::vector<int> region_dvals,
        std::vector<double> region_pops,
        bool split_district_only, bool verbose
);

#endif
