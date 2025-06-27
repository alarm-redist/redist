#pragma once
#ifndef SPLITTING_H
#define SPLITTING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <iostream>
#include <functional>
#include <algorithm>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "redist_types.h"
#include "base_plan_type.h"
#include "tree_splitter_types.h"
#include "splitting_schedule_types.h"





/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(
    const MapParams &map_params, const SplittingSchedule &splitting_schedule,
    RNGState &rng_state,
    int &k, int const last_k, 
    const arma::vec &unnormalized_weights, double thresh,
    double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
    bool split_district_only,
    int const verbosity);




#endif
