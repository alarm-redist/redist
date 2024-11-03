#pragma once
#ifndef MERGING_H
#define MERGING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <set>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "redist_types.h"
#include "splitting.h"
#include "weights.h"


int run_merge_split_step_on_a_plan( 
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    bool const split_district_only,
    int const k_param,
    Plan &plan, int const nsteps_to_run,
    double const lower, double const upper, double const target
);

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    std::vector<Plan> &plans_vec, 
    bool const split_district_only, int const k_param,
    int const nsteps_to_run,
    double const lower, double const upper, double const target,
    std::vector<int> &success_count_vec
);

#endif