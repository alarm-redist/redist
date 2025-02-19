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
#include "gredist_types.h"
#include "splitting.h"
#include "weights.h"


int run_merge_split_step_on_a_plan(
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    bool split_district_only, std::string const merge_prob_type, 
    Plan &plan, Plan &new_plan, 
    TreeSplitter &tree_splitter,
    int const nsteps_to_run
);

void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec, 
    std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
    bool const split_district_only, std::string const merge_prob_type, 
    int const nsteps_to_run,
    Rcpp::IntegerMatrix::Column success_count_vec
);

#endif
