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
#include "gredist_types.h"
#include "base_plan_type.h"
#include "tree_splitter_types.h"





/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(const SplittingSchedule &splitting_schedule,
        const Graph &g, int &k, int const last_k, 
                      const std::vector<double> &unnormalized_weights, double thresh,
                      double tol, std::vector<std::unique_ptr<Plan>> const &plan_ptrs_vec, 
                      const uvec &counties,
                      Multigraph &cg, const uvec &pop, 
                      int const min_region_cut_size, int const max_region_cut_size,
                      bool split_district_only,
                      double const target, int const verbosity);




#endif
