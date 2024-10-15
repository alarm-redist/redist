#pragma once
#ifndef SMC_AND_MCMC_H
#define SMC_AND_MCMC_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "redist_types.h"
#include "gsmc_helpers.h"

bool get_edge_to_cut(Tree &ust, int root,
                     int k_param, bool split_district_only,
                     const uvec &pop, const Plan &plan, const int region_id_to_split,
                     const double lower, const double upper, const double target,
                     int &new_region1_tree_root, int &new_region1_dval, double &new_region1_pop,
                     int &new_region2_tree_root, int &new_region2_dval, double &new_region2_pop
);


#endif
