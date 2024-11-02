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


int run_merge_split_step( 
    Graph const &g, const uvec &counties, Multigraph &cg, const uvec &pop,
    bool split_district_only,
    Tree &ust, int k_param,
    Plan &plan, int nsteps_to_run,
    double &lower, double upper, double target,
    std::vector<bool> &visited, std::vector<bool> &ignore
);

#endif