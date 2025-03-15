#pragma once
#ifndef MANUAL_SPLITTING_H
#define MANUAL_SPLITTING_H


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
#include "splitting.h"
#include "merging.h"
#include "graph_plan_type.h"
#include "weights.h"



// [[Rcpp::export]]
double compute_a_log_optimal_weight(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List control, 
    int ndists, int num_regions, 
    double lower, double target, double upper,
    arma::umat region_ids, arma::umat region_sizes
);


// [[Rcpp::export]]
arma::vec compute_log_unnormalized_plan_target_density(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double pop_temper,  double rho,
    int ndists, int num_regions,
    double lower, double target, double upper,
    arma::umat region_ids, arma::umat region_sizes,
    int num_threads
);

#endif