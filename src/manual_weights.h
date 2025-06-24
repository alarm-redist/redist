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
#include "weights.h"
#include "redist_alg_helpers.h"



// [[Rcpp::export]]
Rcpp::NumericVector compute_log_unnormalized_plan_target_density(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double const pop_temper,  double const rho,
    int const ndists, int const total_seats, int const num_regions,
    Rcpp::IntegerVector const &district_seat_sizes,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int const num_threads
);


// [[Rcpp::export]]
Rcpp::NumericMatrix compute_log_unnormalized_region_target_density(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    List const &constraints, double const pop_temper,  double const rho,
    int const ndists, int const total_seats, int const num_regions,
    Rcpp::IntegerVector const &district_seat_sizes,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int const num_threads
);


// [[Rcpp::export]]
arma::vec compute_plans_log_optimal_weights(
    List const &adj_list, arma::uvec const &counties, arma::uvec const &pop,
    List const &constraints, double const pop_temper,  double const rho,
    std::string const &splitting_schedule_str,
    int const ndists, int const total_seats, Rcpp::IntegerVector const &district_seat_sizes,
    int const num_regions,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int num_threads
);

#endif