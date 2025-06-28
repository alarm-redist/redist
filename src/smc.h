#pragma once
#ifndef SMC_H
#define SMC_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <iostream>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>
#include <atomic>
#include <chrono>


#include "splitting.h"
#include "merging.h"
#include "weights.h"
#include "redist_alg_helpers.h"
#include "tree_splitter_types.h"
#include "scoring.h"
#include "ust_sampler.h"



//' Run Optimalgsmc with Merge Split
//'
//' Uses gsmc method with optimal weights and merge split steps to generate a sample of `M` plans in `c++` 
//' 
//' 
//' Using the procedure outlined in <PAPER HERE> this function uses Sequential
//' Monte Carlo (SMC) methods to generate a sample of `M` plans
//'
//'
//' @param ndists The number of districts the final plans will have
//' @param adj_list A 0-indexed adjacency list representing the undirected graph
//' which represents the underlying map the plans are to be drawn on
//' @param counties Vector of county labels of each vertex in `g`
//' @param pop A vector of the population associated with each vertex in `g`
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param nsims The number of plans (samples) to draw
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param control Named list of additional parameters.
//' @param num_threads The number of threads the threadpool should use
//' @param verbosity What level of detail to print out while the algorithm is
//' running <ADD OPTIONS>
//' @export
//' @keywords internal
// [[Rcpp::export]]
List run_redist_gsmc(
        int const nsims, 
        int const total_seats, int const ndists, Rcpp::IntegerVector const district_seat_sizes,
        int const initial_num_regions, 
        List const &adj_list,
        arma::uvec const &counties, const arma::uvec &pop,
        Rcpp::CharacterVector const &step_types,
        double const target, double const lower, double const upper,
        double const rho, // compactness 
        std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
        List const &control, // control has pop temper, and k parameter value, and splitting method are allowed
        List const &constraints, // constraints 
        int const verbosity, int const diagnostic_level,
        Rcpp::IntegerMatrix const &region_id_mat, 
        Rcpp::IntegerMatrix const &region_sizes_mat
);



#endif
