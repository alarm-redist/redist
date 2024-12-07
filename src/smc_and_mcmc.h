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


#include "gredist_types.h"
#include "splitting.h"
#include "merging.h"
#include "weights.h"



//' Uses gsmc method with optimal weights and merge split steps to generate a sample of `M` plans in `c++`
//'
//' Using the procedure outlined in <PAPER HERE> this function uses Sequential
//' Monte Carlo (SMC) methods to generate a sample of `M` plans
//'
//' @title Run Optimalgsmc with Merge Split
//'
//' @param N The number of districts the final plans will have
//' @param adj_list A 0-indexed adjacency list representing the undirected graph
//' which represents the underlying map the plans are to be drawn on
//' @param counties Vector of county labels of each vertex in `g`
//' @param pop A vector of the population associated with each vertex in `g`
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param M The number of plans (samples) to draw
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param control Named list of additional parameters.
//' @param num_threads The number of threads the threadpool should use
//' @param verbosity What level of detail to print out while the algorithm is
//' running <ADD OPTIONS>
//' @export
// [[Rcpp::export]]
List optimal_gsmc_with_merge_split_plans(
        int N, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        double target, double lower, double upper,
        int M, // M is Number of particles aka number of different plans
        List control, // control has pop temper, and k parameter value, and whether only district splits are allowed
        int num_threads = -1, int verbosity = 3, bool diagnostic_mode = false
);



#endif
