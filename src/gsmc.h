#pragma once
#ifndef GENERALIZED_SMC_H
#define GENERALIZED_SMC_H

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
#include "weights.h"
#include "splitting.h"



//' Uses gsmc method to generate a sample of `M` plans in `c++`
//'
//' Using the procedure outlined in <PAPER HERE> this function uses Sequential
//' Monte Carlo (SMC) methods to generate a sample of `M` plans
//'
//' @title Run redist gsmc
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
List gsmc_plans(
        int N, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        double target, double lower, double upper,
        int M, // Number of particles aka number of different plans
        List control,
        int ncores = -1, int verbosity = 3, bool diagnostic_mode = false);




bool OLD_attempt_region_split(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                           Plan &plan, const int region_id_to_split,
                           std::vector<int> &new_region_ids,
                           std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                           double &lower, double upper, double target, int k_param);

#endif
