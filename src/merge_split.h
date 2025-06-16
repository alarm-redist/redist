#ifndef MERGESPLIT_H
#define MERGESPLIT_H

#include "smc_base.h"

#include <string>
#include <cli/progress.h>
#include <cmath>

// [[Rcpp::depends(redistmetrics)]]

#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include <kirchhoff_inline.h>
#include "mcmc_gibbs.h"
#include "splitting_schedule_types.h"
#include "scoring.h"
#include "redist_alg_helpers.h"
#include "base_plan_type.h"
#include "ust_sampler.h"
#include "merging.h"

/*
 * Main entry point.
 *
 * USING MCMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
Rcpp::List ms_plans(
    int const nsims, int const warmup, int const thin, 
    int const ndists, int const total_seats, Rcpp::IntegerVector const &district_seat_sizes, 
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    double const target, double const lower, double const upper,
    double const rho, // compactness 
    Rcpp::IntegerMatrix const &initial_plan, Rcpp::IntegerMatrix const &initial_region_sizes,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &merge_prob_type, // method for setting probability of picking a pair to merge
    List const &control, // control has pop temper, and k parameter value, and whether only district splits are allowed
    List const &constraints, // constraints 
    int const verbosity = 3, bool const diagnostic_mode = false
);


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
int adapt_ms_parameters(const Graph &g, int n_distr, double thresh,
                         double tol, PlanVector const &region_ids, const uvec &counties,
                         Multigraph const &cg, const uvec &pop, double target,
                         RNGState &rng_state);

#endif