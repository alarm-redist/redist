#ifndef MERGESPLIT_H
#define MERGESPLIT_H

#include "smc_base.h"

#include <cli/progress.h>
#include <cmath>
#include <string>

// [[Rcpp::depends(redistmetrics)]]

#include "base_plan_type.h"
#include "map_calc.h"
#include "mcmc_gibbs.h"
#include "merging.h"
#include "redist_alg_helpers.h"
#include "scoring.h"
#include "splitting_schedule_types.h"
#include "tree_op.h"
#include "ust_sampler.h"
#include "wilson.h"
#include <kirchhoff_inline.h>

/*
 * Main entry point.
 *
 * USING MCMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
Rcpp::List ms_plans(
    int const nsims, int const warmup, int const thin, int const ndists, int const total_seats,
    Rcpp::IntegerVector const &district_seat_sizes, List const &adj_list,
    const arma::uvec &counties, const arma::uvec &pop, double const target, double const lower,
    double const upper,
    double const rho, // compactness
    Rcpp::IntegerMatrix const &init_plan, Rcpp::IntegerMatrix const &init_seats,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &pair_rule, // method for setting probability of picking a pair to merge
    List const &control,     // control has pop temper, and k parameter value, and whether only
                             // district splits are allowed
    List const &constraints, // constraints
    int const verbosity = 3, bool const diagnostic_mode = false);

#endif
