/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Main entry point and MCMC loop
 *
 * This file contains the main entry point for the CycleWalk algorithm,
 * which is called from R via Rcpp.
 ********************************************************/

#ifndef CW_MAIN_H
#define CW_MAIN_H

#include "smc_base.h"
#include "tree_op.h"

#include <string>
#include <cli/progress.h>

/*
 * Main entry point for CycleWalk MCMC sampler.
 *
 * Sample `N` redistricting plans on map `g`, ensuring that the
 * population deviation is between `lower` and `upper` (with target `target`).
 *
 * Returns a List with:
 *   - plans: umat of district assignments (V x n_samples)
 *   - mhdecisions: vector of MH accept/reject decisions
 *   - (future: additional diagnostics)
 *
 * Note: Rcpp::export is in cw_main.cpp, not here, to avoid duplicate entries
 */
Rcpp::List cyclewalk_plans(
    int N,                          // Number of MCMC iterations
    Rcpp::List l,                   // Adjacency list
    const arma::uvec init,          // Initial plan
    const arma::uvec& counties,     // County assignments
    const arma::uvec& pop,          // Population vector
    int n_distr,                    // Number of districts
    double target,                  // Target population per district
    double lower,                   // Lower population bound
    double upper,                   // Upper population bound
    double compactness,             // Compactness parameter (0=uniform over partitions)
    Rcpp::List constraints,         // Constraint specifications
    Rcpp::List control,             // Control parameters
    Rcpp::List edge_weights,        // Edge weights (list of list(edge=c(u,v), weight=w))
    int thin,                       // Thinning interval
    int instep,                     // MCMC iterations per recorded sample
    double cycle_walk_frac,         // Fraction of cycle walk vs forest walk
    int verbosity                   // Verbosity level (0=silent, 1=normal, 3=verbose)
);

#endif // CW_MAIN_H
