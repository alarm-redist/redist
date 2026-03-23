/********************************************************
 * BUD (Balanced Up-Down Walk) MCMC Redistricting Sampler
 * Main entry point and MCMC loop
 ********************************************************/

#ifndef BUD_MAIN_H
#define BUD_MAIN_H

#include "smc_base.h"
#include "tree_op.h"

#include <string>
#include <cli/progress.h>

/*
 * Main entry point for BUD MCMC sampler.
 *
 * Sample `N` redistricting plans on map `g`, using the Balanced Up-Down Walk
 * algorithm. Maintains a spanning forest with district-level tree and marked
 * edges, proposing multi-district cycle-based changes.
 */
Rcpp::List bud_plans(
    int N,                          // Number of MCMC iterations
    Rcpp::List l,                   // Adjacency list
    const arma::uvec init,          // Initial plan
    const arma::uvec& counties,     // County assignments
    const arma::uvec& pop,          // Population vector
    int n_distr,                    // Number of districts
    double target,                  // Target population per district
    double lower,                   // Lower population bound
    double upper,                   // Upper population bound
    double compactness,             // Compactness parameter
    Rcpp::List constraints,         // Constraint specifications
    Rcpp::List control,             // Control parameters
    Rcpp::List edge_weights,        // Edge weights
    int thin,                       // Thinning interval
    int instep,                     // Number of MCMC steps per sample
    int verbosity                   // Verbosity level
);

#endif // BUD_MAIN_H
