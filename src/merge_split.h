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
#include "graph_plan_type.h"
#include "forest_plan_type.h"
#include "linking_edge_plan_type.h"
#include "ust_sampler.h"
#include "merging.h"

/*
 * Main entry point.
 *
 * USING MCMMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
// [[Rcpp::export]]
Rcpp::List ms_plans(
    int nsims, int warmup, int thin, 
    int const ndists, List const &adj_list,
    const arma::uvec &counties, const arma::uvec &pop,
    double const target, double const lower, double const upper,
    double rho, // compactness 
    arma::umat region_id_mat, arma::umat region_sizes_mat,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &merge_prob_type, // method for setting probability of picking a pair to merge
    List const &control, // control has pop temper, and k parameter value, and whether only district splits are allowed
    List const &constraints, // constraints 
    int verbosity = 3, bool diagnostic_mode = false
);


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
int adapt_ms_parameters(const Graph &g, int n_distr, double thresh,
                         double tol, arma::subview_col<arma::uword> const &plan, const uvec &counties,
                         Multigraph const &cg, const uvec &pop, double target,
                         RNGState &rng_state);

#endif