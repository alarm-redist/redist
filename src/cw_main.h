#ifndef CW_MAIN_H
#define CW_MAIN_H

#include "smc_base.h"
#include "tree_op.h"

#include <cli/progress.h>
#include <string>

/*
 * Cyclewalk MCMC sampler. Mirrors the modular shape used by ms_plans /
 * run_redist_smc: a single chain over a PlanEnsemble of size 1, MapParams
 * for graph + bounds + MMD bookkeeping, ScoringFunction for constraints,
 * USTSampler for initial-tree construction, and a per-chain RNGState.
 *
 * Returns a list with at minimum: plans, mhdecisions, diagnostics. For MMD
 * also returns plan_sizes (ndists x nsims) with seat counts.
 */
Rcpp::List
cyclewalk_plans(int N, int warmup, int thin,
                int ndists, int total_seats,
                Rcpp::IntegerVector const &district_seat_sizes,
                Rcpp::List const &adj_list,
                arma::uvec const &counties, arma::uvec const &pop,
                double target, double lower, double upper, double compactness,
                Rcpp::IntegerMatrix const &init_plan,
                Rcpp::IntegerMatrix const &init_seats,
                Rcpp::List const &control, Rcpp::List const &constraints,
                Rcpp::List const &edge_weights,
                int instep, double cycle_walk_frac,
                int verbosity, bool diagnostic_mode);

#endif
