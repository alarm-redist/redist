/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Main entry point and MCMC loop
 ********************************************************/

#include "cw_main.h"
#include "cw_lct.h"
#include "cw_partition.h"
#include "cw_forest_walk.h"
#include "cw_proposal.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List cyclewalk_plans(
    int N,
    Rcpp::List l,
    const arma::uvec init,
    const arma::uvec& counties,
    const arma::uvec& pop,
    int n_distr,
    double target,
    double lower,
    double upper,
    double compactness,
    Rcpp::List constraints,
    Rcpp::List control,
    Rcpp::List edge_weights,
    int thin,
    int instep,
    double cycle_walk_frac,
    int verbosity
) {
    // Re-seed RNG for reproducibility
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // Convert adjacency list to Graph
    Graph g = list_to_graph(l);
    int V = g.size();

    // Calculate number of output samples
    int n_out = N / thin + 2;  // +2 for init and final

    // Initialize output matrix
    umat plans(V, n_out, fill::zeros);
    plans.col(0) = init;
    plans.col(1) = init;

    // Initialize MH decisions
    IntegerVector mh_decisions(N / thin + 1);

    // Initialize diagnostic vectors
    std::vector<double> diag_accept_prob;
    std::vector<int> diag_cycle_length;
    std::vector<int> diag_n_valid_cuts;
    diag_accept_prob.reserve(N);
    diag_cycle_length.reserve(N);
    diag_n_valid_cuts.reserve(N);

    // Progress bar setup
    RObject bar = R_NilValue;
    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << "CYCLEWALK MCMC\n";
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
    }

    // Initialize LCT Partition
    if (verbosity >= 1) {
        Rcout << "\nInitializing partition with spanning trees...\n";
    }
    LCTPartition partition(V, n_distr);
    int init_result = partition.init_from_plan(g, init, pop, counties, lower, upper);

    if (init_result != 0) {
        Rcpp::stop("Failed to initialize partition - Wilson's algorithm failed.");
    }

    // Set edge weights if provided
    if (edge_weights.size() > 0) {
        partition.set_edge_weights(edge_weights);
        if (verbosity >= 1) {
            Rcout << "Using " << edge_weights.size() << " custom edge weights.\n";
        }
    }

    if (verbosity >= 1) {
        partition.print_state(verbosity);
        Rcout << "\n";
        bar = cli_progress_bar(N, cli_config(false));
    }

    // Main MCMC loop
    int forest_walk_success = 0;
    int forest_walk_fail = 0;
    int cycle_walk_accept = 0;
    int cycle_walk_reject = 0;
    int cycle_walk_fail_no_adj = 0;
    int cycle_walk_fail_few_edges = 0;
    int cycle_walk_fail_no_path = 0;
    int cycle_walk_fail_no_cuts = 0;

    int idx = 1;
    try {
    for (int i = 1; i <= N; i++) {
        // Track whether any cycle_walk was accepted during instep iterations
        bool any_cycle_walk_accepted = false;

        // Run instep MCMC iterations per recorded sample
        for (int step = 0; step < instep; step++) {
            int result;

            if (r_unif() < cycle_walk_frac) {
                // Cycle walk proposal - can change districts
                CycleWalkDiagnostics diag;
                result = cycle_walk(partition, lower, upper, target, compactness, counties, constraints, diag);

                diag_accept_prob.push_back(diag.accept_prob);
                diag_cycle_length.push_back(diag.cycle_length);
                diag_n_valid_cuts.push_back(diag.n_valid_cuts);

                if (result == 1) {
                    cycle_walk_accept++;
                    any_cycle_walk_accepted = true;
                } else if (result == 0) {
                    cycle_walk_reject++;
                } else if (result == -1) {
                    cycle_walk_fail_no_adj++;
                } else if (result == -2) {
                    cycle_walk_fail_few_edges++;
                } else if (result == -3) {
                    cycle_walk_fail_no_path++;
                } else if (result == -4) {
                    cycle_walk_fail_no_cuts++;
                }
            } else {
                // Internal forest walk - only changes spanning trees
                result = internal_forest_walk(partition);

                diag_accept_prob.push_back(NA_REAL);
                diag_cycle_length.push_back(0);
                diag_n_valid_cuts.push_back(0);

                if (result == 0) {
                    forest_walk_success++;
                } else {
                    forest_walk_fail++;
                }
            }
        }

        // Copy current plan to output at thinning intervals
        if (i % thin == 0) {
            // Record whether any cycle_walk was accepted during these instep iterations
            mh_decisions(idx - 1) = any_cycle_walk_accepted ? 1 : 0;
            plans.col(idx + 1) = partition.get_plan();
            idx++;
        }

        // Update progress bar
        if (verbosity >= 1 && CLI_SHOULD_TICK) {
            cli_progress_set(bar, i - 1);
        }

        // Check for user interrupt
        if (i % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }

        // Break if we've filled the output
        if (idx == n_out - 1) {
            if (verbosity >= 1) cli_progress_set(bar, N);
            break;
        }
    }
    } catch (Rcpp::internal::InterruptedException& e) {
        if (verbosity >= 1) cli_progress_done(bar);
        throw;  // Re-throw to let R handle it
    }

    if (verbosity >= 2) {
        Rcout << "\n[Forest Walk] Success: " << forest_walk_success
              << ", Fail: " << forest_walk_fail << "\n";
        Rcout << "[Cycle Walk] Accept: " << cycle_walk_accept
              << ", Reject: " << cycle_walk_reject << "\n";
        int total_fail = cycle_walk_fail_no_adj + cycle_walk_fail_few_edges +
                         cycle_walk_fail_no_path + cycle_walk_fail_no_cuts;
        if (total_fail > 0) {
            Rcout << "[Cycle Walk Failures] Total: " << total_fail;
            if (cycle_walk_fail_few_edges > 0)
                Rcout << ", <2 edges: " << cycle_walk_fail_few_edges;
            if (cycle_walk_fail_no_cuts > 0)
                Rcout << ", no valid cuts: " << cycle_walk_fail_no_cuts;
            if (cycle_walk_fail_no_path > 0)
                Rcout << ", no path: " << cycle_walk_fail_no_path;
            if (cycle_walk_fail_no_adj > 0)
                Rcout << ", no adj: " << cycle_walk_fail_no_adj;
            Rcout << "\n";
        }
    }

    if (verbosity >= 1) {
        cli_progress_done(bar);
    }

    // Return results as a List
    // Create diagnostic list
    List diagnostics = List::create(
        Named("accept_prob") = NumericVector(diag_accept_prob.begin(), diag_accept_prob.end()),
        Named("cycle_length") = IntegerVector(diag_cycle_length.begin(), diag_cycle_length.end()),
        Named("n_valid_cuts") = IntegerVector(diag_n_valid_cuts.begin(), diag_n_valid_cuts.end()),
        Named("failure_modes") = List::create(
            Named("no_adj") = cycle_walk_fail_no_adj,
            Named("few_edges") = cycle_walk_fail_few_edges,
            Named("no_path") = cycle_walk_fail_no_path,
            Named("no_cuts") = cycle_walk_fail_no_cuts
        )
    );

    // Always return a List for easy iteration and extension
    List out = List::create(
        Named("plans") = plans,
        Named("mhdecisions") = mh_decisions,
        Named("diagnostics") = diagnostics,
        Named("algorithm") = "cyclewalk"
    );

    return out;
}
