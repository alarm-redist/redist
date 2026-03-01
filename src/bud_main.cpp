/********************************************************
 * BUD MCMC Redistricting Sampler
 * Main entry point and MCMC loop
 *
 * Modeled on cw_main.cpp (CycleWalk MCMC loop)
 ********************************************************/

#include "bud_main.h"
#include "bud_partition.h"
#include "bud_proposal.h"

#ifdef _WIN32
#include <R.h>
extern "C" void R_ProcessEvents(void);
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List bud_plans(
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

    // Progress bar setup
    RObject bar = R_NilValue;
    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << "Balanced Up-Down Walk\n";
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
    }

    // Initialize BUD Partition
    if (verbosity >= 1) {
        Rcout << "\nInitializing partition with spanning trees...\n";
    }
    BUDPartition partition(V, n_distr);
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
        partition.print_bud_state(verbosity);
        Rcout << "\n";
        bar = cli_progress_bar(N, cli_config(false));
    }

    // Counters
    int bud_accept = 0;
    int bud_reject = 0;
    int bud_trivial = 0;

    int idx = 1;
    try {
    for (int i = 1; i <= N; i++) {
        bool any_accepted = false;

        if (verbosity >= 3) {
            Rcout << "[D] iter " << i << " start\n";
        }

        for (int step = 0; step < instep; step++) {
            double accept_prob = 0.0;
            int result = bud_step(partition, lower, upper, target,
                                   compactness, counties, constraints,
                                   accept_prob);

            if (result == 1) {
                bud_accept++;
                any_accepted = true;
            } else if (result == 0) {
                if (accept_prob >= 0.99) {
                    bud_trivial++;
                    any_accepted = true; // trivial swaps count
                } else {
                    bud_reject++;
                }
            }
        }

        // Copy current plan at thinning intervals
        if (i % thin == 0) {
            mh_decisions(idx - 1) = any_accepted ? 1 : 0;
            plans.col(idx + 1) = partition.get_plan();
            idx++;
        }

        // Update progress bar and check for user interrupt
        if (verbosity >= 1) {
            if (CLI_SHOULD_TICK) {
                cli_progress_set(bar, i - 1);
            }
        }
        Rcpp::checkUserInterrupt();
#ifdef _WIN32
        R_ProcessEvents();
#endif

        // Break if we've filled the output
        if (idx == n_out - 1) {
            if (verbosity >= 1) cli_progress_set(bar, N);
            break;
        }
    }
    } catch (Rcpp::internal::InterruptedException& e) {
        if (verbosity >= 1) cli_progress_done(bar);
        throw;
    }

    if (verbosity >= 2) {
        Rcout << "\n[BUD] Accept: " << bud_accept
              << ", Reject: " << bud_reject
              << ", Trivial: " << bud_trivial << "\n";
    }

    if (verbosity >= 1) {
        cli_progress_done(bar);
    }

    // Create diagnostics
    List diagnostics = List::create(
        Named("accept") = bud_accept,
        Named("reject") = bud_reject,
        Named("trivial") = bud_trivial
    );

    List out = List::create(
        Named("plans") = plans,
        Named("mhdecisions") = mh_decisions,
        Named("diagnostics") = diagnostics,
        Named("algorithm") = "bud"
    );

    return out;
}
