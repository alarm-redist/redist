/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Main entry point and MCMC loop
 *
 * ITERATION 5: Constraint Integration
 * - Full cycle walk proposal with MH acceptance
 * - Constraints integrated via calc_gibbs_tgt
 * - Mixed with internal forest walk for better mixing
 ********************************************************/

#include "cw_main.h"
#include "cw_lct.h"
#include "cw_partition.h"
#include "cw_forest_walk.h"
#include "cw_proposal.h"

using namespace Rcpp;
using namespace arma;

/*
 * Main entry point for CycleWalk MCMC sampler.
 *
 * STUB: Currently just returns the initial plan N/thin times.
 * This allows the R interface to be tested before implementing
 * the actual algorithm.
 */
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
    Rcpp::List edge_weights,  // NEW: Edge weights
    int thin,
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

    // === ITERATION 2: Initialize LCT Partition ===
    LCTPartition partition(V, n_distr);
    int init_result = partition.init_from_plan(g, init, pop, counties, lower, upper);

    if (init_result != 0) {
        Rcpp::stop("Failed to initialize partition - Wilson's algorithm failed.");
    }

    // Set edge weights if provided
    if (edge_weights.size() > 0) {
        partition.set_edge_weights(edge_weights);
    }

    if (verbosity >= 1) {
        bar = cli_progress_bar(N, cli_config(false));
    }

    // Main MCMC loop
    // Mix cycle walk (10%) with internal forest walk (90%)
    int forest_walk_success = 0;
    int forest_walk_fail = 0;
    int cycle_walk_accept = 0;
    int cycle_walk_reject = 0;
    int cycle_walk_fail_no_adj = 0;      // -1: no adjacent districts
    int cycle_walk_fail_few_edges = 0;   // -2: <2 boundary edges
    int cycle_walk_fail_no_path = 0;     // -3: couldn't get paths
    int cycle_walk_fail_no_cuts = 0;     // -4: no valid cut pairs

    double cycle_walk_prob = 0.1;  // 10% cycle walk, 90% forest walk

    int idx = 1;
    for (int i = 1; i <= N; i++) {
        int result;

        if (r_unif() < cycle_walk_prob) {
            // Cycle walk proposal - can change districts
            CycleWalkDiagnostics diag;
            result = cycle_walk(partition, lower, upper, target, constraints, diag);
            
            // Collect diagnostics
            diag_accept_prob.push_back(diag.accept_prob);
            diag_cycle_length.push_back(diag.cycle_length);
            diag_n_valid_cuts.push_back(diag.n_valid_cuts);
            
            // Update failure mode counters
            if (result == 1) {
                cycle_walk_accept++;
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
            
            // Forest walk doesn't generate cycle walk diagnostics - use NA/0
            diag_accept_prob.push_back(NA_REAL);
            diag_cycle_length.push_back(0);
            diag_n_valid_cuts.push_back(0);
            
            if (result == 0) {
                forest_walk_success++;
            } else {
                forest_walk_fail++;
            }
        }

        // Record MH decision (1 = accept, 0 = reject)
        mh_decisions(idx - 1) = (result >= 0) ? 1 : 0;

        // Copy current plan to output at thinning intervals
        // Note: District assignments don't change with internal forest walk
        // but we still output the plan for consistency
        if (i % thin == 0) {
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
        Named("algorithm") = "cyclewalk",
        Named("version") = "0.5-constraints"
    );

    return out;
}

/*
 * Test function for Link-Cut Tree operations.
 * Used to validate the LCT implementation.
 *
 * Creates a simple tree:
 *     0
 *    / \
 *   1   2
 *  / \
 * 3   4
 *
 * Then tests link, cut, evert, find_root, find_path, same_tree.
 */
// [[Rcpp::export]]
Rcpp::List test_lct() {
    LinkCutTree lct(5);
    std::vector<std::string> tests;
    std::vector<bool> passed;

    // Test 1: Initially, each node is its own root
    tests.push_back("Initial: each node is own root");
    bool test1 = (lct.find_root(0) == 0 && lct.find_root(1) == 1 &&
                  lct.find_root(2) == 2 && lct.find_root(3) == 3 &&
                  lct.find_root(4) == 4);
    passed.push_back(test1);

    // Test 2: Initially, nodes are in different trees
    tests.push_back("Initial: nodes in different trees");
    bool test2 = (!lct.same_tree(0, 1) && !lct.same_tree(1, 2));
    passed.push_back(test2);

    // Build tree: link 1 to 0, link 2 to 0, link 3 to 1, link 4 to 1
    lct.evert(1);
    lct.link(1, 0);
    lct.evert(2);
    lct.link(2, 0);
    lct.evert(3);
    lct.link(3, 1);
    lct.evert(4);
    lct.link(4, 1);

    // Test 3: After linking, root should be 0
    tests.push_back("After linking: root is 0");
    bool test3 = (lct.find_root(0) == 0 && lct.find_root(1) == 0 &&
                  lct.find_root(2) == 0 && lct.find_root(3) == 0 &&
                  lct.find_root(4) == 0);
    passed.push_back(test3);

    // Test 4: All nodes in same tree
    tests.push_back("After linking: all in same tree");
    bool test4 = (lct.same_tree(0, 1) && lct.same_tree(1, 2) &&
                  lct.same_tree(2, 3) && lct.same_tree(3, 4));
    passed.push_back(test4);

    // Test 5: Path from root to node 3 should be [0, 1, 3]
    tests.push_back("Path from root to 3 is [0,1,3]");
    std::vector<int> path3 = lct.find_path(3);
    bool test5 = (path3.size() == 3 && path3[0] == 0 && path3[1] == 1 && path3[2] == 3);
    passed.push_back(test5);

    // Test 6: Path from root to node 2 should be [0, 2]
    tests.push_back("Path from root to 2 is [0,2]");
    std::vector<int> path2 = lct.find_path(2);
    bool test6 = (path2.size() == 2 && path2[0] == 0 && path2[1] == 2);
    passed.push_back(test6);

    // Test 7: Evert node 3 to make it root
    tests.push_back("Evert 3: new root is 3");
    lct.evert(3);
    bool test7 = (lct.find_root(0) == 3 && lct.find_root(1) == 3 &&
                  lct.find_root(3) == 3);
    passed.push_back(test7);

    // Test 8: Path from new root (3) to node 0 should be [3, 1, 0]
    tests.push_back("After evert: path to 0 is [3,1,0]");
    std::vector<int> path0 = lct.find_path(0);
    bool test8 = (path0.size() == 3 && path0[0] == 3 && path0[1] == 1 && path0[2] == 0);
    passed.push_back(test8);

    // Test 9: Cut node 1 from tree (separates 0,2 from 1,3,4)
    tests.push_back("Cut 0: separates trees");
    lct.cut(0);  // Cut 0 from 1
    bool test9 = (!lct.same_tree(0, 1) && !lct.same_tree(0, 3) &&
                  lct.same_tree(1, 3) && lct.same_tree(1, 4));
    passed.push_back(test9);

    // Test 10: Node 0 is now its own root (was cut)
    tests.push_back("After cut: 0 is own root");
    bool test10 = (lct.find_root(0) == 0);
    passed.push_back(test10);

    // Test 11: Node 2 is now connected to 0 only
    // Wait, 2 was linked to 0. After cutting 0, is 2 with 0?
    tests.push_back("After cut: 2 with 0 or separate");
    bool test11 = lct.same_tree(0, 2);
    passed.push_back(test11);

    // Count passes
    int n_passed = 0;
    for (bool p : passed) if (p) n_passed++;

    // Convert path3 and path0 to R vectors for debugging
    IntegerVector r_path3(path3.size());
    for (size_t i = 0; i < path3.size(); i++) r_path3[i] = path3[i];

    IntegerVector r_path0(path0.size());
    for (size_t i = 0; i < path0.size(); i++) r_path0[i] = path0[i];

    IntegerVector r_path2(path2.size());
    for (size_t i = 0; i < path2.size(); i++) r_path2[i] = path2[i];

    return List::create(
        Named("tests") = tests,
        Named("passed") = passed,
        Named("n_passed") = n_passed,
        Named("n_tests") = (int)tests.size(),
        Named("path_to_3") = r_path3,
        Named("path_to_2") = r_path2,
        Named("path_to_0_after_evert") = r_path0
    );
}
