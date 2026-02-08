/*
 * MEW (Marked Edge Walk) MCMC Redistricting Algorithm
 *
 * Based on: McWhorter, A., & DeFord, D. (2024)
 * "The Marked Edge Walk: A Novel MCMC Algorithm for Sampling of Graph Partitions"
 *
 * Implementation for the redist R package
 */

#include "mew.h"
#include "mew_helpers.h"
#include <algorithm>
#include <cmath>

/*
 * Main entry point for MEW algorithm
 */
// [[Rcpp::export]]
Rcpp::List mew_plans(int nsims, List adj, const arma::uvec &init,
                     const arma::uvec &counties, const arma::uvec &pop,
                     int n_distr, double target, double lower, double upper,
                     double rho, List constraints, List control,
                     int thin, int verbosity) {

    // Re-seed RNG for reproducibility across parallel chains
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // Convert adjacency list to Graph
    Graph g = list_to_graph(adj);
    int V = g.size();
    int n_cty = max(counties);

    // Initialize spanning tree and marked edges from partition
    // This is the correct approach: build tree structure from the partition
    // rather than trying to find marked edges in a random tree
    if (verbosity > 1) {
        Rcpp::Rcout << "Building spanning tree from initial partition..." << std::endl;
    }

    auto init_result = partition_to_tree_marked_edges(g, init, n_distr);
    Tree tree = init_result.first;
    MarkedEdgeSet marked_edges = init_result.second;

    // Extract initial partition from tree + marked edges
    // Note: This may not exactly match init plan due to randomness in Wilson's algorithm
    // within each district, but it will be a valid partition with the correct number of districts
    uvec partition_current = tree_to_partition(tree, marked_edges, V, n_distr);

    // Verify we have the correct number of districts
    std::set<int> unique_districts;
    for (int i = 0; i < V; i++) {
        unique_districts.insert(partition_current(i));
    }

    if ((int)unique_districts.size() != n_distr) {
        std::ostringstream msg;
        msg << "Initialization produced " << unique_districts.size()
            << " districts instead of " << n_distr;
        Rcpp::stop(msg.str());
    }

    if (verbosity > 1) {
        Rcpp::Rcout << "Initialization successful: " << marked_edges.size()
                    << " marked edges, " << n_distr << " districts" << std::endl;
    }

    // Storage for results
    arma::umat plans(V, nsims);
    std::vector<double> accept_rate_vec;
    std::vector<double> cycle_intersect_rate_vec;
    std::vector<double> avg_proposal_tries_vec;

    // Gibbs constraint infrastructure
    Rcpp::NumericVector new_psi;
    if (constraints.size() > 0) {
        Rcpp::CharacterVector constr_names = constraints.names();
        new_psi = Rcpp::NumericVector(constr_names.size());
        new_psi.names() = constr_names;
    }

    // 2-column umat for kirchhoff log_st and calc_gibbs_tgt calls
    // Column 0 = old plan, Column 1 = new plan
    arma::umat plans_pair(V, 2);

    // Tracking variables
    int n_accept = 0;
    int n_cycle_intersect = 0;
    int total_tries = 0;
    int n_failed_proposals = 0;  // Count proposals that hit MAX_TRIES

    // Progress bar
    SEXP pb = R_NilValue;
    bool pb_active = false;
    if (verbosity > 0) {
        List pb_conf = cli_config(false,
            "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | {cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta}");
        pb = cli_progress_bar(nsims, pb_conf);
        pb_active = true;
    }

    // Main MCMC loop
    try {
        for (int iter = 0; iter < nsims; iter++) {
            // Save current state for potential rejection
            Tree tree_old = tree;
            MarkedEdgeSet marked_old = marked_edges;

            // Make proposal
            MEWProposal proposal = mew_proposal(g, tree, marked_edges, pop, n_distr,
                                               target, lower, upper);
            total_tries += proposal.n_rejects + 1;

            // Only accept/reject if proposal is valid
            // Invalid proposals (timed out after MAX_TRIES) are automatically rejected
            if (proposal.valid) {
                // Compute transition probability ratio
                double trans_prob = transition_probability(
                    proposal.cycle.cycle_edges,
                    proposal.cycle.edge_plus,
                    proposal.marked.old_edge,
                    proposal.marked.new_edge,
                    marked_old,
                    proposal.marked.marked_new,
                    tree_old,
                    proposal.cycle.tree_new
                );
                double log_q_ratio = std::log(trans_prob);

                // 2. Get proposed partition (cached from population check)
                uvec partition_new = proposal.partition;

                // Pack into 2-column matrix: col 0 = old, col 1 = new
                plans_pair.col(0) = partition_current;
                plans_pair.col(1) = partition_new;

                // 3. Detect which districts changed
                std::vector<int> changed_districts;
                {
                    std::vector<bool> district_changed(n_distr + 1, false);
                    for (int v = 0; v < V; v++) {
                        if (partition_current(v) != partition_new(v)) {
                            district_changed[partition_current(v)] = true;
                            district_changed[partition_new(v)] = true;
                        }
                    }
                    for (int d = 1; d <= n_distr; d++) {
                        if (district_changed[d]) {
                            changed_districts.push_back(d);
                        }
                    }
                }

                // 4. Log spanning tree correction (Kirchhoff matrix-tree theorem)
                double log_st = 0.0;
                if (rho != 1) {
                    for (int d : changed_districts) {
                        for (int j = 1; j <= n_cty; j++) {
                            log_st += log_st_distr(g, plans_pair, counties, 0, d, j);
                            log_st -= log_st_distr(g, plans_pair, counties, 1, d, j);
                        }
                        log_st += log_st_contr(g, plans_pair, counties, n_cty, 0, d);
                        log_st -= log_st_contr(g, plans_pair, counties, n_cty, 1, d);
                    }
                }

                // 5. Gibbs target constraints
                double log_energy_ratio = 0.0;
                if (constraints.size() > 0) {
                    std::fill(new_psi.begin(), new_psi.end(), 0.0);
                    double old_tgt = calc_gibbs_tgt(
                        plans_pair.col(0), n_distr, V,
                        changed_districts, new_psi, pop, target, g, constraints);

                    std::fill(new_psi.begin(), new_psi.end(), 0.0);
                    double new_tgt = calc_gibbs_tgt(
                        plans_pair.col(1), n_distr, V,
                        changed_districts, new_psi, pop, target, g, constraints);

                    // Positive when new plan is better (lower cost)
                    log_energy_ratio = old_tgt - new_tgt;
                }

                // 6. Combined log acceptance probability
                double log_alpha = log_q_ratio + (1.0 - rho) * log_st + log_energy_ratio;

                // Accept/reject
                bool accept = (std::log(r_unif()) < log_alpha);

                if (accept) {
                    tree = proposal.cycle.tree_new;
                    marked_edges = proposal.marked.marked_new;
                    partition_current = partition_new;
                    n_accept++;

                    // Track cycle-marked edge intersection
                    if (log_q_ratio != 0.0) {
                        n_cycle_intersect++;
                    }
                }
                // If reject, tree and marked_edges stay as tree_old and marked_old
            } else {
                // Proposal failed (hit MAX_TRIES)
                n_failed_proposals++;
                if (n_failed_proposals >= 50 && verbosity > 0 && n_failed_proposals % 50 == 0) {
                    Rcpp::Rcout << "Warning: " << n_failed_proposals
                               << " proposals have failed to meet population constraints. "
                               << "Chain may be stuck." << std::endl;
                }
            }
            // If !proposal.valid, automatically reject (keep current state)

            // Store plan
            plans.col(iter) = partition_current;

            // Update diagnostics periodically
            if ((iter + 1) % 100 == 0) {
                double curr_accept = (double)n_accept / (iter + 1);
                accept_rate_vec.push_back(curr_accept);

                double curr_cycle_intersect = (double)n_cycle_intersect / (iter + 1);
                cycle_intersect_rate_vec.push_back(curr_cycle_intersect);

                double curr_avg_tries = (double)total_tries / (iter + 1);
                avg_proposal_tries_vec.push_back(curr_avg_tries);
            }

            // Update progress
            if (pb_active && (iter + 1) % 10 == 0) {
                cli_progress_set(pb, iter + 1);
            }

            // Check for user interrupt
            if ((iter + 1) % 100 == 0) {
                Rcpp::checkUserInterrupt();
            }
        }
    } catch (Rcpp::internal::InterruptedException& e) {
        if (pb_active) {
            cli_progress_done(pb);
        }
        throw;  // Re-throw to let R handle it
    }

    // Clean up progress bar before returning
    if (pb_active) {
        cli_progress_done(pb);
        pb_active = false;
    }

    // Final diagnostics
    if (nsims > 0) {
        double final_accept = (double)n_accept / nsims;
        accept_rate_vec.push_back(final_accept);

        double final_cycle_intersect = (double)n_cycle_intersect / nsims;
        cycle_intersect_rate_vec.push_back(final_cycle_intersect);

        double final_avg_tries = (double)total_tries / nsims;
        avg_proposal_tries_vec.push_back(final_avg_tries);
    }

    // Apply thinning if requested
    arma::umat plans_thinned;
    if (thin > 1) {
        int n_keep = nsims / thin;
        plans_thinned = arma::umat(V, n_keep);
        for (int i = 0; i < n_keep; i++) {
            plans_thinned.col(i) = plans.col(i * thin);
        }
    } else {
        plans_thinned = plans;
    }

    // Return results
    return List::create(
        Named("plans") = plans_thinned,
        Named("accept_rate") = accept_rate_vec,
        Named("cycle_intersect_rate") = cycle_intersect_rate_vec,
        Named("avg_proposal_tries") = avg_proposal_tries_vec
    );
}
