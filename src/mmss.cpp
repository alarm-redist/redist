#include "mmss.h"

/*
 * Select `l` connected districts from the district graph via random BFS.
 * Returns a vector of 1-indexed district labels.
 * `log_prob` is set to log(probability of this selection).
 */
std::vector<int> select_l_districts(int n_distr, const Graph &dist_g,
                                    int l, double &log_prob) {
    int start = r_int(n_distr);
    std::vector<int> selected;
    selected.push_back(start + 1);

    log_prob = -std::log((double) n_distr);

    std::vector<bool> in_set(n_distr, false);
    in_set[start] = true;

    for (int step = 1; step < l; step++) {
        std::vector<std::pair<int,int>> frontier;
        for (int d : selected) {
            int d0 = d - 1;
            for (int nbor : dist_g[d0]) {
                if (!in_set[nbor]) {
                    frontier.push_back({d0, nbor});
                }
            }
        }

        if (frontier.empty()) {
            log_prob = -std::numeric_limits<double>::infinity();
            return selected;
        }

        int pick = r_int(frontier.size());
        int chosen = frontier[pick].second;
        selected.push_back(chosen + 1);

        int edges_to_chosen = 0;
        for (auto &e : frontier) {
            if (e.second == chosen) edges_to_chosen++;
        }
        log_prob += std::log((double) edges_to_chosen / frontier.size());

        in_set[chosen] = true;
    }

    return selected;
}


/*
 * Compute log probability of selecting districts in a SPECIFIC order.
 * For the MH ratio with l>2, different orderings produce different
 * splits, so we must use the permutation-specific probability (not
 * the set probability summed over all orderings).
 */
static double log_prob_perm(const std::vector<int> &districts, int n_distr,
                            const Graph &dist_g) {
    int l = districts.size();
    double log_prob = -std::log((double) n_distr);

    std::vector<bool> in_set(n_distr, false);
    in_set[districts[0] - 1] = true;

    for (int step = 1; step < l; step++) {
        int target_d = districts[step] - 1;

        int frontier_size = 0;
        int edges_to_target = 0;
        for (int s = 0; s < step; s++) {
            int d0 = districts[s] - 1;
            for (int nbor : dist_g[d0]) {
                if (!in_set[nbor]) {
                    frontier_size++;
                    if (nbor == target_d) edges_to_target++;
                }
            }
        }

        if (edges_to_target == 0 || frontier_size == 0) {
            return -std::numeric_limits<double>::infinity();
        }
        log_prob += std::log((double) edges_to_target / frontier_size);
        in_set[target_d] = true;
    }
    return log_prob;
}


/*
 * Cut a spanning tree to peel off one district from a two-label region.
 * Uses a RANDOM EDGE mechanism: pick a uniformly random non-root vertex
 * from the spanning tree and cut the edge to its parent.
 *
 * This gives q(A,B) ∝ T(A)*T(B)*B(A,B) by the matrix-tree theorem,
 * regardless of whether the population bounds are symmetric or asymmetric.
 * The T factors cancel in the MH ratio, leaving only the boundary ratio.
 *
 * For l=2 (symmetric bounds), this is equivalent to the original merge-split
 * mechanism.  For l>2 with asymmetric intermediate steps, this is correct
 * whereas root-district-only is biased.
 */
static bool cut_one_mms(Tree &ust, int k, int root,
                        subview_col<uword> &districts,
                        int peel_label, int remain_label,
                        const uvec &pop, double total_pop,
                        double peel_lower, double peel_upper,
                        double remain_lower, double remain_upper,
                        double peel_target) {
    int V = ust.size();

    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
    parent[root] = -1;
    tree_pop(ust, root, pop, pop_below, parent);

    // Collect non-root vertices in the merged region
    std::vector<int> region_verts;
    for (int i = 0; i < V; i++) {
        if (i == root) continue;
        if ((int) districts(i) == peel_label || (int) districts(i) == remain_label)
            region_verts.push_back(i);
    }

    if (region_verts.empty()) return false;

    // Pick a random vertex (= random edge in the tree)
    int pick = r_int(region_verts.size());
    int cut_at = region_verts[pick];
    double below = pop_below[cut_at];
    double above = total_pop - below;

    // Check which assignment of labels is valid
    bool subtree_is_peel = (peel_lower <= below && below <= peel_upper &&
                            remain_lower <= above && above <= remain_upper);
    bool subtree_is_remain = (remain_lower <= below && below <= remain_upper &&
                              peel_lower <= above && above <= peel_upper);

    if (!subtree_is_peel && !subtree_is_remain) return false;

    // Remove edge from tree
    std::vector<int> *siblings = &ust[parent[cut_at]];
    for (int j = 0; j < (int) siblings->size(); j++) {
        if ((*siblings)[j] == cut_at) {
            siblings->erase(siblings->begin() + j);
            break;
        }
    }
    parent[cut_at] = -1;

    // Assign labels based on which assignment is valid
    if (subtree_is_peel) {
        assign_district(ust, districts, cut_at, peel_label);
        assign_district(ust, districts, root, remain_label);
    } else {
        assign_district(ust, districts, cut_at, remain_label);
        assign_district(ust, districts, root, peel_label);
    }

    return true;
}


/*
 * Main MMSS MCMC loop.
 */
// [[Rcpp::export]]
Rcpp::List mmss_plans(int N, List l, const arma::uvec init, const arma::uvec &counties,
                     const arma::uvec &pop, int n_distr, double target, double lower,
                     double upper, double rho, List constraints, List control,
                     int k, int thin, int l_merge, int verbosity) {
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    double thresh = (double) control["adapt_k_thresh"];
    bool do_mh = (bool) control["do_mh"];

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int n_cty = max(counties);

    int n_out = N / thin + 2;
    umat districts(V, n_out, fill::zeros);
    districts.col(0) = init;
    districts.col(1) = init;

    Rcpp::IntegerVector mh_decisions(N / thin + 1);
    double mha;

    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << "MARKOV CHAIN MONTE CARLO (Multi-Merge-Split, l=" << l_merge << ")\n";
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
    }

    if (k <= 0) {
        adapt_ms_parameters(g, n_distr, k, thresh, tol, init, counties, cg, pop, target);
    }
    if (verbosity >= 3)
        Rcout << "Using k = " << k << "\n";

    Graph dist_g = district_graph(g, init, n_distr);
    Graph new_dist_g;
    int n_accept = 0;

    CharacterVector psi_names = CharacterVector::create(
        "pop_dev", "splits", "multisplits", "total_splits",
        "segregation", "grp_pow", "grp_hinge", "grp_inv_hinge",
        "compet", "status_quo", "incumbency",
        "polsby", "fry_hold", "log_st", "edges_removed",
        "qps", "custom"
    );
    NumericVector new_psi(psi_names.size());
    new_psi.names() = psi_names;

    RObject bar = cli_progress_bar(N, cli_config(false));
    int idx = 1;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V);

    for (int i = 1; i <= N; i++) {
        double prop_lp = 0.0;
        mh_decisions(idx - 1) = 0;
        districts.col(idx + 1) = districts.col(idx);

        // 1. Select l connected districts
        double log_prob_fwd = 0.0;
        std::vector<int> sel_districts = select_l_districts(
            n_distr, dist_g, l_merge, log_prob_fwd);

        if (!std::isfinite(log_prob_fwd)) {
            districts.col(idx + 1) = districts.col(idx);
            if (i % thin == 0) idx++;
            continue;
        }

        // Save the current plan for label restoration
        uvec saved_plan = districts.col(idx + 1);

        // 2. Sequential split with retries of the ENTIRE sequence.
        // Retrying the complete split as one unit preserves q ∝ T*T*B
        // because each complete attempt is independent and the overall
        // success probability p_any depends only on the merged region
        // (same for forward and reverse), so the retry constant cancels
        // in the MH ratio.
        bool split_failed = true;
        double fwd_boundary_lp = 0.0;
        const int MAX_SPLIT_ATTEMPTS = 500;

        for (int attempt = 0; attempt < MAX_SPLIT_ATTEMPTS; attempt++) {
            // Reset to pre-merge state for each complete attempt
            districts.col(idx + 1) = saved_plan;
            fwd_boundary_lp = 0.0;
            bool attempt_ok = true;

            for (int s = 0; s < l_merge - 1; s++) {
                int peel = sel_districts[s];
                int remain = sel_districts[l_merge - 1];

                // Temporarily relabel "middle" remaining districts -> remain label
                for (int v = 0; v < V; v++) {
                    for (int t = s + 1; t < l_merge - 1; t++) {
                        if ((int) districts(v, idx + 1) == sel_districts[t]) {
                            districts(v, idx + 1) = remain;
                            break;
                        }
                    }
                }

                // set ignore for vertices not in the two-label region
                double region_pop = 0.0;
                for (int v = 0; v < V; v++) {
                    ignore[v] = ((int) districts(v, idx + 1) != peel &&
                                 (int) districts(v, idx + 1) != remain);
                    if (!ignore[v]) region_pop += pop(v);
                }

                // population bounds for this split step
                int remaining_splits = l_merge - 1 - s;
                double peel_lower = std::max(lower, region_pop - remaining_splits * upper);
                double peel_upper = std::min(upper, region_pop - remaining_splits * lower);
                double remain_lower, remain_upper;
                if (remaining_splits > 1) {
                    remain_lower = remaining_splits * lower;
                    remain_upper = remaining_splits * upper;
                } else {
                    remain_lower = lower;
                    remain_upper = upper;
                }

                if (peel_lower > peel_upper) {
                    attempt_ok = false;
                    break;
                }

                // Sample UST and try random edge cut (single try per step)
                double ust_lower = std::min(peel_lower, remain_lower);
                double ust_upper = std::max(peel_upper, remain_upper);
                int root;
                clear_tree(ust);
                int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                            pop, ust_lower, ust_upper, counties, cg);
                if (result != 0) { attempt_ok = false; break; }

                auto col_ref = districts.col(idx + 1);
                if (!cut_one_mms(ust, k, root, col_ref,
                                 peel, remain, pop, region_pop,
                                 peel_lower, peel_upper,
                                 remain_lower, remain_upper, target)) {
                    attempt_ok = false;
                    break;
                }

                // forward boundary AFTER cut
                fwd_boundary_lp += log_boundary(g, districts.col(idx + 1), peel, remain);

                // Restore remaining vertices' original labels for next step
                for (int v = 0; v < V; v++) {
                    if ((int) districts(v, idx + 1) == remain) {
                        int orig = saved_plan(v);
                        for (int t = s + 1; t < l_merge; t++) {
                            if (orig == sel_districts[t]) {
                                districts(v, idx + 1) = orig;
                                break;
                            }
                        }
                    }
                }
            } // end steps loop

            if (attempt_ok) {
                split_failed = false;
                break;
            }
        } // end retry loop

        if (split_failed) {
            districts.col(idx + 1) = districts.col(idx);
            if (i % thin == 0) idx++;
            continue;
        }

        // Ensure all vertices in the selected set are properly assigned
        int last_label = sel_districts[l_merge - 1];
        for (int v = 0; v < V; v++) {
            bool in_selected = false;
            for (int t = 0; t < l_merge; t++) {
                if ((int) saved_plan(v) == sel_districts[t]) {
                    in_selected = true;
                    break;
                }
            }
            if (!in_selected) continue;
            bool peeled = false;
            for (int t = 0; t < l_merge - 1; t++) {
                if ((int) districts(v, idx + 1) == sel_districts[t]) {
                    peeled = true;
                    break;
                }
            }
            if (!peeled) {
                districts(v, idx + 1) = last_label;
            }
        }

        // 3. Reverse proposal boundary terms
        double rev_boundary_lp = 0.0;
        {
            uvec old_plan = districts.col(idx);
            umat work_mat(V, 1);
            work_mat.col(0) = old_plan;
            // merge old districts to 0
            for (int v = 0; v < V; v++) {
                for (int d : sel_districts) {
                    if ((int) work_mat(v, 0) == d) {
                        work_mat(v, 0) = 0;
                        break;
                    }
                }
            }
            // reveal old plan's districts sequentially
            for (int s = 0; s < l_merge - 1; s++) {
                int dist_label = sel_districts[s];
                for (int v = 0; v < V; v++) {
                    if ((int) old_plan(v) == dist_label && work_mat(v, 0) == 0) {
                        work_mat(v, 0) = dist_label;
                    }
                }
                rev_boundary_lp += log_boundary(g, work_mat.col(0), 0, dist_label);
            }
        }

        // proposal ratio = log q(y->x) - log q(x->y)
        prop_lp = rev_boundary_lp - fwd_boundary_lp;

        // 4. Compactness (tau)
        if (rho != 1) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                for (int d : sel_districts) {
                    log_st += log_st_distr(g, districts, counties, idx, d, j);
                    log_st -= log_st_distr(g, districts, counties, idx + 1, d, j);
                }
            }
            for (int d : sel_districts) {
                log_st += log_st_contr(g, districts, counties, n_cty, idx, d);
                log_st -= log_st_contr(g, districts, counties, n_cty, idx + 1, d);
            }
            prop_lp += (1 - rho) * log_st;
        }

        // 5. Gibbs constraint target
        std::vector<int> distr_vec(sel_districts.begin(), sel_districts.end());
        double gibbs_new = calc_gibbs_tgt(districts.col(idx + 1), n_distr, V, distr_vec,
                                  new_psi, pop, target, g, constraints);
        double gibbs_old = calc_gibbs_tgt(districts.col(idx), n_distr, V, distr_vec,
                                  new_psi, pop, target, g, constraints);
        prop_lp += gibbs_old - gibbs_new;

        // 6. Selection probability correction
        // The last two selected districts can be swapped at the final step
        // (symmetric population bounds), producing the same plan via opposite
        // orientation. Must sum over both orderings for correctness.
        new_dist_g = district_graph(g, districts.col(idx + 1), n_distr);

        std::vector<int> sel_swapped = sel_districts;
        std::swap(sel_swapped[l_merge - 2], sel_swapped[l_merge - 1]);

        double lp_fwd_swap = log_prob_perm(sel_swapped, n_distr, dist_g);
        double mx_fwd = std::max(log_prob_fwd, lp_fwd_swap);
        double log_psum_fwd = mx_fwd + std::log(
            std::exp(log_prob_fwd - mx_fwd) + std::exp(lp_fwd_swap - mx_fwd));

        double lp_rev_orig = log_prob_perm(sel_districts, n_distr, new_dist_g);
        double lp_rev_swap = log_prob_perm(sel_swapped, n_distr, new_dist_g);
        double mx_rev = std::max(lp_rev_orig, lp_rev_swap);
        double log_psum_rev = mx_rev + std::log(
            std::exp(lp_rev_orig - mx_rev) + std::exp(lp_rev_swap - mx_rev));

        prop_lp += log_psum_rev - log_psum_fwd;

        // 7. Accept/reject
        if (!do_mh || prop_lp >= 0 || std::log(r_unif()) <= prop_lp) {
            n_accept++;
            districts.col(idx) = districts.col(idx + 1);
            dist_g = new_dist_g;
            mh_decisions(idx - 1) = 1;
        } else {
            districts.col(idx + 1) = districts.col(idx);
        }

        if (i % thin == 0) idx++;

        if (verbosity >= 1 && CLI_SHOULD_TICK) {
            cli_progress_set(bar, i - 1);
            mha = (double) n_accept / (i - 1);
            cli_progress_set_format(bar,
                "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} | MH Acceptance: %.2f", mha);
        }
        if (idx == n_out - 1) {
            cli_progress_set(bar, N);
            break;
        }
        Rcpp::checkUserInterrupt();
    }
    cli_progress_done(bar);

    if (verbosity >= 1) {
        Rcout << "Acceptance rate: " << std::setprecision(2)
              << (100.0 * n_accept) / (N - 1) << "%\n";
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["est_k"] = k;
    out["mhdecisions"] = mh_decisions;

    return out;
}
