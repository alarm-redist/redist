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
 *
 * Two picking strategies are supported via `from_valid_only`:
 *
 * FALSE (exact path): pick a uniformly random non-root vertex from ALL region
 *   vertices and return false if it is not a valid cut.  This gives
 *   q(A,B) ∝ T(A)*T(B)*B(A,B) by the matrix-tree theorem, so the T factors
 *   cancel in the MH ratio and only log_boundary() is needed.  The caller
 *   retries the entire sequence on failure (whole-sequence retry).
 *
 * TRUE (approximate path): enumerate valid-cut vertices first, then pick
 *   uniformly from them only — never returning false on a tree that has at
 *   least one valid cut.  This changes the proposal distribution by a factor
 *   of 1/p_s^edge relative to the exact proposal, which the caller corrects
 *   via log_p_edge_estimate() for steps s >= 1 (per-step retry).
 */
static bool cut_one_mms(Tree &ust, int root,
                        subview_col<uword> &districts,
                        int peel_label, int remain_label,
                        const uvec &pop, double total_pop,
                        double peel_lower, double peel_upper,
                        double remain_lower, double remain_upper,
                        double peel_target,
                        bool from_valid_only) {
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

    int cut_at;
    if (from_valid_only) {
        // Approximate path: enumerate valid cuts, pick uniformly from them.
        std::vector<int> valid_verts;
        for (int v : region_verts) {
            double below = pop_below[v];
            double above = total_pop - below;
            bool ok = (peel_lower <= below && below <= peel_upper &&
                       remain_lower <= above && above <= remain_upper) ||
                      (remain_lower <= below && below <= remain_upper &&
                       peel_lower <= above && above <= peel_upper);
            if (ok) valid_verts.push_back(v);
        }
        if (valid_verts.empty()) return false;
        cut_at = valid_verts[r_int(valid_verts.size())];
    } else {
        // Exact path: pick a random vertex from all region vertices.
        // Return false if it is not a valid cut; the caller retries the whole sequence.
        cut_at = region_verts[r_int(region_verts.size())];
        double below = pop_below[cut_at];
        double above = total_pop - below;
        bool ok = (peel_lower <= below && below <= peel_upper &&
                   remain_lower <= above && above <= remain_upper) ||
                  (remain_lower <= below && below <= remain_upper &&
                   peel_lower <= above && above <= peel_upper);
        if (!ok) return false;
    }

    double below = pop_below[cut_at];
    double above = total_pop - below;
    bool subtree_is_peel = (peel_lower <= below && below <= peel_upper &&
                            remain_lower <= above && above <= remain_upper);

    // Remove edge from tree
    std::vector<int> *siblings = &ust[parent[cut_at]];
    for (int j = 0; j < (int) siblings->size(); j++) {
        if ((*siblings)[j] == cut_at) {
            siblings->erase(siblings->begin() + j);
            break;
        }
    }
    parent[cut_at] = -1;

    // Assign labels
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
 * Estimate p_s^edge = E[fraction of valid-cut edges per UST] on a subgraph.
 * Draws K spanning trees, counts valid cuts in each, returns log(mean fraction).
 * Used to correct the MH ratio for per-step retry (steps s >= 1).
 */
static double log_p_edge_estimate(const Graph &g, const uvec &pop,
                                  const std::vector<bool> &ignore,
                                  double region_pop,
                                  double peel_lower, double peel_upper,
                                  double remain_lower, double remain_upper,
                                  const uvec &counties, Multigraph &cg,
                                  int K_est) {
    int V = g.size();
    int region_size = 0;
    for (int v = 0; v < V; v++) if (!ignore[v]) region_size++;
    if (region_size <= 1) return 0.0;

    int total_valid = 0;
    int trees_drawn = 0;
    Tree ust_tmp = init_tree(V);
    std::vector<bool> visited_tmp(V);

    double ust_lo = std::min(peel_lower, remain_lower);
    double ust_hi = std::max(peel_upper, remain_upper);

    for (int kk = 0; kk < K_est; kk++) {
        clear_tree(ust_tmp);
        int root_tmp;
        int result = sample_sub_ust(g, ust_tmp, V, root_tmp, visited_tmp, ignore,
                                    pop, ust_lo, ust_hi, counties, cg);
        if (result != 0) continue;
        trees_drawn++;

        std::vector<int> pop_below_tmp(V, 0);
        std::vector<int> parent_tmp(V);
        parent_tmp[root_tmp] = -1;
        tree_pop(ust_tmp, root_tmp, pop, pop_below_tmp, parent_tmp);

        for (int v = 0; v < V; v++) {
            if (ignore[v] || v == root_tmp) continue;
            double below = pop_below_tmp[v];
            double above = region_pop - below;
            bool ok = (peel_lower <= below && below <= peel_upper &&
                       remain_lower <= above && above <= remain_upper) ||
                      (remain_lower <= below && below <= remain_upper &&
                       peel_lower <= above && above <= peel_upper);
            if (ok) total_valid++;
        }
    }

    if (trees_drawn == 0) return -23.0; // ~log(1e-10)

    // Return log of expected *count* of valid-cut vertices per tree (= p_hat * (region_size-1)).
    // The log-ratio log_p_fwd[s] - log_p_rev[s] is the paper's approximate correction (eq. 163).
    double q_hat = (double) total_valid / (double) trees_drawn;
    return std::log(std::max(q_hat, 1e-10));
}


/*
 * Main MMSS MCMC loop.
 */
// [[Rcpp::export]]
Rcpp::List mmss_plans(int N, List l, const arma::uvec init, const arma::uvec &counties,
                     const arma::uvec &pop, int n_distr, double target, double lower,
                     double upper, double rho, List constraints, List control,
                     int thin, int l_merge, int verbosity) {
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    int max_retries = 200;
    if (control.containsElementNamed("max_retries")) {
        max_retries = (int) control["max_retries"];
    }
    int K_est = 25;
    if (control.containsElementNamed("k_est")) {
        K_est = (int) control["k_est"];
    }
    bool exact_mh = true;
    if (control.containsElementNamed("exact_mh")) {
        exact_mh = (bool) control["exact_mh"];
    }
    bool valid_cuts_only = !exact_mh; // correct default for each path
    if (control.containsElementNamed("valid_cuts_only")) {
        valid_cuts_only = (bool) control["valid_cuts_only"];
    }

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
        Rcout << "MARKOV CHAIN MONTE CARLO (Multiple Merge Sequential Split, l=" << l_merge << ")\n";
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
    }

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

        bool split_failed = false;
        double fwd_boundary_lp = 0.0;
        double rev_boundary_lp = 0.0;
        double prop_correction = 0.0;

        if (exact_mh) {
        // ====== EXACT PATH: whole-sequence retry ======
        // Retrying the complete split as one unit preserves q ∝ T*T*B.
        // For l_merge == 2 the retry factor cancels exactly, since the success
        // probability depends only on the merged region (same for forward and
        // reverse). For l_merge >= 3 the cancellation is approximate because
        // intermediate subgraph shapes can differ; the residual is corrected
        // below via the region-size factor (lines 629-648).
        split_failed = true;

        for (int attempt = 0; attempt < max_retries; attempt++) {
            districts.col(idx + 1) = saved_plan;
            fwd_boundary_lp = 0.0;
            bool attempt_ok = true;

            for (int s = 0; s < l_merge - 1; s++) {
                int peel = sel_districts[s];
                int remain = sel_districts[l_merge - 1];

                for (int v = 0; v < V; v++) {
                    for (int t = s + 1; t < l_merge - 1; t++) {
                        if ((int) districts(v, idx + 1) == sel_districts[t]) {
                            districts(v, idx + 1) = remain;
                            break;
                        }
                    }
                }

                double region_pop = 0.0;
                for (int v = 0; v < V; v++) {
                    ignore[v] = ((int) districts(v, idx + 1) != peel &&
                                 (int) districts(v, idx + 1) != remain);
                    if (!ignore[v]) region_pop += pop(v);
                }

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

                if (peel_lower > peel_upper) { attempt_ok = false; break; }

                double ust_lower = std::min(peel_lower, remain_lower);
                double ust_upper = std::max(peel_upper, remain_upper);

                clear_tree(ust);
                int root;
                int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                            pop, ust_lower, ust_upper, counties, cg);
                if (result != 0) { attempt_ok = false; break; }

                auto col_ref = districts.col(idx + 1);
                if (!cut_one_mms(ust, root, col_ref,
                                 peel, remain, pop, region_pop,
                                 peel_lower, peel_upper,
                                 remain_lower, remain_upper,
                                 target, /*from_valid_only=*/valid_cuts_only)) {
                    attempt_ok = false;
                    break;
                }

                fwd_boundary_lp += log_boundary(g, districts.col(idx + 1), peel, remain);

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

            if (attempt_ok) { split_failed = false; break; }
        } // end retry loop

        } else {
        // ====== APPROXIMATE PATH: per-step retry with MH correction ======
        // Each split step retries independently until success. For steps s>=1,
        // we estimate the per-step success probability on both forward and
        // reverse subgraphs and include a correction in the MH ratio.
        std::vector<double> log_p_fwd(l_merge - 1, 0.0);

        for (int s = 0; s < l_merge - 1 && !split_failed; s++) {
            int peel = sel_districts[s];
            int remain = sel_districts[l_merge - 1];

            for (int v = 0; v < V; v++) {
                for (int t = s + 1; t < l_merge - 1; t++) {
                    if ((int) districts(v, idx + 1) == sel_districts[t]) {
                        districts(v, idx + 1) = remain;
                        break;
                    }
                }
            }

            double region_pop = 0.0;
            for (int v = 0; v < V; v++) {
                ignore[v] = ((int) districts(v, idx + 1) != peel &&
                             (int) districts(v, idx + 1) != remain);
                if (!ignore[v]) region_pop += pop(v);
            }

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

            if (peel_lower > peel_upper) { split_failed = true; break; }

            double ust_lower = std::min(peel_lower, remain_lower);
            double ust_upper = std::max(peel_upper, remain_upper);

            uvec step_state(V);
            for (int v = 0; v < V; v++) step_state(v) = districts(v, idx + 1);

            bool step_ok = false;
            for (int tree_attempt = 0; tree_attempt < max_retries; tree_attempt++) {
                for (int v = 0; v < V; v++) {
                    if (!ignore[v]) districts(v, idx + 1) = step_state(v);
                }

                clear_tree(ust);
                int root;
                int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                            pop, ust_lower, ust_upper, counties, cg);
                if (result != 0) continue;

                auto col_ref = districts.col(idx + 1);
                if (cut_one_mms(ust, root, col_ref,
                                peel, remain, pop, region_pop,
                                peel_lower, peel_upper,
                                remain_lower, remain_upper,
                                target, /*from_valid_only=*/valid_cuts_only)) {
                    step_ok = true;
                    break;
                }
            }

            if (!step_ok) { split_failed = true; break; }

            fwd_boundary_lp += log_boundary(g, districts.col(idx + 1), peel, remain);

            if (s >= 1) {
                log_p_fwd[s] = log_p_edge_estimate(g, pop, ignore, region_pop,
                                                    peel_lower, peel_upper,
                                                    remain_lower, remain_upper,
                                                    counties, cg, K_est);
            }

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
        } // end per-step loop

        // Compute p_edge correction for the approximate path
        if (!split_failed) {
            std::vector<double> log_p_rev(l_merge - 1, 0.0);
            uvec old_plan = districts.col(idx);
            umat work_mat(V, 1);
            work_mat.col(0) = old_plan;
            for (int v = 0; v < V; v++) {
                for (int d : sel_districts) {
                    if ((int) work_mat(v, 0) == d) {
                        work_mat(v, 0) = 0;
                        break;
                    }
                }
            }
            for (int s = 0; s < l_merge - 1; s++) {
                if (s >= 1) {
                    std::vector<bool> ignore_rev(V, true);
                    double region_pop_rev = 0.0;
                    for (int v = 0; v < V; v++) {
                        if ((int) work_mat(v, 0) == 0) {
                            ignore_rev[v] = false;
                            region_pop_rev += pop(v);
                        }
                    }
                    int remaining_splits = l_merge - 1 - s;
                    double p_lo = std::max(lower, region_pop_rev - remaining_splits * upper);
                    double p_hi = std::min(upper, region_pop_rev - remaining_splits * lower);
                    double r_lo = (remaining_splits > 1) ? remaining_splits * lower : lower;
                    double r_hi = (remaining_splits > 1) ? remaining_splits * upper : upper;
                    log_p_rev[s] = log_p_edge_estimate(g, pop, ignore_rev, region_pop_rev,
                                                        p_lo, p_hi, r_lo, r_hi,
                                                        counties, cg, K_est);
                }
                int dist_label = sel_districts[s];
                for (int v = 0; v < V; v++) {
                    if ((int) old_plan(v) == dist_label && work_mat(v, 0) == 0) {
                        work_mat(v, 0) = dist_label;
                    }
                }
            }
            for (int s = 1; s < l_merge - 1; s++) {
                prop_correction += log_p_fwd[s] - log_p_rev[s];
            }
        }

        } // end approximate path

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

        // 3. Reverse proposal boundary terms (needed for both paths)
        {
            uvec old_plan = districts.col(idx);
            umat work_mat(V, 1);
            work_mat.col(0) = old_plan;
            for (int v = 0; v < V; v++) {
                for (int d : sel_districts) {
                    if ((int) work_mat(v, 0) == d) {
                        work_mat(v, 0) = 0;
                        break;
                    }
                }
            }
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

        // For the exact MH path with l_merge >= 3: apply the region-size correction
        // from eq. 135. For s >= 1, the forward region size
        // |Ṽ_s^fwd| = |R| - sum_{j<s}|V_{d_j}^new| differs from the reverse
        // |Ṽ_s^rev| = |R| - sum_{j<s}|V_{d_j}^old| whenever district sizes change.
        // The correction factor is prod_{s=1}^{l-2} (|Ṽ_s^fwd|-1)/(|Ṽ_s^rev|-1).
        if (exact_mh && l_merge >= 3) {
            int total_region = 0;
            for (int v = 0; v < V; v++) {
                for (int d : sel_districts) {
                    if ((int) saved_plan(v) == d) { total_region++; break; }
                }
            }
            int fwd_cumsize = 0, rev_cumsize = 0;
            for (int s = 1; s < l_merge - 1; s++) {
                int fwd_s = 0, rev_s = 0;
                for (int v = 0; v < V; v++) {
                    if ((int) districts(v, idx + 1) == sel_districts[s - 1]) fwd_s++;
                    if ((int) districts(v, idx)     == sel_districts[s - 1]) rev_s++;
                }
                fwd_cumsize += fwd_s;
                rev_cumsize += rev_s;
                prop_correction += std::log((double)(total_region - fwd_cumsize - 1) /
                                            (double)(total_region - rev_cumsize - 1));
            }
        }

        // MH proposal log-ratio: log q(y->x) - log q(x->y) + path-specific correction
        prop_lp = rev_boundary_lp - fwd_boundary_lp + prop_correction;

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
        if (prop_lp >= 0 || std::log(r_unif()) <= prop_lp) {
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
    out["mhdecisions"] = mh_decisions;

    return out;
}
