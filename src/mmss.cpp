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


static double total_pop_of(const std::vector<int> &verts, const uvec &pop) {
    double total = 0.0;
    for (int v : verts) total += pop(v);
    return total;
}


static void collect_subtree_vertices(const Tree &tree, int root,
                                     std::vector<int> &verts) {
    std::vector<int> stack = {root};
    while (!stack.empty()) {
        int v = stack.back();
        stack.pop_back();
        verts.push_back(v);
        for (int child : tree[v]) stack.push_back(child);
    }
}


static Tree restrict_tree_to_vertices(const Tree &tree,
                                      const std::vector<bool> &keep) {
    int V = tree.size();
    Tree restricted = init_tree(V);
    for (int v = 0; v < V; v++) {
        if (!keep[v]) continue;
        for (int child : tree[v]) {
            if (keep[child]) restricted[v].push_back(child);
        }
    }
    return restricted;
}


/*
 * Estimate the top-k sequence for the exact path.  Draws K_est spanning trees
 * on G[R], walks each through the split sequence by cutting at the best valid
 * vertex and restricting to the remainder.  For each step s, also draws one
 * fresh G̃_s UST to count valid cuts in that subgraph directly — this gives
 * the correct distribution for k_seq[s] since the actual forward proposal
 * draws fresh G̃_s USTs, not restricted G[R]-trees.
 */
static std::vector<int> estimate_k_sequence(const Graph &g, const uvec &pop,
                                            const std::vector<int> &region_verts,
                                            int l_merge, int K_est,
                                            double total_pop_full,
                                            int n_districts, double pop_dev,
                                            const uvec &counties,
                                            Multigraph &cg) {
    std::vector<int> k_seq(l_merge - 1, 1);
    if (l_merge <= 1 || K_est <= 0 || region_verts.size() <= 1) return k_seq;

    int V = g.size();
    double pbar = total_pop_full / n_districts;
    double region_total = total_pop_of(region_verts, pop);

    std::vector<bool> ignore(V, true);
    std::vector<bool> visited(V, false);
    for (int v : region_verts) ignore[v] = false;

    for (int kk = 0; kk < K_est; kk++) {
        Tree current_tree = init_tree(V);
        int current_root;
        int result = sample_sub_ust(g, current_tree, V, current_root, visited, ignore,
                                    pop, 0.0, region_total, counties, cg);
        if (result != 0) continue;

        std::vector<int> current_region = region_verts;
        double current_pop = region_total;
        bool tree_alive = true;

        for (int s = 0; s < l_merge - 1 && tree_alive; s++) {
            int remaining_splits = l_merge - 1 - s;
            double peel_lower = std::max(pbar * (1 - pop_dev),
                                         current_pop - remaining_splits * pbar * (1 + pop_dev));
            double peel_upper = std::min(pbar * (1 + pop_dev),
                                         current_pop - remaining_splits * pbar * (1 - pop_dev));
            if (peel_lower > peel_upper) { tree_alive = false; break; }

            std::vector<int> pop_below(V, 0);
            std::vector<int> parent(V, -1);
            tree_pop(current_tree, current_root, pop, pop_below, parent);

            int best_v = -1;
            double best_dev = std::numeric_limits<double>::infinity();
            for (int v : current_region) {
                if (v == current_root) continue;
                double below = pop_below[v];
                double above = current_pop - below;
                bool ok = (peel_lower <= below && below <= peel_upper) ||
                          (peel_lower <= above && above <= peel_upper);
                if (!ok) continue;
                double dev = std::min(std::abs(below - pbar), std::abs(above - pbar));
                if (dev < best_dev) { best_dev = dev; best_v = v; }
            }

            if (best_v < 0) { tree_alive = false; break; }

            // Draw a fresh G̃_s UST to count valid cuts from the correct distribution.
            // The restricted G[R]-tree has a different distribution than the fresh
            // G̃_s USTs that cut_one_mms actually draws in the forward proposal.
            {
                std::vector<bool> ignore_s(V, true);
                for (int v : current_region) ignore_s[v] = false;

                Tree fresh_tree = init_tree(V);
                std::vector<bool> fresh_visited(V, false);
                int fresh_root;
                int fresh_result = sample_sub_ust(g, fresh_tree, V, fresh_root,
                                                  fresh_visited, ignore_s, pop,
                                                  0.0, current_pop, counties, cg);
                if (fresh_result == 0) {
                    std::vector<int> fresh_pop_below(V, 0);
                    std::vector<int> fresh_parent(V, -1);
                    tree_pop(fresh_tree, fresh_root, pop, fresh_pop_below, fresh_parent);

                    int fresh_valid = 0;
                    for (int v : current_region) {
                        if (v == fresh_root) continue;
                        double below = fresh_pop_below[v];
                        double above = current_pop - below;
                        bool ok = (peel_lower <= below && below <= peel_upper) ||
                                  (peel_lower <= above && above <= peel_upper);
                        if (ok) fresh_valid++;
                    }
                    k_seq[s] = std::max(k_seq[s], fresh_valid);
                }
            }

            double below_best = pop_below[best_v];
            double above_best = current_pop - below_best;
            bool peel_is_below = (peel_lower <= below_best && below_best <= peel_upper);

            std::vector<int> subtree_vertices;
            collect_subtree_vertices(current_tree, best_v, subtree_vertices);

            std::vector<bool> keep(V, false);
            if (peel_is_below) {
                for (int v : current_region) keep[v] = true;
                for (int v : subtree_vertices) keep[v] = false;
                current_pop = above_best;
            } else {
                for (int v : subtree_vertices) keep[v] = true;
                current_root = best_v;
                current_pop = below_best;
            }

            current_tree = restrict_tree_to_vertices(current_tree, keep);
            current_region.clear();
            current_region.reserve(V);
            for (int v = 0; v < V; v++) {
                if (keep[v]) current_region.push_back(v);
            }
            if (current_region.size() <= 1) tree_alive = false;
        }
    }

    return k_seq;
}


/*
 * Cut a spanning tree to peel off one district from a two-label region.
 *
 * Two picking strategies via `from_valid_only` and `k_topk`:
 *
 * from_valid_only=TRUE, k_topk=0: pick uniformly from valid cuts only.
 *
 * from_valid_only=TRUE, k_topk>0 (both exact and approximate paths):
 *   Sort all region vertices by |pop_below - target|, pick uniformly from
 *   the top-k_topk, return false if chosen vertex is not a valid cut (caller
 *   redraws the spanning tree).  Both paths use k_topk = k_seq[s] estimated
 *   from G[R] via estimate_k_sequence, so k_topk is identical for forward
 *   and reverse proposals and cancels in the MH ratio (Lemma 1, paper).
 *
 * from_valid_only=FALSE: pick a uniformly random non-root vertex; return
 *   false if not a valid cut.
 */
static bool cut_one_mms(Tree &ust, int root,
                        subview_col<uword> &districts,
                        int peel_label, int remain_label,
                        const uvec &pop, double total_pop,
                        double peel_lower, double peel_upper,
                        double remain_lower, double remain_upper,
                        double peel_target,
                        bool from_valid_only,
                        int &n_valid_cuts,
                        int k_topk = 0) {
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
        // Enumerate valid cuts (needed by both branches for n_valid_cuts and approx path)
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
        n_valid_cuts = (int) valid_verts.size();

        if (k_topk == 0) {
            // Approximate path (s=0): pick uniformly from valid cuts.
            if (valid_verts.empty()) return false;
            cut_at = valid_verts[r_int(valid_verts.size())];
        } else {
            // Top-k (exact path, or approximate path at s>=1): sort all
            // region vertices by |pop_below - target|, pick uniformly from
            // the top-k, fail (caller redraws tree) if chosen vertex is
            // not a valid cut.  With k_topk >= K_s = max_T[k_T], all valid
            // cuts of any tree fall within the top-k window and this never
            // fails when k_T > 0.
            std::vector<std::pair<double, int>> dev_verts;
            dev_verts.reserve(region_verts.size());
            for (int v : region_verts) {
                double below = pop_below[v];
                double dev = std::min(std::abs(below - peel_target),
                                      std::abs(total_pop - below - peel_target));
                dev_verts.push_back({dev, v});
            }
            std::sort(dev_verts.begin(), dev_verts.end());
            int k_actual = std::min(k_topk, (int) dev_verts.size());
            if (k_actual <= 0) return false;
            cut_at = dev_verts[r_int(k_actual)].second;
            double below = pop_below[cut_at];
            double above = total_pop - below;
            bool ok = (peel_lower <= below && below <= peel_upper &&
                       remain_lower <= above && above <= remain_upper) ||
                      (remain_lower <= below && below <= remain_upper &&
                       peel_lower <= above && above <= peel_upper);
            if (!ok) return false;
        }
    } else {
        n_valid_cuts = 0;
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
    bool exact_mh = true;
    if (control.containsElementNamed("exact_mh")) {
        exact_mh = (bool) control["exact_mh"];
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
    double total_pop = arma::accu(pop);

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

    // Fixed top-k sequence ONCE before the main loop.
    // k_seq must be identical for forward and reverse proposals for reversibility,
    // so it must NOT be estimated from the selected region R (which varies per proposal).
    // Default: k_seq[s] = l_merge - s (top-l at step 0, top-2 at last step).
    // User can override by passing k_seq as an integer vector in control.
    std::vector<int> fixed_k_seq(std::max(l_merge - 1, 0), 1);
    if (l_merge > 1) {
        if (control.containsElementNamed("k_seq")) {
            Rcpp::IntegerVector ks = control["k_seq"];
            for (int s = 0; s < (int) fixed_k_seq.size() && s < ks.size(); s++) {
                fixed_k_seq[s] = ks[s];
            }
        } else {
            for (int s = 0; s < (int) fixed_k_seq.size(); s++) {
                fixed_k_seq[s] = l_merge - s;
            }
        }
        if (verbosity >= 1) {
            Rcout << "Fixed top-k sequence (k_seq):";
            for (int s = 0; s < (int) fixed_k_seq.size(); s++) {
                Rcout << " " << fixed_k_seq[s];
            }
            Rcout << "\n";
        }
    }

    int n_accept = 0;
    // Diagnostics for approximate path.
    // n_cuts_dist[k] counts successfully-sampled USTs with exactly k valid cuts
    // (k=3 means 3 or more). Index 0 captures trees where cut_one_mms failed
    // (no valid cuts found), i.e. the rejected attempts before a winning tree.
    int n_m_hit = 0;
    // Step s=0 spans the full merged region G[R]; valid-cut multiplicity
    // cancels between forward and reverse, so it doesn't bias the chain.
    int n_cuts_dist_s0[4] = {0, 0, 0, 0};
    int max_valid_cuts_s0 = 0;
    long long n_valid_trees_s0 = 0, n_valid_sum_s0 = 0;
    // Steps s>=1 span peeled subregions; extra valid cuts here create bias
    // because the single-valid-cut assumption may not hold for these trees.
    int n_cuts_dist_s1[4] = {0, 0, 0, 0};
    int max_valid_cuts_s1 = 0;
    long long n_valid_trees_s1 = 0, n_valid_sum_s1 = 0;

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
                int dummy_nvalid = 0;
                if (!cut_one_mms(ust, root, col_ref,
                                 peel, remain, pop, region_pop,
                                 peel_lower, peel_upper,
                                 remain_lower, remain_upper,
                                 target, /*from_valid_only=*/true,
                                 dummy_nvalid, fixed_k_seq[s])) {
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
        // ====== APPROXIMATE PATH: per-step retry, pre-fixed top-k from G[R] ======
        // Uses the same k_seq estimated from G[R] as the exact path, so k_s is
        // identical for forward and reverse proposals and cancels in the MH ratio.
        // The approximation is that the per-step retry success probability is
        // symmetric (E[k_T | G̃_s^fwd] ≈ E[k_T | G̃_s^rev]), validated empirically.
        // No MH correction term is added (prop_correction stays 0 for this path).

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
                int nvc = 0;
                if (cut_one_mms(ust, root, col_ref,
                                peel, remain, pop, region_pop,
                                peel_lower, peel_upper,
                                remain_lower, remain_upper,
                                target, /*from_valid_only=*/true,
                                nvc, fixed_k_seq[s])) {
                    step_ok = true;
                    if (s == 0) {
                        n_cuts_dist_s0[std::min(nvc, 3)]++;
                        n_valid_trees_s0++;
                        n_valid_sum_s0 += nvc;
                        if (nvc > max_valid_cuts_s0) max_valid_cuts_s0 = nvc;
                    } else {
                        n_cuts_dist_s1[std::min(nvc, 3)]++;
                        n_valid_trees_s1++;
                        n_valid_sum_s1 += nvc;
                        if (nvc > max_valid_cuts_s1) max_valid_cuts_s1 = nvc;
                    }
                    break;
                } else {
                    // Record nvc: 0 = no valid cuts; >0 = k_seq underestimated K_s
                    if (s == 0) n_cuts_dist_s0[std::min(nvc, 3)]++;
                    else        n_cuts_dist_s1[std::min(nvc, 3)]++;
                }
            }

            if (!step_ok) { n_m_hit++; split_failed = true; break; }

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
        } // end per-step loop

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
        // |Γ_s^fwd| = |R| - sum_{j<s}|V_{d_j}^new| differs from the reverse
        // |Γ_s^rev| = |R| - sum_{j<s}|V_{d_j}^old| whenever district sizes change.
        // The correction factor is prod_{s=1}^{l-2} (|Γ_s^fwd|-1)/(|Γ_s^rev|-1).
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
        if (!exact_mh) {
            if (n_m_hit > 0) {
                Rcout << "WARNING: " << n_m_hit
                      << " proposal step(s) exhausted max_retries=" << max_retries
                      << ". Consider increasing max_retries or k_est.\n";
            }
            long long n_s0 = n_cuts_dist_s0[0] + n_cuts_dist_s0[1] +
                             n_cuts_dist_s0[2] + n_cuts_dist_s0[3];
            long long n_s1 = n_cuts_dist_s1[0] + n_cuts_dist_s1[1] +
                             n_cuts_dist_s1[2] + n_cuts_dist_s1[3];
            if (n_s0 > 0) {
                double mean_vc0 = n_valid_trees_s0 > 0
                    ? (double) n_valid_sum_s0 / n_valid_trees_s0 : 0.0;
                Rcout << std::setprecision(2)
                      << "Valid cuts in proposal trees (s=0, full region; cancels in MH): "
                      << "0=" << n_cuts_dist_s0[0]
                      << ", 1=" << n_cuts_dist_s0[1]
                      << ", 2=" << n_cuts_dist_s0[2]
                      << ", 3+=" << n_cuts_dist_s0[3]
                      << "; mean=" << mean_vc0
                      << ", max=" << max_valid_cuts_s0 << ".\n";
            }
            if (n_s1 > 0) {
                double mean_vc1 = n_valid_trees_s1 > 0
                    ? (double) n_valid_sum_s1 / n_valid_trees_s1 : 0.0;
                Rcout << std::setprecision(2)
                      << "Valid cuts in proposal trees (s>=1, subregions): "
                      << "0=" << n_cuts_dist_s1[0]
                      << ", 1=" << n_cuts_dist_s1[1]
                      << ", 2=" << n_cuts_dist_s1[2]
                      << ", 3+=" << n_cuts_dist_s1[3]
                      << "; mean=" << mean_vc1
                      << ", max=" << max_valid_cuts_s1 << ".\n";
            }
        }
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["mhdecisions"] = mh_decisions;
    out["n_m_hit"] = n_m_hit;
    out["max_valid_cuts_s0"] = max_valid_cuts_s0;
    out["max_valid_cuts_s1"] = max_valid_cuts_s1;
    out["mean_valid_cuts_s0"] = n_valid_trees_s0 > 0 ? (double) n_valid_sum_s0 / n_valid_trees_s0 : 0.0;
    out["mean_valid_cuts_s1"] = n_valid_trees_s1 > 0 ? (double) n_valid_sum_s1 / n_valid_trees_s1 : 0.0;
    out["valid_cuts_dist_s0"] = Rcpp::IntegerVector::create(
        Rcpp::Named("0")    = n_cuts_dist_s0[0],
        Rcpp::Named("1")    = n_cuts_dist_s0[1],
        Rcpp::Named("2")    = n_cuts_dist_s0[2],
        Rcpp::Named("3+")   = n_cuts_dist_s0[3]
    );
    out["valid_cuts_dist_s1"] = Rcpp::IntegerVector::create(
        Rcpp::Named("0")    = n_cuts_dist_s1[0],
        Rcpp::Named("1")    = n_cuts_dist_s1[1],
        Rcpp::Named("2")    = n_cuts_dist_s1[2],
        Rcpp::Named("3+")   = n_cuts_dist_s1[3]
    );
    // Fixed top-k sequence (same for all proposals; k_s cancels in MH ratio)
    out["k_seq"] = Rcpp::IntegerVector(fixed_k_seq.begin(), fixed_k_seq.end());

    return out;
}
