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
 * Cut a spanning tree to peel off one district from a two-label region.
 *
 * Two picking strategies via `from_valid_only` and `k_topk`:
 *
 * from_valid_only=TRUE, k_topk=0: pick uniformly from valid cuts only.
 *
 * from_valid_only=TRUE, k_topk>0 (top-k proposal):
 *   Sort all region vertices by |pop_below - target|, pick uniformly from
 *   the top-k_topk, return false if chosen vertex is not a valid cut (caller
 *   redraws the spanning tree).  k_topk = k_seq[s] is fixed from G[R] before
 *   the main loop, so it is identical for forward and reverse proposals and
 *   cancels in the MH ratio (Lemma 1, paper).
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
                        int k_topk,
                        std::vector<int> &pop_below,
                        std::vector<int> &parent,
                        std::vector<int> &region_verts,
                        std::vector<std::pair<double, int>> &dev_verts,
                        const std::vector<int> &iter_region) {

    parent[root] = -1;
    tree_pop(ust, root, pop, pop_below, parent);

    // Collect non-root vertices in the active region (peel or remain)
    region_verts.clear();
    for (int i : iter_region) {
        if (i == root) continue;
        if ((int) districts(i) == peel_label || (int) districts(i) == remain_label)
            region_verts.push_back(i);
    }

    if (region_verts.empty()) return false;

    int cut_at;
    if (from_valid_only) {
        if (k_topk == 0) {
            // Pick uniformly from valid cuts via reservoir sampling (O(1) memory).
            int count = 0;
            cut_at = -1;
            for (int v : region_verts) {
                double below = pop_below[v];
                double above = total_pop - below;
                bool ok = (peel_lower <= below && below <= peel_upper &&
                           remain_lower <= above && above <= remain_upper) ||
                          (remain_lower <= below && below <= remain_upper &&
                           peel_lower <= above && above <= peel_upper);
                if (ok) {
                    count++;
                    if (r_int(count) == 0) cut_at = v;
                }
            }
            n_valid_cuts = count;
            if (count == 0) return false;
        } else {
            // Top-k: use partial_sort on region vertices by |pop_below - target|
            dev_verts.clear();
            int n_valid = 0;
            for (int v : region_verts) {
                double below = pop_below[v];
                double dev = std::min(std::abs(below - peel_target),
                                      std::abs(total_pop - below - peel_target));
                dev_verts.push_back({dev, v});
                double above = total_pop - below;
                bool ok = (peel_lower <= below && below <= peel_upper &&
                           remain_lower <= above && above <= remain_upper) ||
                          (remain_lower <= below && below <= remain_upper &&
                           peel_lower <= above && above <= peel_upper);
                if (ok) n_valid++;
            }
            n_valid_cuts = n_valid;
            int k_actual = std::min(k_topk, (int) dev_verts.size());
            if (k_actual <= 0) return false;
            if (k_actual < (int) dev_verts.size()) {
                std::nth_element(dev_verts.begin(), dev_verts.begin() + k_actual,
                                 dev_verts.end());
            }
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
 * Compute log boundary between peel and remain, scanning only region vertices.
 */
static double log_boundary_region(const Graph &g,
                                  const subview_col<uword> &districts,
                                  int distr_root, int distr_other,
                                  const std::vector<int> &region_verts) {
    double count = 0;
    for (int v : region_verts) {
        if ((int) districts(v) != distr_root) continue;
        for (int nbor : g[v]) {
            if ((int) districts(nbor) == distr_other) count += 1.0;
        }
    }
    return std::log(count);
}


/*
 * Random walk along `g` from `root` until something in `visited` is hit.
 * No county check — for single-county or county-free regions.
 * Uses CSR flat graph and int8_t status array for fast random access.
 * Status: 0=unvisited, 1=in-tree, 2=ignored
 */
static int walk_until_nc(const std::vector<int> &flat_adj,
                         const std::vector<int> &flat_off,
                         int root,
                         std::vector<int> &path, int MAX,
                         std::vector<int8_t> &status) {
    path[0] = root;
    int curr = root;
    int added = 1;
    int i;
    for (i = 0; i < MAX; i++) {
        int off = flat_off[curr];
        int deg = flat_off[curr + 1] - off;
        int proposal = flat_adj[off + r_int(deg)];
        int8_t s = status[proposal];
        if (s == 2) {
            continue;
        } else if (s == 0) {
            for (int j = added - 1; j >= 0; j--) {
                if (path[j] == proposal) {
                    added = j;
                    break;
                }
            }
            path[added++] = proposal;
        } else {
            path[added++] = proposal;
            break;
        }
        curr = proposal;
    }
    if (i == MAX) {
        added = 0;
    }
    return added;
}


/*
 * Simplified UST sampler for no-county regions.
 * Uses CSR flat graph, int8_t status array, and pre-allocated unvisited list.
 */
static int sample_sub_ust_nc(const std::vector<int> &flat_adj,
                             const std::vector<int> &flat_off,
                             Tree &tree, int V, int &root,
                             std::vector<int8_t> &status,
                             const std::vector<bool> &ignore,
                             std::vector<int> &walk_buf,
                             std::vector<int> &unvisited,
                             std::vector<int> &unvis_pos) {
    // Build status and unvisited list from ignore array
    // Status: 0=unvisited, 1=in-tree, 2=ignored
    unvisited.clear();
    for (int i = 0; i < V; i++) {
        if (ignore[i]) {
            status[i] = 2;
        } else {
            status[i] = 0;
            unvis_pos[i] = (int) unvisited.size();
            unvisited.push_back(i);
        }
    }

    int remaining = (int) unvisited.size();
    if (remaining == 0) return 1;

    // Pick root via O(1) random selection
    int pick = r_int(remaining);
    root = unvisited[pick];
    // Swap-remove from unvisited
    int last_v = unvisited[remaining - 1];
    unvisited[pick] = last_v;
    unvis_pos[last_v] = pick;
    unvisited.pop_back();
    status[root] = 1;
    remaining--;

    if (remaining > 0) {
        int max_try = std::max(50, 50 * remaining * ((int) std::log(remaining)));
        while (remaining > 0) {
            // O(1) random unvisited vertex selection
            pick = r_int(remaining);
            int add = unvisited[pick];
            int added = walk_until_nc(flat_adj, flat_off, add, walk_buf,
                                       max_try, status);
            if (added == 0) return 1;
            // Mark in-tree and remove from unvisited list
            for (int i = 0; i < added - 1; i++) {
                int v = walk_buf[i];
                status[v] = 1;
                tree[walk_buf[i+1]].push_back(v);
                // Swap-remove v from unvisited
                int pos = unvis_pos[v];
                remaining--;
                if (remaining > 0 && pos < remaining) {
                    int last = unvisited[remaining];
                    unvisited[pos] = last;
                    unvis_pos[last] = pos;
                }
                unvisited.pop_back();
            }
        }
    }

    return 0;
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
    Graph g= list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int n_cty = max(counties);
    std::vector<uvec> tree_levels;
    bool use_multi_hierarchy = false;
    if (control.containsElementNamed("tree_levels")) {
        IntegerMatrix level_mat = control["tree_levels"];
        for (int j = 0; j < level_mat.ncol(); j++) {
            uvec lvl(V);
            for (int i = 0; i < V; i++) {
                lvl(i) = level_mat(i, j);
            }
            tree_levels.push_back(lvl);
        }
        use_multi_hierarchy = tree_levels.size() > 1;
    }
    std::vector<uvec> county_tree_levels;
    if (!use_multi_hierarchy && n_cty > 1) {
        county_tree_levels.push_back(counties);
    }

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
        if (use_multi_hierarchy) {
            Rcout << "Sampling hierarchically across " << tree_levels.size()
                  << " administrative levels.\n";
        } else if (cg.size() > 1) {
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
        }
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
            if (ks.size() == 1) {
                std::fill(fixed_k_seq.begin(), fixed_k_seq.end(), ks[0]);
            } else {
                for (int s = 0; s < (int) fixed_k_seq.size() && s < ks.size(); s++) {
                    fixed_k_seq[s] = ks[s];
                }
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

    // Diagnostics for the exact path.
    // n_m_hit counts proposals where all max_retries whole-sequence attempts
    // failed — these force staying at the current state and can bias the chain
    // if they happen frequently.
    int n_m_hit = 0;
    long long n_attempts = 0;
    long long n_ust_draws = 0;
    long long n_ust_fail = 0;
    long long n_cut_fail = 0;
    long long n_proposal_success = 0;
    // Per-step valid cut distributions from successful cuts, indexed [s][0..3]
    // where 3 means "3 or more". Records K_T^s (total valid cuts in the drawn
    // UST) whenever cut_one_mms succeeds. If K_T^s > k_seq[s] for any tree,
    // some valid cuts fall outside the top-k window — this is the key
    // assumption-violation diagnostic.
    int n_steps = std::max(l_merge - 1, 1);
    std::vector<std::array<int,4>> n_cuts_dist(n_steps);
    std::vector<int> max_valid_cuts(n_steps, 0);
    for (auto &a : n_cuts_dist) a.fill(0);

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

    // Pre-allocated buffers for cut_one_mms (avoid allocation in retry loop)
    std::vector<int> pop_below_buf(V, 0);
    std::vector<int> parent_buf(V);
    std::vector<int> region_verts_buf;
    region_verts_buf.reserve(V);
    std::vector<std::pair<double, int>> dev_verts_buf;
    dev_verts_buf.reserve(V);

    // Pre-allocated walk buffer for UST sampling (avoids allocation per retry)
    std::vector<int> walk_buf(V + 2);

    // Detect single-county case for fast-path UST sampling
    bool use_fast_ust = (n_cty == 1 && !use_multi_hierarchy);
    dev_verts_buf.reserve(V);

    // CSR flat graph for cache-friendly random walks in no-county UST
    std::vector<int> flat_adj;
    std::vector<int> flat_off(V + 1);
    if (use_fast_ust) {
        int total_edges = 0;
        for (int i = 0; i < V; i++) total_edges += (int) g[i].size();
        flat_adj.reserve(total_edges);
        for (int i = 0; i < V; i++) {
            flat_off[i] = (int) flat_adj.size();
            for (int nbor : g[i]) flat_adj.push_back(nbor);
        }
        flat_off[V] = (int) flat_adj.size();
    }

    // Pre-allocated unvisited vertex list and position index for O(1) selection
    std::vector<int> unvisited_buf;
    unvisited_buf.reserve(V);
    std::vector<int> unvis_pos(V);

    // Pre-allocated status array for fast-path UST (0=unvisited, 1=in-tree, 2=ignored)
    std::vector<int8_t> ust_status(V);

    // Pre-allocated buffer for reverse boundary computation
    std::vector<int> rev_label(V, 0);

    // Pre-allocated buffer for iteration region vertices
    std::vector<int> iter_region;
    iter_region.reserve(V);
    std::vector<int> active_region;
    active_region.reserve(V);

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

        // Pre-collect region vertices (those in any selected district)
        // districts.col(idx) is unchanged during proposal and serves as the
        // "saved plan" for restoring on retry.
        iter_region.clear();
        for (int v = 0; v < V; v++) {
            int lbl = (int) districts(v, idx);
            for (int d : sel_districts) {
                if (lbl == d) { iter_region.push_back(v); break; }
            }
        }
        // Pre-set ignore for all non-region vertices (constant across retries)
        std::fill(ignore.begin(), ignore.end(), true);

        bool split_failed = false;
        double fwd_boundary_lp = 0.0;
        double rev_boundary_lp = 0.0;
        double prop_correction = 0.0;

        // Whole-sequence retry: draw full split sequence, retry from scratch if any step fails.
        split_failed = true;

        for (int attempt = 0; attempt < max_retries; attempt++) {
            n_attempts++;
            // Restore only region vertices (non-region are unchanged)
            for (int v : iter_region) {
                districts(v, idx + 1) = districts(v, idx);
            }
            fwd_boundary_lp = 0.0;
            bool attempt_ok = true;

            for (int s = 0; s < l_merge - 1; s++) {
                int peel = sel_districts[s];
                int remain = sel_districts[l_merge - 1];

                // Merge: relabel intermediate districts to remain (region only)
                // For l_merge=2, no intermediate districts exist — skip.
                if (l_merge > 2) {
                    for (int v : iter_region) {
                        int lbl = (int) districts(v, idx + 1);
                        for (int t = s + 1; t < l_merge - 1; t++) {
                            if (lbl == sel_districts[t]) {
                                districts(v, idx + 1) = remain;
                                break;
                            }
                        }
                    }
                }

                // Set ignore and compute region_pop (region only)
                double region_pop = 0.0;
                active_region.clear();
                for (int v : iter_region) {
                    int lbl = (int) districts(v, idx + 1);
                    bool active = (lbl == peel || lbl == remain);
                    ignore[v] = !active;
                    if (active) {
                        region_pop += pop(v);
                        active_region.push_back(v);
                    }
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

                // Clear tree entries for region vertices only (O(R) instead of O(V))
                for (int v : iter_region) {
                    ust[v].clear();
                }
                int root;
                int result;
                n_ust_draws++;
                if (use_fast_ust) {
                    result = sample_sub_ust_nc(flat_adj, flat_off,
                                               ust, V, root, ust_status,
                                               ignore, walk_buf,
                                               unvisited_buf, unvis_pos);
                } else if (use_multi_hierarchy) {
                    result = sample_sub_ust_hier(g, ust, V, root, visited,
                                                 ignore, active_region,
                                                 pop, ust_lower, ust_upper,
                                                 tree_levels);
                } else if (n_cty > 1) {
                    result = sample_sub_ust_hier(g, ust, V, root, visited,
                                                 ignore, active_region,
                                                 pop, ust_lower, ust_upper,
                                                 county_tree_levels);
                } else {
                    result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                             pop, ust_lower, ust_upper, counties, cg);
                }
                if (result != 0) {
                    n_ust_fail++;
                    attempt_ok = false;
                    break;
                }

                auto col_ref = districts.col(idx + 1);
                int nvc = 0;
                if (!cut_one_mms(ust, root, col_ref,
                                 peel, remain, pop, region_pop,
                                 peel_lower, peel_upper,
                                 remain_lower, remain_upper,
                                 target, /*from_valid_only=*/true,
                                 nvc, fixed_k_seq[s],
                                 pop_below_buf, parent_buf,
                                 region_verts_buf, dev_verts_buf,
                                 iter_region)) {
                    n_cut_fail++;
                    attempt_ok = false;
                    break;
                }

                n_cuts_dist[s][std::min(nvc, 3)]++;
                if (nvc > max_valid_cuts[s]) max_valid_cuts[s] = nvc;

                fwd_boundary_lp += log_boundary_region(g, districts.col(idx + 1), peel, remain, iter_region);

                // Restore: un-merge vertices that were relabeled to remain
                // For l_merge=2, no intermediate districts were merged — skip.
                if (l_merge > 2) {
                    for (int v : iter_region) {
                        if ((int) districts(v, idx + 1) == remain) {
                            int orig = districts(v, idx);
                            for (int t = s + 1; t < l_merge; t++) {
                                if (orig == sel_districts[t]) {
                                    districts(v, idx + 1) = orig;
                                    break;
                                }
                            }
                        }
                    }
                }
            } // end steps loop

            if (attempt_ok) { split_failed = false; break; }
        } // end retry loop

        if (split_failed) {
            n_m_hit++;
            districts.col(idx + 1) = districts.col(idx);
            if (i % thin == 0) idx++;
            continue;
        }
        n_proposal_success++;

        // Ensure all vertices in the selected set are properly assigned.
        // For l=2, cut_one_mms already assigns all vertices via assign_district — skip.
        if (l_merge > 2) {
            int last_label = sel_districts[l_merge - 1];
            for (int v : iter_region) {
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
        }

        // 3. Reverse proposal boundary terms (needed for both paths)
        {
            const uvec &old_plan = districts.unsafe_col(idx);
            for (int v : iter_region) {
                rev_label[v] = 0;
            }
            for (int s = 0; s < l_merge - 1; s++) {
                int dist_label = sel_districts[s];
                for (int v : iter_region) {
                    if ((int) old_plan(v) == dist_label && rev_label[v] == 0) {
                        rev_label[v] = dist_label;
                    }
                }
                // Compute boundary between label 0 and dist_label, region only
                // Note: non-region neighbors always have rev_label = 0 (cleaned up after each iter)
                double count = 0;
                for (int v : iter_region) {
                    if (rev_label[v] != 0) continue;
                    for (int nbor : g[v]) {
                        if (rev_label[nbor] == dist_label) count += 1.0;
                    }
                }
                rev_boundary_lp += std::log(count);
            }
            // Clean up rev_label for region vertices so next iteration starts clean
            for (int v : iter_region) {
                rev_label[v] = 0;
            }
        }

        // For l_merge >= 3: apply the region-size correction
        if (l_merge >= 3) {
            int total_region = (int) iter_region.size();
            int fwd_cumsize = 0, rev_cumsize = 0;
            for (int s = 1; s < l_merge - 1; s++) {
                int fwd_s = 0, rev_s = 0;
                for (int v : iter_region) {
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
        new_dist_g = district_graph(g, districts.col(idx + 1), n_distr);

        if (l_merge == 2) {
            // Fast path for l=2: equivalent to MS pair selection correction.
            // log_prob_perm([d1,d2]) = -log(n) - log(deg(d1))
            // log_prob_perm([d2,d1]) = -log(n) - log(deg(d2))
            // log_psum = -log(n) + log(1/deg(d1) + 1/deg(d2))
            // The -log(n) cancels in the difference.
            int d1 = sel_districts[0] - 1;
            int d2 = sel_districts[1] - 1;
            prop_lp -= std::log(1.0 / dist_g[d1].size() + 1.0 / dist_g[d2].size());
            prop_lp += std::log(1.0 / new_dist_g[d1].size() + 1.0 / new_dist_g[d2].size());
        } else {
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
        }

        // 7. Accept/reject
        if (prop_lp >= 0 || std::log(r_unif()) <= prop_lp) {
            n_accept++;
            districts.col(idx) = districts.col(idx + 1);
            dist_g = std::move(new_dist_g);
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
        if (n_m_hit > 0) {
            Rcout << "WARNING: " << n_m_hit
                  << " proposal(s) exhausted max_retries=" << max_retries
                  << " without a valid whole-sequence split. Consider increasing max_retries.\n";
        }
        if (l_merge > 1) {
            for (int s = 0; s < n_steps; s++) {
                long long n_s = n_cuts_dist[s][0] + n_cuts_dist[s][1] +
                                n_cuts_dist[s][2] + n_cuts_dist[s][3];
                if (n_s > 0) {
                    Rcout << "Valid cuts in successful USTs (s=" << s << ", k_seq=" << fixed_k_seq[s] << "): "
                          << "1=" << n_cuts_dist[s][1]
                          << ", 2=" << n_cuts_dist[s][2]
                          << ", 3+=" << n_cuts_dist[s][3]
                          << "; max=" << max_valid_cuts[s];
                    if (max_valid_cuts[s] > fixed_k_seq[s])
                        Rcout << " [EXCEEDS k_seq!]";
                    Rcout << ".\n";
                }
            }
        }
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["mhdecisions"] = mh_decisions;
    out["n_m_hit"] = n_m_hit;
    out["n_attempts"] = (double) n_attempts;
    out["n_ust_draws"] = (double) n_ust_draws;
    out["n_ust_fail"] = (double) n_ust_fail;
    out["n_cut_fail"] = (double) n_cut_fail;
    out["n_proposal_success"] = (double) n_proposal_success;
    // Per-step valid cut distribution matrix: rows = steps, cols = [1, 2, 3+, max]
    // "0" column omitted since cut_one_mms only records on success (nvc >= 1)
    Rcpp::IntegerMatrix valid_cuts_mat(n_steps, 4);
    Rcpp::CharacterVector col_names = {"1", "2", "3+", "max"};
    Rcpp::CharacterVector row_names(n_steps);
    for (int s = 0; s < n_steps; s++) {
        row_names[s] = "s" + std::to_string(s);
        valid_cuts_mat(s, 0) = n_cuts_dist[s][1];
        valid_cuts_mat(s, 1) = n_cuts_dist[s][2];
        valid_cuts_mat(s, 2) = n_cuts_dist[s][3];
        valid_cuts_mat(s, 3) = max_valid_cuts[s];
    }
    Rcpp::colnames(valid_cuts_mat) = col_names;
    Rcpp::rownames(valid_cuts_mat) = row_names;
    out["valid_cuts_by_step"] = valid_cuts_mat;
    // Fixed top-k sequence (same for all proposals; k_s cancels in MH ratio)
    out["k_seq"] = Rcpp::IntegerVector(fixed_k_seq.begin(), fixed_k_seq.end());

    return out;
}
