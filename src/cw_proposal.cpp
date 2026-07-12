#include "cw_proposal.h"

#include <algorithm>
#include <cmath>
#include <map>

// ---------------------------------------------------------------------
// Subtree population helpers (operate on the LCT)
// ---------------------------------------------------------------------

static int topological_sort_helper(std::map<int, int> &cut_pop, LCTNode *node, LCTNode *source,
                                   arma::uvec const &pop, bool reversed, int mass);

static std::map<int, int> compute_subtree_pops(LCTGraphPlan &plan, arma::uvec const &pop,
                                               int root_vertex) {
    std::map<int, int> cut_pop;
    LinkCutTree &lct = *plan.lct;

    lct.evert(root_vertex);
    LCTNode *root = lct.node(root_vertex);

    bool rev = root->reversed;
    int lc = rev ? 1 : 0;
    int rc = 1 - lc;

    int total = 0;
    if (root->children[rc] != nullptr) {
        total += topological_sort_helper(cut_pop, root->children[rc], root, pop, rev, 0);
    }
    for (LCTNode *child : root->path_children) {
        total += topological_sort_helper(cut_pop, child, root, pop, false, 0);
    }

    int root_pop = pop(root_vertex);
    cut_pop[root_vertex] = root_pop + total;
    return cut_pop;
}

static int topological_sort_helper(std::map<int, int> &cut_pop, LCTNode *node, LCTNode *source,
                                   arma::uvec const &pop, bool reversed, int mass) {
    if (node == nullptr)
        return 0;

    int remainder = 0;
    reversed = reversed ^ node->reversed;
    int lc = reversed ? 1 : 0;
    int rc = 1 - lc;

    if (node->children[rc] != nullptr) {
        remainder +=
            topological_sort_helper(cut_pop, node->children[rc], node, pop, reversed, mass);
    }
    for (LCTNode *child : node->path_children) {
        remainder += topological_sort_helper(cut_pop, child, node, pop, false, 0);
    }

    int node_pop = pop(node->vertex);
    cut_pop[node->vertex] = remainder + node_pop + mass;

    if (node->children[lc] != nullptr && node->children[lc] != source) {
        remainder += topological_sort_helper(cut_pop, node->children[lc], node, pop, reversed,
                                             cut_pop[node->vertex]);
    }

    return remainder + node_pop;
}

// ---------------------------------------------------------------------
// Step 1: pick a random adjacent district pair
// ---------------------------------------------------------------------

bool get_random_adjacent_districts(LCTGraphPlan const &plan, RNGState &rng_state, int &d1,
                                   int &d2) {
    auto pairs = plan.get_adjacent_district_pairs();
    if (pairs.empty())
        return false;
    int idx = rng_state.r_int(static_cast<int>(pairs.size()));
    d1 = pairs[idx].first;
    d2 = pairs[idx].second;
    return true;
}

// ---------------------------------------------------------------------
// Step 2: pick two distinct boundary edges between d1 and d2
// ---------------------------------------------------------------------

bool get_random_edge_pair(LCTGraphPlan const &plan, int d1, int d2, RNGState &rng_state,
                          CWEdge &e1, CWEdge &e2) {
    EdgeSet const &edges = plan.get_cross_edges(d1, d2);
    int n_edges = static_cast<int>(edges.size());
    if (n_edges < 2)
        return false;

    int idx1 = rng_state.r_int(n_edges);
    int idx2 = rng_state.r_int(n_edges - 1);
    if (idx2 >= idx1)
        idx2++;

    auto it = edges.begin();
    std::advance(it, idx1);
    e1 = *it;
    it = edges.begin();
    std::advance(it, idx2);
    e2 = *it;
    return true;
}

// ---------------------------------------------------------------------
// Step 3: compute the cycle paths in each district
// ---------------------------------------------------------------------

bool get_cycle_paths(LCTGraphPlan &plan, CWEdge const &e1, CWEdge const &e2,
                     std::vector<int> &path1, std::vector<int> &path2) {
    int u1 = e1.u, v1 = e1.v;
    int u2 = e2.u, v2 = e2.v;

    int d1 = plan.region_ids[u1];
    if (plan.region_ids[u1] != d1)
        std::swap(u1, v1);
    if (plan.region_ids[u2] != d1)
        std::swap(u2, v2);

    LinkCutTree &lct = *plan.lct;

    lct.evert(u1);
    path1 = lct.find_path(u2);
    std::reverse(path1.begin(), path1.end());

    lct.evert(v1);
    path2 = lct.find_path(v2);
    std::reverse(path2.begin(), path2.end());

    return !path1.empty() && !path2.empty();
}

// ---------------------------------------------------------------------
// Step 4: collapsed cycle weights
// ---------------------------------------------------------------------

std::vector<int> get_collapsed_cycle_weights(LCTGraphPlan &plan, std::vector<int> const &path1,
                                             std::vector<int> const &path2,
                                             arma::uvec const &pop) {
    int cycle_len = static_cast<int>(path1.size() + path2.size());
    std::vector<int> collapsed_weights(cycle_len);

    int u1 = path1.back();
    std::map<int, int> u_cut_pop = compute_subtree_pops(plan, pop, u1);

    int v1 = path2.back();
    std::map<int, int> v_cut_pop = compute_subtree_pops(plan, pop, v1);

    for (size_t ii = 0; ii < path1.size(); ii++) {
        int vertex = path1[ii];
        collapsed_weights[ii] = u_cut_pop[vertex];
        if (ii > 0) {
            collapsed_weights[ii] -= u_cut_pop[path1[ii - 1]];
        }
    }

    for (size_t ii = 0; ii < path2.size(); ii++) {
        int pos = cycle_len - 1 - static_cast<int>(ii);
        int vertex = path2[ii];
        collapsed_weights[pos] = v_cut_pop[vertex];
        if (ii > 0) {
            collapsed_weights[pos] -= v_cut_pop[path2[ii - 1]];
        }
    }

    return collapsed_weights;
}

// ---------------------------------------------------------------------
// Step 5: enumerate valid (cut1, cut2) pairs
// ---------------------------------------------------------------------

std::vector<std::pair<int, int>> find_valid_cut_pairs(std::vector<int> const &cycle_pops,
                                                      int initial_cut, int total_pop,
                                                      double lower1, double upper1,
                                                      double lower2, double upper2) {
    std::vector<std::pair<int, int>> valid_pairs;
    int n = static_cast<int>(cycle_pops.size());

    std::vector<int> prefix(n + 1, 0);
    for (int i = 0; i < n; i++)
        prefix[i + 1] = prefix[i] + cycle_pops[i];

    for (int cut1 = 1; cut1 <= n; cut1++) {
        for (int cut2 = cut1; cut2 <= n - 1; cut2++) {
            if (cut1 == 1 && cut2 == initial_cut)
                continue; // identity

            int pop1 = prefix[cut2] - prefix[cut1 - 1];
            int pop2 = total_pop - pop1;

            // Accept if EITHER assignment of (pop1,pop2) to (d1,d2) is valid.
            bool ok_direct =
                (pop1 >= lower1 && pop1 <= upper1) && (pop2 >= lower2 && pop2 <= upper2);
            bool ok_swap =
                (pop1 >= lower2 && pop1 <= upper2) && (pop2 >= lower1 && pop2 <= upper1);
            if (ok_direct || ok_swap)
                valid_pairs.push_back({cut1, cut2});
        }
    }
    return valid_pairs;
}

// ---------------------------------------------------------------------
// Link-vertex position helpers (matches Julia's get_link_path_ind)
// ---------------------------------------------------------------------

static int get_link_path_ind(int link_vertex, std::vector<int> const &path1,
                             std::vector<int> const &path2) {
    int path1_len = static_cast<int>(path1.size());
    int path2_len = static_cast<int>(path2.size());

    if (path1[path1_len - 1] == link_vertex)
        return 1;
    if (path1[0] == link_vertex)
        return path1_len;
    if (path2[path2_len - 1] == link_vertex)
        return path1_len + path2_len;
    if (path2[0] == link_vertex)
        return path1_len + 1;

    throw std::runtime_error("Couldn't find link vertex in paths");
}

static bool swap_assignment_check(int path_ind, int cut1, int cut2,
                                  std::vector<int> const &path1, std::vector<int> const &path2,
                                  std::vector<int> const &cycle_weights) {
    int path1_len = static_cast<int>(path1.size());
    int cycle_len = static_cast<int>(cycle_weights.size());

    int overlap1 = 0;
    int tot_pop = 0;
    for (int w : cycle_weights)
        tot_pop += w;

    if (cut1 <= path1_len) {
        for (int i = cut1; i <= std::min(path1_len, cut2); i++)
            overlap1 += cycle_weights[i - 1];
    } else if (cut1 > path1_len + 1) {
        for (int i = path1_len + 1; i < cut1; i++)
            overlap1 += cycle_weights[i - 1];
    }
    if (cut2 < cycle_len) {
        for (int i = std::max(path1_len + 1, cut2 + 1); i <= cycle_len; i++)
            overlap1 += cycle_weights[i - 1];
    }

    bool uPathToInterval = (2 * overlap1 > tot_pop);
    bool l11_in_interval = (cut1 <= path_ind && path_ind <= cut2);
    bool l11_in_uPath = (path_ind <= path1_len);

    return (l11_in_uPath != l11_in_interval) == uPathToInterval;
}

// ---------------------------------------------------------------------
// Apply an update to plan/LCT
// ---------------------------------------------------------------------

void apply_update(LCTGraphPlan &plan, MapParams const &map_params,
                  CycleWalkUpdate const &update) {
    if (!update.valid)
        return;

    LinkCutTree &lct = *plan.lct;

    for (auto const &cut : update.cuts)
        lct.cut(cut.second);
    for (auto const &link : update.links) {
        lct.evert(link.first);
        lct.link(link.first, link.second);
    }

    int d1 = update.changed_districts.first;
    int d2 = update.changed_districts.second;

    if (!update.cuts.empty()) {
        int new_root1 = lct.find_root(update.cuts[0].first);
        int new_root2 = lct.find_root(update.cuts[0].second);

        if (!update.links.empty()) {
            int link11_vertex = update.links[0].first;
            int link11_dist_cur = plan.region_ids[link11_vertex];
            int r11_new_root = lct.find_root(link11_vertex);
            int r11_root_ind_new = (r11_new_root != new_root1) ? 2 : 1;
            int r11_root_ind_cur = (link11_dist_cur != d1) ? 2 : 1;
            bool should_swap = (r11_root_ind_new == r11_root_ind_cur) == update.swap_link11;
            if (should_swap)
                std::swap(new_root1, new_root2);
        }

        plan.district_roots[d1] = new_root1;
        plan.district_roots[d2] = new_root2;
    }

    // Reassign vertices by LCT root.
    int const V = map_params.V;
    for (int v = 0; v < V; v++) {
        int root = lct.find_root(v);
        if (root == plan.district_roots[d1]) {
            plan.region_ids[v] = static_cast<RegionID>(d1);
        } else if (root == plan.district_roots[d2]) {
            plan.region_ids[v] = static_cast<RegionID>(d2);
        }
    }

    // Recompute populations for the two affected regions.
    plan.region_pops[d1] = 0;
    plan.region_pops[d2] = 0;
    for (int v = 0; v < V; v++) {
        int d = plan.region_ids[v];
        if (d == d1)
            plan.region_pops[d1] += map_params.pop(v);
        else if (d == d2)
            plan.region_pops[d2] += map_params.pop(v);
    }

    plan.rebuild_cross_edges(map_params.g);
}

// ---------------------------------------------------------------------
// Main entry: cycle_walk
// ---------------------------------------------------------------------

int cycle_walk(LCTGraphPlan &plan, MapParams const &map_params,
               ScoringFunction const &scoring_function, double const compactness,
               RNGState &rng_state, CycleWalkDiagnostics &diagnostics) {
    diagnostics = CycleWalkDiagnostics();

    int d1, d2;
    if (!get_random_adjacent_districts(plan, rng_state, d1, d2)) {
        diagnostics.status = -1;
        return -1;
    }

    CWEdge e1(0, 0), e2(0, 0);
    if (!get_random_edge_pair(plan, d1, d2, rng_state, e1, e2)) {
        diagnostics.status = -2;
        return -2;
    }

    std::vector<int> path1, path2;
    if (!get_cycle_paths(plan, e1, e2, path1, path2)) {
        diagnostics.status = -3;
        return -3;
    }

    std::vector<int> cycle_pops =
        get_collapsed_cycle_weights(plan, path1, path2, map_params.pop);
    diagnostics.cycle_length = static_cast<int>(cycle_pops.size());

    int total_pop = plan.region_pops[d1] + plan.region_pops[d2];
    int cycle_pop_sum = 0;
    for (int p : cycle_pops)
        cycle_pop_sum += p;
    if (cycle_pop_sum != total_pop) {
        Rcpp::Rcout << "[ERROR] Cycle pop sum mismatch: sum=" << cycle_pop_sum
                    << ", total=" << total_pop << "\n";
    }

    int initial_cut = static_cast<int>(path1.size());

    int s1 = plan.region_sizes[d1];
    int s2 = plan.region_sizes[d2];
    double lo1 = map_params.lower * s1;
    double up1 = map_params.upper * s1;
    double lo2 = map_params.lower * s2;
    double up2 = map_params.upper * s2;
    std::vector<std::pair<int, int>> valid_pairs =
        find_valid_cut_pairs(cycle_pops, initial_cut, total_pop, lo1, up1, lo2, up2);
    diagnostics.n_valid_cuts = static_cast<int>(valid_pairs.size());

    if (valid_pairs.empty()) {
        diagnostics.status = -4;
        plan.lct->evert(plan.district_roots[d1]);
        plan.lct->evert(plan.district_roots[d2]);
        return -4;
    }

    int n_valid_pairs_fwd = static_cast<int>(valid_pairs.size());

    auto get_edge_at_position = [&](int edge_ind) -> std::pair<int, int> {
        int path1_len = static_cast<int>(path1.size());
        int path2_len = static_cast<int>(path2.size());
        if (edge_ind == 1)
            return {path1[0], path2[0]};
        if (edge_ind <= path1_len)
            return {path1[edge_ind - 1], path1[edge_ind - 2]};
        if (edge_ind == path1_len + 1)
            return {path1[path1_len - 1], path2[path2_len - 1]};
        int ind = edge_ind - path1_len - 1;
        return {path2[path2_len - ind], path2[path2_len - ind - 1]};
    };

    std::vector<double> cum_edge_weight_product(n_valid_pairs_fwd);
    double cumsum = 0.0;
    for (int i = 0; i < n_valid_pairs_fwd; i++) {
        auto [c1, c2] = valid_pairs[i];
        auto [e1u, e1v] = get_edge_at_position(c1);
        auto [e2u, e2v] = get_edge_at_position(c2 + 1);
        double w1 = plan.get_edge_weight(e1u, e1v);
        double w2 = plan.get_edge_weight(e2u, e2v);
        cumsum += 1.0 / (w1 * w2);
        cum_edge_weight_product[i] = cumsum;
    }

    double rand_samp = rng_state.r_unif() * cum_edge_weight_product[n_valid_pairs_fwd - 1];
    int sample_idx = 0;
    while (sample_idx < n_valid_pairs_fwd - 1 &&
           rand_samp > cum_edge_weight_product[sample_idx]) {
        sample_idx++;
    }
    auto [cut1, cut2] = valid_pairs[sample_idx];

    auto [sel_e1u, sel_e1v] = get_edge_at_position(cut1);
    auto [sel_e2u, sel_e2v] = get_edge_at_position(cut2 + 1);
    double w1_cuts = plan.get_edge_weight(sel_e1u, sel_e1v);
    double w2_cuts = plan.get_edge_weight(sel_e2u, sel_e2v);
    double w1w2_cuts_inv = 1.0 / (w1_cuts * w2_cuts);

    double w1_links = plan.get_edge_weight(e1.u, e1.v);
    double w2_links = plan.get_edge_weight(e2.u, e2.v);
    double w1w2_links_inv = 1.0 / (w1_links * w2_links);

    double sum_edge_weight_products = cum_edge_weight_product[n_valid_pairs_fwd - 1];

    int old_boundary = static_cast<int>(plan.get_cross_edges(d1, d2).size());

    // ---- snapshot state for revert ----
    std::vector<RegionID> old_region_ids(plan.region_ids.begin(), plan.region_ids.end());
    std::vector<int> old_region_pops(plan.region_pops.begin(), plan.region_pops.end());
    std::vector<int> old_district_roots = plan.district_roots;
    CrossEdgeMap old_cross_edges = plan.cross_edges;

    // ---- scoring on the OLD state ----
    bool const is_final = true; // cyclewalk plans are always fully districted
    double old_soft_score_d1 = scoring_function.compute_region_soft_score(plan, d1);
    double old_soft_score_d2 = scoring_function.compute_region_soft_score(plan, d2);
    auto old_plan_score = scoring_function.compute_plan_score(plan);
    double old_log_st = 0.0;
    if (compactness != 1.0) {
        old_log_st = plan.compute_log_region_spanning_trees(map_params, d1) +
                     plan.compute_log_region_spanning_trees(map_params, d2);
    }

    // ---- construct update ----
    CycleWalkUpdate update;
    update.changed_districts = {d1, d2};
    update.valid = true;

    auto [fe1u, fe1v] = get_edge_at_position(cut1);
    auto [fe2u, fe2v] = get_edge_at_position(cut2 + 1);

    int e1_d1, e1_d2, e2_d1, e2_d2;
    if (plan.region_ids[e1.u] == d1) {
        e1_d1 = e1.u;
        e1_d2 = e1.v;
    } else {
        e1_d1 = e1.v;
        e1_d2 = e1.u;
    }
    if (plan.region_ids[e2.u] == d1) {
        e2_d1 = e2.u;
        e2_d2 = e2.v;
    } else {
        e2_d1 = e2.v;
        e2_d2 = e2.u;
    }

    std::vector<std::pair<int, int>> links = {{e1_d1, e1_d2}, {e2_d1, e2_d2}};
    std::vector<std::pair<int, int>> cuts = {{fe1u, fe1v}, {fe2u, fe2v}};

    auto edges_equal = [](std::pair<int, int> a, std::pair<int, int> b) {
        return (a.first == b.first && a.second == b.second) ||
               (a.first == b.second && a.second == b.first);
    };

    if (edges_equal(links[0], cuts[0]) || edges_equal(links[0], cuts[1])) {
        if (edges_equal(links[0], cuts[0]))
            cuts.erase(cuts.begin());
        else
            cuts.erase(cuts.begin() + 1);
        links.erase(links.begin());
    } else if (edges_equal(links[1], cuts[0]) || edges_equal(links[1], cuts[1])) {
        if (edges_equal(links[1], cuts[0]))
            cuts.erase(cuts.begin());
        else
            cuts.erase(cuts.begin() + 1);
        links.erase(links.begin() + 1);
    }

    update.cuts = cuts;
    update.links = links;

    if (!links.empty()) {
        int link11_vertex = links[0].first;
        int path_ind_l11 = get_link_path_ind(link11_vertex, path1, path2);
        update.swap_link11 =
            swap_assignment_check(path_ind_l11, cut1, cut2, path1, path2, cycle_pops);
    } else {
        update.swap_link11 = false;
    }

    apply_update(plan, map_params, update);

    // ---- scoring on the NEW state ----
    double new_soft_score_d1 = scoring_function.compute_region_soft_score(plan, d1);
    double new_soft_score_d2 = scoring_function.compute_region_soft_score(plan, d2);
    auto new_plan_score = scoring_function.compute_plan_score(plan);
    double new_log_st = 0.0;
    if (compactness != 1.0) {
        new_log_st = plan.compute_log_region_spanning_trees(map_params, d1) +
                     plan.compute_log_region_spanning_trees(map_params, d2);
    }

    int new_boundary = static_cast<int>(plan.get_cross_edges(d1, d2).size());

    auto count_adj_dists_involving = [](CrossEdgeMap const &cross_edges, int a, int b) {
        int count = 0;
        for (auto const &[key, edges] : cross_edges) {
            if (key.first == a || key.first == b || key.second == a || key.second == b) {
                if (!edges.empty())
                    count++;
            }
        }
        return count;
    };

    int old_adj_dists_d1d2 = count_adj_dists_involving(old_cross_edges, d1, d2);
    int new_adj_dists_d1d2 = count_adj_dists_involving(plan.cross_edges, d1, d2);
    int delta_adj_dists = new_adj_dists_d1d2 - old_adj_dists_d1d2;

    int old_adj_dists_total = 0;
    for (auto const &[key, edges] : old_cross_edges) {
        if (!edges.empty())
            old_adj_dists_total++;
    }

    // ---- log-MH ratio ----
    bool new_state_ok = old_plan_score.first && new_plan_score.first;
    double log_mh_ratio = 0.0;

    if (old_adj_dists_total > 0 && old_adj_dists_total + delta_adj_dists > 0) {
        log_mh_ratio += std::log(static_cast<double>(old_adj_dists_total)) -
                        std::log(static_cast<double>(old_adj_dists_total + delta_adj_dists));
    }
    if (new_boundary > 0 && old_boundary > 0) {
        log_mh_ratio += std::log(static_cast<double>(old_boundary * (old_boundary - 1))) -
                        std::log(static_cast<double>(new_boundary * (new_boundary - 1)));
    }
    double new_sum_weights = sum_edge_weight_products + w1w2_links_inv - w1w2_cuts_inv;
    if (sum_edge_weight_products > 0 && new_sum_weights > 0) {
        log_mh_ratio += std::log(sum_edge_weight_products) - std::log(new_sum_weights);
    }
    // Soft constraint ratio: target ∝ exp(-penalty), so log ratio = old - new.
    double old_soft_total = old_soft_score_d1 + old_soft_score_d2 + old_plan_score.second;
    double new_soft_total = new_soft_score_d1 + new_soft_score_d2 + new_plan_score.second;
    log_mh_ratio += old_soft_total - new_soft_total;

    if (compactness != 1.0) {
        log_mh_ratio += (1.0 - compactness) * (old_log_st - new_log_st);
    }

    double accept_ratio = new_state_ok ? std::min(1.0, std::exp(log_mh_ratio)) : 0.0;
    diagnostics.accept_prob = accept_ratio;

    double rand_val = rng_state.r_unif();
    bool do_accept = rand_val < accept_ratio;

    if (do_accept) {
        diagnostics.status = 1;
        return 1;
    }

    // Reject: undo LCT changes (reverse order), restore plan snapshot.
    diagnostics.status = 0;
    LinkCutTree &lct = *plan.lct;
    for (auto const &link : update.links) {
        lct.evert(link.second);
        lct.cut(link.first);
    }
    for (auto const &cut : update.cuts) {
        lct.evert(cut.first);
        lct.link(cut.first, cut.second);
    }

    std::copy(old_region_ids.begin(), old_region_ids.end(), plan.region_ids.begin());
    std::copy(old_region_pops.begin(), old_region_pops.end(), plan.region_pops.begin());
    plan.district_roots = old_district_roots;
    plan.cross_edges = old_cross_edges;

    lct.evert(old_district_roots[d1]);
    lct.evert(old_district_roots[d2]);
    return 0;
}
