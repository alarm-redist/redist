/********************************************************
 * BUD MCMC Redistricting Sampler
 * BUD Proposal Implementation
 *
 * Faithful translation of BUD.jl/src/proposals/balanced_up_down_walk.jl
 * and BUD.jl/src/proposals/update.jl
 ********************************************************/

#include "bud_proposal.h"
#include "cw_forest_walk.h"
#include "wilson.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <queue>

using namespace Rcpp;

// Diagnostic: max sub_dists to allow (0 = no limit)
static int g_max_sub_dists = 0;  // TEMP: 0 = no limit
// Helper: safe_link that checks same_tree before linking
static void safe_link(LinkCutTree& lct, int u, int v, const char* /*site*/) {
    if (lct.same_tree(u, v)) {
        return; // skip — already in the same tree
    }
    lct.evert(u);
    lct.link(u, v);
}

// Helper: check if u and v are actual neighbors (directly connected) in the LCT
// This is the C++ equivalent of Julia's neighbors_huh
static bool neighbors_in_lct(LinkCutTree& lct, int u, int v) {
    if (!lct.same_tree(u, v)) return false;
    lct.evert(u);
    lct.expose(v);
    auto path = lct.find_path(v);
    return path.size() == 2;
}

// ============================================================
// 1. get_edge_for_cycle: pick random non-tree, non-marked edge
//    Julia: get_edge_for_cycle in balanced_up_down_walk.jl:86
// ============================================================
static bool get_edge_for_cycle(BUDPartition& partition,
                                int& u_out, int& v_out,
                                bool& same_district) {
    const Graph& g = *(partition.graph);
    int V = partition.n_vertices;

    // Build edge list
    std::vector<std::pair<int, int>> all_edges;
    for (int i = 0; i < V; i++) {
        for (int j : g[i]) {
            if (j > i) all_edges.push_back({i, j});
        }
    }
    if (all_edges.empty()) return false;

    // Build set of marked edge node pairs for fast lookup
    std::set<std::pair<int, int>> marked_node_edges;
    for (auto& [dp, me] : partition.marked_edges) {
        int a = std::min(me.first, me.second);
        int b = std::max(me.first, me.second);
        marked_node_edges.insert({a, b});
    }

    for (int attempt = 0; attempt < 500; attempt++) {
        int idx = r_int((int)all_edges.size());
        int u = all_edges[idx].first;
        int v = all_edges[idx].second;

        // Check not in spanning tree (not LCT neighbors)
        bool are_neighbors = false;
        partition.lct.evert(u);
        std::vector<int> p = partition.lct.find_path(v);
        if (p.size() == 2 && p[0] == u && p[1] == v) {
            are_neighbors = true;
        }

        if (are_neighbors) continue;

        // Check not a marked edge
        int a = std::min(u, v), b = std::max(u, v);
        if (marked_node_edges.count({a, b})) continue;

        u_out = u;
        v_out = v;
        same_district = (partition.node_to_district[u] == partition.node_to_district[v]);
        return true;
    }
    return false;
}

// ============================================================
// 2. internal_forest_walk_edge: forest walk using given edge
//    Julia: internal_forest_walk!(partition.lcp, rng; edge=(u, v))
//    Also updates forest_edges for BUDPartition.
// ============================================================
static int internal_forest_walk_edge(BUDPartition& partition, int u, int v) {
    LinkCutTree& lct = partition.lct;

    // evert u, find path to v
    lct.evert(u);
    std::vector<int> path = lct.find_path(v);

    // If path length is 2, edge is already in tree
    if (path.size() <= 2) {
        int district = partition.node_to_district[v];
        partition.district_roots[district] = lct.find_root(v);
        return 0;
    }

    // Compute cumulative weights for path edges + new edge
    int path_len = (int)path.size();
    std::vector<double> cum_weights(path_len);
    cum_weights[0] = 0.0;
    double cum = 0.0;
    for (int i = 1; i < path_len; i++) {
        double w = partition.get_edge_weight(path[i - 1], path[i]);
        cum += 1.0 / w;
        cum_weights[i] = cum;
    }
    double new_w = partition.get_edge_weight(u, v);
    double total = cum + 1.0 / new_w;

    double r = r_unif() * total;
    if (r > cum_weights[path_len - 1]) {
        // Selected new edge - no change
        int district = partition.node_to_district[v];
        partition.district_roots[district] = lct.find_root(v);
        return 0;
    }

    // Find which edge to cut
    int edge_idx = -1;
    for (int i = 0; i < path_len - 1; i++) {
        if (r > cum_weights[i] && r <= cum_weights[i + 1]) {
            edge_idx = i;
            break;
        }
    }
    if (edge_idx < 0) return 1;

    // Cut and link
    int cut_parent = path[edge_idx];
    int cut_child = path[edge_idx + 1];
    lct.cut(cut_child);
    safe_link(lct, u, v, "forest_walk_edge");

    // Update forest_edges: remove cut edge, add new edge
    auto& fe = partition.forest_edges;
    fe.erase(std::remove_if(fe.begin(), fe.end(), [&](const std::pair<int,int>& e) {
        return (e.first == cut_parent && e.second == cut_child) ||
               (e.first == cut_child && e.second == cut_parent);
    }), fe.end());
    fe.push_back({u, v});

    int district = partition.node_to_district[v];
    partition.district_roots[district] = lct.find_root(v);
    return 0;
}

// ============================================================
// 3. link_district_path: link marked edges in main LCT
//    Julia: link_district_path! in balanced_up_down_walk.jl:104
// ============================================================
static void link_district_path(BUDPartition& partition,
                                const std::vector<std::pair<int,int>>& district_path) {
    for (auto& dp : district_path) {
        DistrictPair key = {std::min(dp.first, dp.second), std::max(dp.first, dp.second)};
        auto it = partition.marked_edges.find(key);
        if (it == partition.marked_edges.end()) continue;
        int n1 = it->second.first;
        int n2 = it->second.second;
        if (partition.lct.same_tree(n1, n2)) continue; // skip if already connected
        partition.lct.evert(n1);
        partition.lct.link(n1, n2);
    }
}

// ============================================================
// 4. cut_district_path: cut marked edges, update roots
//    Julia: cut_district_path! in balanced_up_down_walk.jl:122
// ============================================================
static void cut_district_path(BUDPartition& partition,
                               const std::vector<std::pair<int,int>>& district_path) {
    std::set<int> dists;
    for (auto& dp : district_path) {
        dists.insert(dp.first);
        dists.insert(dp.second);
        DistrictPair key = {std::min(dp.first, dp.second), std::max(dp.first, dp.second)};
        auto it = partition.marked_edges.find(key);
        if (it == partition.marked_edges.end()) continue;
        int n1 = it->second.first;
        int n2 = it->second.second;
        // Only cut if actually neighbors (may not have been linked)
        if (neighbors_in_lct(partition.lct, n1, n2)) {
            partition.lct.evert(n1);
            partition.lct.cut(n2);
        }
    }

    // Fix district roots
    for (int d : dists) {
        int root = partition.district_roots[d];
        int new_root = partition.lct.find_root(root);
        partition.district_roots[d] = new_root;
    }
}

// ============================================================
// 5. compute_collapsed_pops: fragment populations for cycle path
//    Julia: get_collapsed_cycle_weights! in balanced_up_down_walk.jl:171
//    Also computes fragment_map: for each vertex in affected districts,
//    which cycle_path vertex is its fragment root.
// ============================================================
static std::vector<int> compute_collapsed_pops(
    BUDPartition& partition,
    const std::vector<int>& path,
    std::map<int, int>& fragment_map
) {
    int V = partition.n_vertices;
    int n = (int)path.size();
    if (n == 0) return {};

    // After evert(path[0]), cut all consecutive path edges
    partition.lct.evert(path[0]);
    for (int i = n - 1; i >= 1; i--) {
        partition.lct.cut(path[i]);
    }

    // Build set of path vertices for quick lookup
    std::set<int> path_set(path.begin(), path.end());

    // Count population per fragment and build fragment_map
    // Each vertex belongs to exactly one path node's fragment
    std::vector<int> fragment_pop(V, 0);
    fragment_map.clear();
    for (int v = 0; v < V; v++) {
        int root = partition.lct.find_root(v);
        fragment_pop[root] += partition.pop_values[v];
        // Map this vertex to the path vertex that is its fragment root
        // (root might not be a path vertex if it's deep in a subtree,
        //  but actually after cutting path edges, each component has
        //  exactly one path vertex as root or somewhere in it)
        fragment_map[v] = root;
    }

    // Collect collapsed pops indexed by vertex
    std::vector<int> collapsed_pop(V, 0);
    for (int i = 0; i < n; i++) {
        int root = partition.lct.find_root(path[i]);
        collapsed_pop[path[i]] = fragment_pop[root];
    }

    // Re-link all path edges
    for (int i = 1; i < n; i++) {
        safe_link(partition.lct, path[i - 1], path[i], "collapsed_pops_relink");
    }

    return collapsed_pop;
}

// ============================================================
// 6. setup_cycle_lct: build tmp_cycle from path with collapsed pops
//    Julia: part of get_collapsed_cycle_weights!
// ============================================================
static std::vector<int> setup_cycle_lct(
    BUDPartition& partition,
    const std::vector<int>& path,
    const std::vector<int>& collapsed_pop
) {
    LinkCutTree& tmp = partition.tmp_cycle_lct;
    int n = (int)path.size();

    // Clear nodes on path
    for (int v : path) {
        LCTNode* node = tmp.node(v);
        node->parent = nullptr;
        node->children[0] = nullptr;
        node->children[1] = nullptr;
        node->reversed = false;
        node->path_parent = nullptr;
        node->path_children.clear();
    }

    // Build chain: path[0] -> path[1] -> ... -> path[n-1]
    // Link in reverse so path[0] is root
    for (int i = n - 2; i >= 0; i--) {
        safe_link(tmp, path[i + 1], path[i], "setup_cycle_lct");
    }

    // Return cycle_path = path (same order, rooted at path[0])
    tmp.evert(path[0]);
    tmp.expose(path[n - 1]);
    return tmp.find_path(path[n - 1]);
}

// ============================================================
// 7. cycle_cut_and_link: cut/link on tmp cycle LCT
//    Julia: cut_and_link! in tree_operations.jl:260
// ============================================================
static void cycle_cut_and_link(
    BUDPartition& partition,
    int cut_to,       // child node to cut
    int link_from,    // root of its subtree (after cut)
    int link_to,      // node to link to
    const std::vector<int>& collapsed_pop
) {
    LinkCutTree& tmp = partition.tmp_cycle_lct;

    // Make link_from the root before cutting (Julia: make_root!(cut[1]))
    tmp.evert(link_from);

    // Cut
    tmp.cut(cut_to);

    // Link (link_from should be root of its tree after the cut)
    safe_link(tmp, link_from, link_to, "cycle_cut_and_link");

    // Recompute cuttable info for the whole cycle tree (it's small)
    int root = tmp.find_root(link_to);
    set_cuttable_tree_info(tmp, root, collapsed_pop,
                           partition.pop_lower, partition.pop_upper,
                           partition.cuttable_info);
}

// ============================================================
// ============================================================
// 8. get_balanced_cuts: find valid balanced cut positions
//    Julia: getBalancedCuts in balanced_up_down_walk.jl:2
// ============================================================
static BalancedCutsResult get_balanced_cuts(
    BUDPartition& partition,
    const std::vector<int>& cycle_path,
    int sub_dists,
    const std::vector<int>& collapsed_pop
) {
    BalancedCutsResult result;
    result.cumWeight = 0;
    result.pathWeights.push_back(0);

    int n = (int)cycle_path.size();

    // Padded path: [cycle_path[n-1], cycle_path[0], ..., cycle_path[n-1], cycle_path[0]]
    std::vector<int> ext_path;
    ext_path.push_back(cycle_path[n - 1]);
    for (int v : cycle_path) ext_path.push_back(v);
    ext_path.push_back(cycle_path[0]);

    // Initial cuttable info on cycle
    LinkCutTree& tmp = partition.tmp_cycle_lct;
    int root = tmp.find_root(cycle_path[0]);
    set_cuttable_tree_info(tmp, root, collapsed_pop,
                           partition.pop_lower, partition.pop_upper,
                           partition.cuttable_info);

    // Walk around cycle
    int ext_len = (int)ext_path.size();
    for (int ii = 1; ii <= n; ii++) {
        // Julia: ii goes from 2 to length(ext_path)-1 (1-based)
        // C++: ii goes from 1 to n (0-based into ext_path at ii)
        int cut_to = ext_path[ii + 1];     // ext_path[ii+1]
        int link_from = ext_path[ii];      // ext_path[ii]
        int link_to = ext_path[ii - 1];    // ext_path[ii-1]

        cycle_cut_and_link(partition, cut_to, link_from, link_to, collapsed_pop);

        Rcpp::checkUserInterrupt();

        // Check if cuttable
        int new_root = tmp.find_root(cut_to);
        bool closure_says = cuttable_huh(partition.cuttable_info.node_closures[new_root], sub_dists);

        if (closure_says) {
            int edge_ind = ii - 1; // 0-based index into cycle_path
            result.edge_inds.push_back(edge_ind);

            // Weight = 1/edge_weight
            double w = partition.get_edge_weight(cycle_path[edge_ind],
                                                  cycle_path[(edge_ind + 1) % n]);
            result.cumWeight += 1.0 / w;
            result.pathWeights.push_back(result.cumWeight);
        }
    }

    return result;
}

// ============================================================
// 9. find_km1_cuttable_positions: find cut positions where first
//    part is balanced and rest is (k-1)-cuttable
//    Julia: inside get_path_cut_log_prob! and mark_edges!
// ============================================================
static std::vector<int> find_km1_cuttable_positions(
    CuttableInfo& ci,
    const std::vector<int>& edge_inds,
    const std::vector<int>& cycle_path,
    const std::vector<int>& collapsed_pop,
    double lower, double upper,
    int num_dists
) {
    std::vector<int> result;
    int cur_pop = 0;
    int cur_ind = 0;

    for (int i = 0; i < (int)edge_inds.size(); i++) {
        int ei = edge_inds[i];
        for (int j = cur_ind; j <= ei; j++) {
            cur_pop += collapsed_pop[cycle_path[j]];
        }

        int next_v = cycle_path[ei + 1];
        bool has_closure = (next_v >= 0 && next_v < (int)ci.node_closures.size() &&
                           !ci.node_closures[next_v].empty());
        bool is_cuttable = has_closure &&
                           cuttable_huh(ci.node_closures[next_v], num_dists - 1);

        if ((double)cur_pop > upper) break;

        if ((double)cur_pop >= lower && is_cuttable) {
            result.push_back(ei);
        }

        cur_ind = ei + 1;
    }

    return result;
}

// ============================================================
// 10. get_path_cut_log_prob_recursive: backward probability
//     Julia: get_path_cut_log_prob! in balanced_up_down_walk.jl:245
//     C++ uses copies, so we work in local coordinates.
// ============================================================
static double get_path_cut_log_prob_recursive(
    BUDPartition& partition,
    std::vector<int> edge_inds,
    std::vector<int> marked_edge_inds,
    std::vector<int> cycle_path,
    std::vector<int> node_perm,
    const std::vector<int>& collapsed_pop,
    double lower, double upper,
    int num_dists
) {
    if (num_dists == 1) return 0.0;

    int n = (int)cycle_path.size();
    bool reversed = (node_perm[0] > node_perm[n - 1]);

    if (reversed) {
        for (auto& ei : edge_inds) ei = (n - 2) - ei;
        std::reverse(edge_inds.begin(), edge_inds.end());
        for (auto& mei : marked_edge_inds) mei = (n - 2) - mei;
        std::reverse(marked_edge_inds.begin(), marked_edge_inds.end());
        std::reverse(cycle_path.begin(), cycle_path.end());
        std::reverse(node_perm.begin(), node_perm.end());

        partition.tmp_cycle_lct.evert(cycle_path[0]);
        int root = partition.tmp_cycle_lct.find_root(cycle_path[0]);
        set_cuttable_tree_info(partition.tmp_cycle_lct, root, collapsed_pop,
                               lower, upper, partition.cuttable_info);
    }

    // Find km1-cuttable positions
    std::vector<int> km1_positions = find_km1_cuttable_positions(
        partition.cuttable_info, edge_inds, cycle_path,
        collapsed_pop, lower, upper, num_dists);

    if (km1_positions.empty()) return -1e18;

    int marked_ei = marked_edge_inds[0];

    // The marked edge must be in km1_positions

    int low_edge_ind = -1;
    for (int i = 0; i < (int)edge_inds.size(); i++) {
        if (edge_inds[i] == marked_ei) { low_edge_ind = i; break; }
    }
    if (low_edge_ind < 0) return -1e18;

    // Cut the edge (Julia: just cut, no evert or recompute — existing closures valid)
    int cut_v1 = cycle_path[marked_ei];
    int cut_v2 = cycle_path[marked_ei + 1];
    partition.tmp_cycle_lct.cut(cut_v2);

    // Build suffix, shifted to local coordinates
    std::vector<int> r_edge_inds;
    for (int i = low_edge_ind + 1; i < (int)edge_inds.size(); i++) {
        r_edge_inds.push_back(edge_inds[i] - (marked_ei + 1));
    }
    std::vector<int> r_marked_edge_inds;
    for (int i = 1; i < (int)marked_edge_inds.size(); i++) {
        r_marked_edge_inds.push_back(marked_edge_inds[i] - (marked_ei + 1));
    }
    std::vector<int> r_cycle_path(cycle_path.begin() + marked_ei + 1, cycle_path.end());
    std::vector<int> r_node_perm(node_perm.begin() + marked_ei + 1, node_perm.end());

    double log_prob = get_path_cut_log_prob_recursive(
        partition, r_edge_inds, r_marked_edge_inds, r_cycle_path, r_node_perm,
        collapsed_pop, lower, upper, num_dists - 1);

    // Re-link (Julia: evert!(cut_edge[1]); link!(cut_edge[1], cut_edge[2]))
    safe_link(partition.tmp_cycle_lct, cut_v1, cut_v2, "get_path_cut_log_prob_relink");

    // Julia only undoes array reversals (using views). Since C++ uses copies,
    // no array undo is needed. No evert or recompute in Julia's undo step.

    log_prob += -std::log((double)km1_positions.size());
    return log_prob;
}

// ============================================================
// 11. get_path_cut_log_prob: wrapper
//     Julia: get_path_cut_log_prob in balanced_up_down_walk.jl:225
// ============================================================
static double get_path_cut_log_prob(
    BUDPartition& partition,
    const std::vector<int>& edge_inds,
    const std::vector<int>& marked_edge_inds,
    const std::vector<int>& cycle_path,
    const std::vector<int>& node_perm,
    const std::vector<int>& collapsed_pop,
    double lower, double upper,
    int num_dists
) {
    std::vector<int> v_edge_inds(edge_inds.begin(), edge_inds.end() - 1);

    double log_prob = get_path_cut_log_prob_recursive(
        partition, v_edge_inds, std::vector<int>(marked_edge_inds),
        std::vector<int>(cycle_path), std::vector<int>(node_perm),
        collapsed_pop, lower, upper, num_dists);

    // Restore cuttable info
    partition.tmp_cycle_lct.evert(cycle_path[0]);
    int root = partition.tmp_cycle_lct.find_root(cycle_path[0]);
    set_cuttable_tree_info(partition.tmp_cycle_lct, root, collapsed_pop,
                           lower, upper, partition.cuttable_info);

    return log_prob;
}

// ============================================================
// 12. mark_edges_recursive: forward sampling
//     Julia: mark_edges! in balanced_up_down_walk.jl:354
//     C++ uses copies (not views), so we work in local coordinates
//     without accumulating offsets. Edge positions are relative to
//     the current cycle_path.
// ============================================================
static std::pair<std::vector<int>, double> mark_edges_recursive(
    BUDPartition& partition,
    std::vector<int> edge_inds,
    std::vector<int> cycle_path,
    std::vector<int> node_perm,
    const std::vector<int>& collapsed_pop,
    double lower, double upper,
    int num_dists
) {
    if (num_dists == 1) {
        return {{}, 0.0};
    }

    int n = (int)cycle_path.size();
    bool reversed = (node_perm[0] > node_perm[n - 1]);

    if (reversed) {
        for (auto& ei : edge_inds) ei = (n - 2) - ei;
        std::reverse(edge_inds.begin(), edge_inds.end());
        std::reverse(cycle_path.begin(), cycle_path.end());
        std::reverse(node_perm.begin(), node_perm.end());

        partition.tmp_cycle_lct.evert(cycle_path[0]);
        int root = partition.tmp_cycle_lct.find_root(cycle_path[0]);
        set_cuttable_tree_info(partition.tmp_cycle_lct, root, collapsed_pop,
                               lower, upper, partition.cuttable_info);
    }

    // Find km1-cuttable positions
    std::vector<int> km1_positions = find_km1_cuttable_positions(
        partition.cuttable_info, edge_inds, cycle_path,
        collapsed_pop, lower, upper, num_dists);

    if (km1_positions.empty()) {
        return {{}, -1e18};
    }

    // Sample uniformly
    int sample_idx = r_int((int)km1_positions.size());
    int marked_edge = km1_positions[sample_idx];

    // Find position in edge_inds
    int low_edge_ind = -1;
    for (int i = 0; i < (int)edge_inds.size(); i++) {
        if (edge_inds[i] == marked_edge) { low_edge_ind = i; break; }
    }

    // Cut (Julia: just cut, no evert or recompute)
    int cut_v1 = cycle_path[marked_edge];
    int cut_v2 = cycle_path[marked_edge + 1];
    partition.tmp_cycle_lct.cut(cut_v2);

    // Build suffix edge_inds, shifted to be local to the suffix
    std::vector<int> r_edge_inds;
    for (int i = low_edge_ind + 1; i < (int)edge_inds.size(); i++) {
        r_edge_inds.push_back(edge_inds[i] - (marked_edge + 1));
    }
    std::vector<int> r_cycle_path(cycle_path.begin() + marked_edge + 1, cycle_path.end());
    std::vector<int> r_node_perm(node_perm.begin() + marked_edge + 1, node_perm.end());

    auto [child_marked, child_log_prob] = mark_edges_recursive(
        partition, r_edge_inds, r_cycle_path, r_node_perm,
        collapsed_pop, lower, upper, num_dists - 1);

    // Map child results back to parent's local coordinates
    for (auto& mei : child_marked) mei += (marked_edge + 1);
    child_marked.push_back(marked_edge);

    // Re-link
    safe_link(partition.tmp_cycle_lct, cut_v1, cut_v2, "mark_edges_recursive_relink");

    // Undo reversal
    if (reversed) {
        for (auto& mei : child_marked) mei = (n - 2) - mei;
        std::reverse(child_marked.begin(), child_marked.end());
    }

    double log_prob = child_log_prob + (-std::log((double)km1_positions.size()));
    return {child_marked, log_prob};
}

// ============================================================
// 13. mark_edges_forward: wrapper
//     Julia: mark_edges in balanced_up_down_walk.jl:334
// ============================================================
static std::pair<std::vector<int>, double> mark_edges_forward(
    BUDPartition& partition,
    const std::vector<int>& edge_inds,
    const std::vector<int>& cycle_path,
    const std::vector<int>& node_perm,
    const std::vector<int>& collapsed_pop,
    double lower, double upper,
    int num_dists
) {
    std::vector<int> v_edge_inds(edge_inds.begin(), edge_inds.end() - 1);

    auto [marked, log_prob] = mark_edges_recursive(
        partition, v_edge_inds, std::vector<int>(cycle_path),
        std::vector<int>(node_perm), collapsed_pop,
        lower, upper, num_dists);

    // Restore cuttable info
    partition.tmp_cycle_lct.evert(cycle_path[0]);
    int root = partition.tmp_cycle_lct.find_root(cycle_path[0]);
    set_cuttable_tree_info(partition.tmp_cycle_lct, root, collapsed_pop,
                           lower, upper, partition.cuttable_info);

    return {marked, log_prob};
}

// ============================================================
// 14. swap_marked_edge: trivial swap
//     Julia: swap_marked_edge! in balanced_marked_tree.jl:229
// ============================================================
static void swap_marked_edge(BUDPartition& partition,
                              int cut_from, int cut_to,
                              int link_from, int link_to) {
    int d1n = partition.node_to_district[cut_from];
    int d2n = partition.node_to_district[cut_to];
    int d1 = partition.node_to_district[link_from];
    int d2 = partition.node_to_district[link_to];

    // Remove old marked edge
    DistrictPair old_dp = {std::min(d1n, d2n), std::max(d1n, d2n)};
    partition.marked_edges.erase(old_dp);

    // Add new marked edge
    DistrictPair new_dp = {std::min(d1, d2), std::max(d1, d2)};
    partition.marked_edges[new_dp] = {link_from, link_to};

    // Update district tree
    partition.district_tree[d1n].erase(d2n);
    partition.district_tree[d2n].erase(d1n);
    partition.district_tree[d1].insert(d2);
    partition.district_tree[d2].insert(d1);
}

// ============================================================
// 15. bud_proposal: main proposal function
//     Julia: balanced_up_down_walk! in balanced_up_down_walk.jl:477
//     Returns: (proposal_probability, update)
//     If prob == 0, changes are already applied (forest walk / trivial)
// ============================================================
static std::pair<double, BUDUpdate> bud_proposal(
    BUDPartition& partition,
    double lower, double upper
) {
    BUDUpdate update;

    // Pick random non-tree, non-marked edge
    int u, v;
    bool same_district;
    if (!get_edge_for_cycle(partition, u, v, same_district)) {
        return {0.0, update};
    }

    // Same district → internal forest walk
    if (same_district) {
        internal_forest_walk_edge(partition, u, v);
        return {0.0, update};
    }

    int d1 = partition.node_to_district[u];
    int d2 = partition.node_to_district[v];

    // Get district path
    std::vector<std::pair<int,int>> district_path = partition.get_district_path(d1, d2);
    if (district_path.empty()) {
        return {0.0, update};
    }
    int sub_dists = (int)district_path.size() + 1;

    // Diagnostic: skip proposals with too many sub-districts
    if (g_max_sub_dists > 0 && sub_dists > g_max_sub_dists) {
        return {0.0, update};
    }

    Rcpp::checkUserInterrupt();

    // Link marked edges in LCT
    link_district_path(partition, district_path);

    // Get LCT path from u to v
    partition.lct.evert(u);
    partition.lct.expose(v);
    std::vector<int> path = partition.lct.find_path(v);

    if (path.size() < 2) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    // Compute collapsed populations and fragment map
    std::map<int, int> fragment_map;
    std::vector<int> collapsed_pop = compute_collapsed_pops(partition, path, fragment_map);

    // Setup cycle LCT
    std::vector<int> cycle_path = setup_cycle_lct(partition, path, collapsed_pop);

    // Guard: cycle_path must be sane
    if ((int)cycle_path.size() < 2 || (int)cycle_path.size() > partition.n_vertices) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    Rcpp::checkUserInterrupt();

    // Get balanced cuts
    BalancedCutsResult bcuts = get_balanced_cuts(partition, cycle_path, sub_dists, collapsed_pop);

    Rcpp::checkUserInterrupt();

    // Sample a cut
    if (bcuts.edge_inds.empty() || bcuts.cumWeight <= 0) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    double randSamp = r_unif() * bcuts.cumWeight;
    if (randSamp > bcuts.pathWeights[bcuts.pathWeights.size() - 2]) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    // Find selected edge
    int edge_ind = -1;
    for (int i = 0; i < (int)bcuts.pathWeights.size() - 1; i++) {
        if (randSamp > bcuts.pathWeights[i] && randSamp <= bcuts.pathWeights[i + 1]) {
            edge_ind = bcuts.edge_inds[i];
            break;
        }
    }
    if (edge_ind < 0) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    int n_cycle = (int)cycle_path.size();
    int cut_v1 = cycle_path[edge_ind];
    int cut_v2 = cycle_path[(edge_ind + 1) % n_cycle];
    int link_v1 = cycle_path[0];
    int link_v2 = cycle_path[n_cycle - 1];

    // Trivial case: edge_inds == sub_dists → just swap marked edge
    if ((int)bcuts.edge_inds.size() == sub_dists) {
        cut_district_path(partition, district_path);
        swap_marked_edge(partition, cut_v1, cut_v2, u, v);
        return {0.0, update}; // accepted with probability 1
    }

    // Non-trivial: compute forward/backward MH ratio

    // Random permutation for tie-breaking
    std::vector<int> node_perm(n_cycle);
    std::iota(node_perm.begin(), node_perm.end(), 0);
    for (int i = n_cycle - 1; i > 0; i--) {
        int j = r_int(i + 1);
        std::swap(node_perm[i], node_perm[j]);
    }

    // Find marked_edge_inds: which edge_inds (except last) are cross-district
    std::vector<int> marked_edge_inds;
    for (int i = 0; i < (int)bcuts.edge_inds.size() - 1; i++) {
        int ei = bcuts.edge_inds[i];
        int v1 = path[ei];      // Use original path for district lookup
        int v2 = path[ei + 1];
        if (partition.node_to_district[v1] != partition.node_to_district[v2]) {
            marked_edge_inds.push_back(ei);
        }
    }

    // Compute backward probability
    double backward_prob = get_path_cut_log_prob(
        partition, bcuts.edge_inds, marked_edge_inds, cycle_path, node_perm,
        collapsed_pop, lower, upper, sub_dists);

    // Guard: if backward probability is -inf, reject
    if (backward_prob < -1e17) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    Rcpp::checkUserInterrupt();

    // Cut and link on the cycle (cut selected edge, link endpoints)
    {
        // Restore cycle LCT to known state first
        partition.tmp_cycle_lct.evert(cycle_path[0]);
        set_cuttable_tree_info(partition.tmp_cycle_lct,
                               partition.tmp_cycle_lct.find_root(cycle_path[0]),
                               collapsed_pop, lower, upper, partition.cuttable_info);

        // Now cut the selected edge and link u-v on cycle
        partition.tmp_cycle_lct.evert(cut_v1);
        partition.tmp_cycle_lct.cut(cut_v2);
        // link_v1 should be root of its component
        safe_link(partition.tmp_cycle_lct, link_v1, link_v2, "bud_proposal_cycle_link");

        int new_root = partition.tmp_cycle_lct.find_root(cut_v2);
        set_cuttable_tree_info(partition.tmp_cycle_lct, new_root, collapsed_pop,
                               lower, upper, partition.cuttable_info);
    }

    // Rotate cycle_path, edge_inds, node_perm by edge_ind positions
    // New cycle starts at edge_ind+1 (the vertex after the cut).
    // The selected cut becomes the closing edge of the new cycle.
    std::vector<int> new_cycle_path;
    for (int i = edge_ind + 1; i < n_cycle; i++) new_cycle_path.push_back(cycle_path[i]);
    for (int i = 0; i <= edge_ind; i++) new_cycle_path.push_back(cycle_path[i]);

    std::vector<int> new_node_perm;
    for (int i = edge_ind + 1; i < n_cycle; i++) new_node_perm.push_back(node_perm[i]);
    for (int i = 0; i <= edge_ind; i++) new_node_perm.push_back(node_perm[i]);

    // Map old edge positions to new edge positions using vertex lookup
    // Old edge at position p: between cycle_path[p] and cycle_path[(p+1)%n]
    // Build new_cycle vertex→index map
    std::map<int, int> new_cp_idx;
    for (int i = 0; i < n_cycle; i++) {
        new_cp_idx[new_cycle_path[i]] = i;
    }

    // Map each old edge_ind to its new position in the rotated cycle
    std::vector<int> new_edge_inds;
    for (int i = 0; i < (int)bcuts.edge_inds.size(); i++) {
        int old_p = bcuts.edge_inds[i];
        if (old_p == edge_ind) {
            // Selected cut → closing edge (virtual position n_cycle)
            new_edge_inds.push_back(n_cycle);
        } else {
            // Find position of old edge's first vertex in new cycle
            int v1 = cycle_path[old_p];
            int new_p = new_cp_idx[v1];
            new_edge_inds.push_back(new_p);
        }
    }
    // Sort the new edge_inds (they should be in order, with n_cycle at the end)
    std::sort(new_edge_inds.begin(), new_edge_inds.end());

    // Forward sample: mark new edges on the modified cycle
    Rcpp::checkUserInterrupt();
    auto [new_marked_edges, forward_prob] = mark_edges_forward(
        partition, new_edge_inds, new_cycle_path, new_node_perm,
        collapsed_pop, lower, upper, sub_dists);

    // Guard: if forward marking failed, reject the proposal.
    if (forward_prob < -1e17 || (int)new_marked_edges.size() != sub_dists - 1) {
        cut_district_path(partition, district_path);
        return {0.0, update};
    }

    Rcpp::checkUserInterrupt();

    // Cut district path (restore to "rejected" state)
    cut_district_path(partition, district_path);

    // Build Update struct
    // new_marked_edges are indices in the rotated cycle_path
    // Add the closing edge
    new_marked_edges.push_back(n_cycle); // closing edge index (beyond end)

    // Build original cycle vertex→index map
    std::map<int, int> orig_cp_idx;
    for (int i = 0; i < n_cycle; i++) {
        orig_cp_idx[cycle_path[i]] = i;
    }

    // Map marked edges from rotated positions back to original cycle and build cuts
    update.valid = true;
    update.district_path = district_path;
    update.link_u = u;
    update.link_v = v;
    update.cycle_path = cycle_path;
    update.fragment_map = fragment_map;

    // Build cuts sorted by position in original cycle
    std::vector<std::pair<int, int>> order_entries; // (orig_pos, index_in_new_marked)
    for (int ni = 0; ni < (int)new_marked_edges.size(); ni++) {
        int nme = new_marked_edges[ni];
        int orig_pos;
        if (nme == n_cycle) {
            // Closing edge → selected cut at position edge_ind
            orig_pos = edge_ind;
        } else {
            // Normal edge: find original position via vertex
            int v1 = new_cycle_path[nme];
            orig_pos = orig_cp_idx[v1];
        }
        order_entries.push_back({orig_pos, ni});
    }
    std::sort(order_entries.begin(), order_entries.end());

    for (auto& [orig_pos, ni] : order_entries) {
        int cp_idx = orig_pos;
        int cp_next = (orig_pos + 1) % n_cycle;
        update.cuts.push_back({cycle_path[cp_idx], cycle_path[cp_next]});
        if (new_marked_edges[ni] == n_cycle) {
            // This is the actual cut edge
            update.cut_ind = (int)update.cuts.size() - 1;
        }
    }

    double p = std::exp(backward_prob - forward_prob);

    return {p, update};
}

// ============================================================
// 16a. compute_new_plan: compute new district assignments from update
//      without modifying the partition. Returns new plan as arma::uvec.
// ============================================================
arma::uvec compute_new_plan(const BUDPartition& partition, const BUDUpdate& update) {
    arma::uvec new_plan = partition.get_plan();

    // Collect affected districts (ordered) from district path
    std::vector<int> districts;
    for (auto& dp : update.district_path) {
        if (districts.empty() || districts.back() != dp.first) {
            districts.push_back(dp.first);
        }
    }
    districts.push_back(update.district_path.back().second);
    int sub_dists = (int)districts.size();

    int n_cycle = (int)update.cycle_path.size();
    std::map<int, int> vertex_to_cycle_pos;
    for (int i = 0; i < n_cycle; i++) {
        vertex_to_cycle_pos[update.cycle_path[i]] = i;
    }

    std::vector<int> cut_positions(sub_dists);
    for (int ii = 0; ii < sub_dists; ii++) {
        auto it = vertex_to_cycle_pos.find(update.cuts[ii].first);
        if (it != vertex_to_cycle_pos.end()) {
            cut_positions[ii] = it->second;
        }
    }

    // For each cycle vertex, determine its new district
    std::map<int, int> cycle_vertex_new_district;
    for (int ii = 0; ii < sub_dists; ii++) {
        int prev_cut_pos = cut_positions[(ii - 1 + sub_dists) % sub_dists];
        int this_cut_pos = cut_positions[ii];
        int pos = (prev_cut_pos + 1) % n_cycle;
        while (true) {
            cycle_vertex_new_district[update.cycle_path[pos]] = districts[ii];
            if (pos == this_cut_pos) break;
            pos = (pos + 1) % n_cycle;
        }
    }

    // Build reverse fragment_map
    std::map<int, std::vector<int>> fragments;
    for (auto& [graph_v, frag_root] : update.fragment_map) {
        fragments[frag_root].push_back(graph_v);
    }

    // Assign new districts
    for (auto& [cycle_v, new_d] : cycle_vertex_new_district) {
        auto fit = fragments.find(cycle_v);
        if (fit != fragments.end()) {
            for (int graph_v : fit->second) {
                new_plan(graph_v) = new_d;
            }
        }
    }

    return new_plan;
}

// ============================================================
// 16b. apply_bud_update: apply accepted update to partition
//      Uses incremental LCT updates (Julia-style: cut all first, then link)
//      to preserve the random spanning forest.
// ============================================================
void apply_bud_update(BUDPartition& partition, const BUDUpdate& update) {
    if (!update.valid) return;

    const Graph& g = *(partition.graph);
    const arma::uvec& pop = *(partition.pop);
    int V = partition.n_vertices;

    // Collect affected districts
    std::vector<int> districts;
    for (auto& dp : update.district_path) {
        if (districts.empty() || districts.back() != dp.first) {
            districts.push_back(dp.first);
        }
    }
    districts.push_back(update.district_path.back().second);
    std::set<int> district_set(districts.begin(), districts.end());
    int sub_dists = (int)districts.size();

    // ---- Step 1: Cut all edges in update.cuts from LCT ----
    std::vector<bool> actually_cut(update.cuts.size(), false);
    for (int ci = 0; ci < (int)update.cuts.size(); ci++) {
        auto [n1, n2] = update.cuts[ci];
        if (neighbors_in_lct(partition.lct, n1, n2)) {
            partition.lct.evert(n2);
            partition.lct.cut(n1);
            actually_cut[ci] = true;
        }
    }

    // ---- Step 2: Link old marked edges + (u,v), unless in cuts ----
    std::set<std::pair<int,int>> cuts_set;
    for (auto& c : update.cuts) {
        cuts_set.insert(c);
        cuts_set.insert({c.second, c.first});
    }

    std::vector<std::pair<int,int>> links;
    for (auto& dp : update.district_path) {
        DistrictPair key = {std::min(dp.first, dp.second), std::max(dp.first, dp.second)};
        auto it = partition.marked_edges.find(key);
        if (it != partition.marked_edges.end()) {
            links.push_back(it->second);
        }
    }
    links.push_back({update.link_u, update.link_v});

    std::vector<bool> actually_linked(links.size(), false);
    for (int li = 0; li < (int)links.size(); li++) {
        auto [l1, l2] = links[li];
        if (cuts_set.count({l1, l2}) || cuts_set.count({l2, l1})) continue;
        if (partition.lct.same_tree(l1, l2)) continue;
        partition.lct.evert(l1);
        partition.lct.link(l1, l2);
        actually_linked[li] = true;
    }

    // ---- Step 2b: Update forest_edges ----
    auto& fe = partition.forest_edges;
    // Remove cut edges
    for (int ci = 0; ci < (int)update.cuts.size(); ci++) {
        if (!actually_cut[ci]) continue;
        auto [n1, n2] = update.cuts[ci];
        fe.erase(std::remove_if(fe.begin(), fe.end(), [&](const std::pair<int,int>& e) {
            return (e.first == n1 && e.second == n2) ||
                   (e.first == n2 && e.second == n1);
        }), fe.end());
    }
    // Add linked edges
    for (int li = 0; li < (int)links.size(); li++) {
        if (!actually_linked[li]) continue;
        auto [l1, l2] = links[li];
        fe.push_back({l1, l2});
    }

    // ---- Step 3: Update district assignments ----
    int n_cycle = (int)update.cycle_path.size();
    std::map<int, int> vertex_to_cycle_pos;
    for (int i = 0; i < n_cycle; i++) {
        vertex_to_cycle_pos[update.cycle_path[i]] = i;
    }
    std::vector<int> cut_positions(sub_dists);
    for (int ii = 0; ii < sub_dists; ii++) {
        auto it = vertex_to_cycle_pos.find(update.cuts[ii].first);
        if (it != vertex_to_cycle_pos.end()) cut_positions[ii] = it->second;
    }
    std::map<int, int> cycle_vertex_new_district;
    for (int ii = 0; ii < sub_dists; ii++) {
        int prev_cut_pos = cut_positions[(ii - 1 + sub_dists) % sub_dists];
        int this_cut_pos = cut_positions[ii];
        int pos = (prev_cut_pos + 1) % n_cycle;
        while (true) {
            cycle_vertex_new_district[update.cycle_path[pos]] = districts[ii];
            if (pos == this_cut_pos) break;
            pos = (pos + 1) % n_cycle;
        }
    }
    std::map<int, std::vector<int>> fragments;
    for (auto& [graph_v, frag_root] : update.fragment_map) {
        fragments[frag_root].push_back(graph_v);
    }
    for (auto& [cycle_v, new_d] : cycle_vertex_new_district) {
        auto fit = fragments.find(cycle_v);
        if (fit != fragments.end()) {
            for (int graph_v : fit->second) {
                partition.node_to_district[graph_v] = new_d;
            }
        }
    }

    // ---- Step 4: Update district roots ----
    for (int ii = 0; ii < sub_dists; ii++) {
        int d = districts[ii];
        int root_v = update.cuts[ii].first;
        partition.district_roots[d] = partition.lct.find_root(root_v);
    }

    // ---- Step 5: Recompute district populations ----
    for (int d : districts) {
        partition.district_pop[d] = 0;
    }
    for (int v = 0; v < V; v++) {
        if (district_set.count(partition.node_to_district[v])) {
            partition.district_pop[partition.node_to_district[v]] += pop(v);
        }
    }

    // ---- Step 6: Update marked edges ----
    std::vector<DistrictPair> to_remove;
    std::map<DistrictPair, std::pair<int,int>> to_add;

    for (auto& [dp, edge] : partition.marked_edges) {
        int del_dists = (int)district_set.count(dp.first) +
                        (int)district_set.count(dp.second);
        if (del_dists == 0) continue;
        if (del_dists == 1) {
            int d1 = partition.node_to_district[edge.first];
            int d2 = partition.node_to_district[edge.second];
            DistrictPair new_dp = {std::min(d1, d2), std::max(d1, d2)};
            if (new_dp != dp) {
                to_remove.push_back(dp);
                to_add[new_dp] = edge;
            }
            continue;
        }
        to_remove.push_back(dp);
    }
    for (auto& dp : to_remove) partition.marked_edges.erase(dp);
    for (auto& [dp, edge] : to_add) partition.marked_edges[dp] = edge;

    for (int ci = 0; ci < (int)update.cuts.size(); ci++) {
        if (ci == update.cut_ind) continue;
        int n1 = update.cuts[ci].first;
        int n2 = update.cuts[ci].second;
        int d1 = partition.node_to_district[n1];
        int d2 = partition.node_to_district[n2];
        if (d1 != d2) {
            DistrictPair dp = {std::min(d1, d2), std::max(d1, d2)};
            partition.marked_edges[dp] = {n1, n2};
        }
    }

    // ---- Step 7: Rebuild district tree ----
    partition.district_tree.clear();
    partition.district_tree.resize(partition.n_districts);
    for (auto& [dp, edge] : partition.marked_edges) {
        partition.district_tree[dp.first].insert(dp.second);
        partition.district_tree[dp.second].insert(dp.first);
    }

    // ---- Step 8: Rebuild cross-district edges ----
    partition.update_cross_district_edges(g, districts);
}

// ============================================================
// 16c. revert_bud_update: undo apply_bud_update on rejection
//      Reverse the LCT cuts/links and restore partition state.
// ============================================================
static void revert_bud_update(BUDPartition& partition, const BUDUpdate& update,
                               const std::vector<std::pair<int,int>>& links,
                               const std::vector<bool>& actually_cut,
                               const std::vector<bool>& actually_linked,
                               const std::vector<int>& old_node_to_district,
                               const std::vector<int>& old_district_pop,
                               const std::vector<int>& old_district_roots,
                               const std::map<DistrictPair, std::pair<int,int>>& old_marked_edges,
                               const std::vector<std::set<int>>& old_district_tree,
                               const CrossEdgeMap& old_cross_edges) {
    // Undo LCT links (reverse order)
    for (int li = (int)links.size() - 1; li >= 0; li--) {
        if (!actually_linked[li]) continue;
        auto [l1, l2] = links[li];
        partition.lct.evert(l2);
        partition.lct.cut(l1);
    }
    // Undo LCT cuts (reverse order)  
    for (int ci = (int)update.cuts.size() - 1; ci >= 0; ci--) {
        if (!actually_cut[ci]) continue;
        auto [n1, n2] = update.cuts[ci];
        partition.lct.evert(n1);
        partition.lct.link(n1, n2);
    }

    // Restore partition metadata
    partition.node_to_district = old_node_to_district;
    partition.district_pop = old_district_pop;
    partition.district_roots = old_district_roots;
    partition.marked_edges = old_marked_edges;
    partition.district_tree = old_district_tree;
    partition.cross_edges = old_cross_edges;

    // Restore LCT roots
    for (int d = 0; d < partition.n_districts; d++) {
        partition.lct.evert(old_district_roots[d]);
    }
}

// ============================================================
// 16d. compute_log_wst: log weighted spanning tree count of district graph
//      Uses matrix-tree theorem on the small k×k district graph.
//      Edge weight = number of boundary edges between districts.
//      WST = det(L*) where L* is any (k-1)×(k-1) cofactor of weighted Laplacian.
// ============================================================
static double compute_log_wst(const CrossEdgeMap& cross_edges, int n_districts) {
    if (n_districts <= 1) return 0.0;

    // Build weighted Laplacian (k × k)
    int k = n_districts;
    arma::mat L(k, k, arma::fill::zeros);

    for (auto& [dp, edges] : cross_edges) {
        double w = (double)edges.size();
        if (w <= 0) continue;
        int i = dp.first;
        int j = dp.second;
        L(i, j) -= w;
        L(j, i) -= w;
        L(i, i) += w;
        L(j, j) += w;
    }

    // Remove last row and column to get (k-1) × (k-1) cofactor
    arma::mat L_star = L.submat(0, 0, k - 2, k - 2);

    // Compute log determinant
    double val, sign;
    arma::log_det(val, sign, L_star);

    if (sign <= 0) return -std::numeric_limits<double>::infinity();
    return val;
}

// ============================================================
// 17. bud_step: full BUD step with MH acceptance and energy
//     Forest-preserving approach:
//     - Rebuild LCT from forest_edges at start to ensure clean state
//     - Run proposal (may temporarily corrupt LCT internals)
//     - For p>0: rebuild again, apply incremental update, MH decision
//     - For p==0: forest walk already updated forest_edges
// ============================================================
int bud_step(BUDPartition& partition,
             double lower, double upper,
             double target,
             double compactness,
             const arma::uvec& counties,
             Rcpp::List constraints,
             double& accept_prob_out) {

    accept_prob_out = 0.0;

    // Rebuild LCT from forest_edges to ensure clean state for this step
    partition.rebuild_lct_from_forest();

    // Save old plan BEFORE proposal
    arma::uvec old_plan = partition.get_plan();

    // Save old configuration for potential revert
    auto old_forest_edges = partition.forest_edges;
    std::vector<int> old_node_to_district = partition.node_to_district;
    std::vector<int> old_district_pop = partition.district_pop;
    std::vector<int> old_district_roots = partition.district_roots;
    std::map<DistrictPair, std::pair<int,int>> old_marked_edges = partition.marked_edges;
    std::vector<std::set<int>> old_district_tree = partition.district_tree;
    CrossEdgeMap old_cross_edges = partition.cross_edges;

    // Run proposal
    double p;
    BUDUpdate update;
    try {
        auto [p_val, update_val] = bud_proposal(partition, lower, upper);
        p = p_val;
        update = std::move(update_val);
    } catch (const std::exception& e) {
        // LCT corrupted — restore old state
        partition.node_to_district = old_node_to_district;
        partition.district_pop = old_district_pop;
        partition.district_roots = old_district_roots;
        partition.marked_edges = old_marked_edges;
        partition.district_tree = old_district_tree;
        partition.cross_edges = old_cross_edges;
        partition.forest_edges = old_forest_edges;
        accept_prob_out = 0.0;
        return 0;
    }

    // If p == 0, changes already applied (forest walk or trivial)
    if (p == 0.0) {
        // Forest walk: internal_forest_walk_edge already updated forest_edges
        // Trivial swap / closing edge: forest_edges unchanged
        accept_prob_out = 1.0;
        return 0;
    }

    // Multi-district proposal (p > 0)
    // Rebuild from saved forest edges to get a clean pre-proposal LCT
    partition.node_to_district = old_node_to_district;
    partition.district_pop = old_district_pop;
    partition.district_roots = old_district_roots;
    partition.marked_edges = old_marked_edges;
    partition.district_tree = old_district_tree;
    partition.cross_edges = old_cross_edges;
    partition.forest_edges = old_forest_edges;
    partition.rebuild_lct_from_forest();

    // Apply the incremental update to the clean LCT
    apply_bud_update(partition, update);

    // Check population bounds
    bool pop_valid = true;
    for (int d = 0; d < partition.n_districts; d++) {
        if (partition.district_pop[d] < lower || partition.district_pop[d] > upper) {
            pop_valid = false;
            break;
        }
    }

    if (!pop_valid) {
        // Revert: restore old state
        partition.node_to_district = old_node_to_district;
        partition.district_pop = old_district_pop;
        partition.district_roots = old_district_roots;
        partition.marked_edges = old_marked_edges;
        partition.district_tree = old_district_tree;
        partition.cross_edges = old_cross_edges;
        partition.forest_edges = old_forest_edges;
        accept_prob_out = 0;
        return 0;
    }

    // Collect changed districts
    std::set<int> district_set;
    for (auto& dp : update.district_path) {
        district_set.insert(dp.first);
        district_set.insert(dp.second);
    }
    std::vector<int> changed_districts;
    for (int d : district_set) changed_districts.push_back(d + 1); // 1-indexed

    // Compute constraint penalty for OLD state
    Rcpp::NumericVector psi_vec;
    Rcpp::CharacterVector constr_names;
    if (constraints.size() > 0) {
        constr_names = constraints.names();
        psi_vec = Rcpp::NumericVector(constr_names.size());
        psi_vec.names() = constr_names;
    }

    double old_constraint = 0.0;
    if (constraints.size() > 0) {
        arma::subview_col<arma::uword> old_plan_view = old_plan.subvec(0, partition.n_vertices - 1);
        old_constraint = calc_gibbs_tgt(
            old_plan_view, partition.n_districts, partition.n_vertices,
            changed_districts, psi_vec, *(partition.pop), target,
            *(partition.graph), constraints);
    }

    // Compute constraint penalty for NEW state
    arma::uvec new_plan = partition.get_plan();

    double new_constraint = 0.0;
    if (constraints.size() > 0) {
        psi_vec = Rcpp::NumericVector(constr_names.size());
        psi_vec.names() = constr_names;
        arma::subview_col<arma::uword> new_plan_view = new_plan.subvec(0, partition.n_vertices - 1);
        new_constraint = calc_gibbs_tgt(
            new_plan_view, partition.n_districts, partition.n_vertices,
            changed_districts, psi_vec, *(partition.pop), target,
            *(partition.graph), constraints);
    }

    // Compute log spanning tree ratio (compactness)
    // BUD with preserved forest naturally targets ∝ |ST(x)| (same as CycleWalk)
    double log_st_ratio = 0.0;
    if (compactness != 1.0) {
        double log_st = 0.0;
        int n_cty = arma::max(counties);
        const Graph& g = *(partition.graph);

        arma::umat plans_mat(partition.n_vertices, 2);
        plans_mat.col(0) = old_plan;
        plans_mat.col(1) = new_plan;

        for (int d : district_set) {
            int distr = d + 1; // 1-indexed
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, plans_mat, counties, 0, distr, j);
                log_st -= log_st_distr(g, plans_mat, counties, 1, distr, j);
            }
            log_st += log_st_contr(g, plans_mat, counties, n_cty, 0, distr);
            log_st -= log_st_contr(g, plans_mat, counties, n_cty, 1, distr);
        }

        log_st_ratio = (1.0 - compactness) * log_st;
    }

    // BUD WST correction: BUD naturally targets ∝ |ST(x)| × WST(D_x)
    // To target ∝ |ST(x)|^compactness, we correct by WST_old/WST_new
    double log_wst_old = compute_log_wst(old_cross_edges, partition.n_districts);
    double log_wst_new = compute_log_wst(partition.cross_edges, partition.n_districts);

    // MH acceptance: p * wst_correction * st_correction * energy_ratio
    double log_mh = std::log(p);
    log_mh += log_wst_old - log_wst_new;
    log_mh += log_st_ratio;
    log_mh += old_constraint - new_constraint;

    double accept_prob = std::min(1.0, std::exp(log_mh));
    accept_prob_out = accept_prob;

    if (r_unif() < accept_prob) {
        // Accept — forest_edges already updated by apply_bud_update
        return 1;
    } else {
        // Reject — restore old state
        partition.node_to_district = old_node_to_district;
        partition.district_pop = old_district_pop;
        partition.district_roots = old_district_roots;
        partition.marked_edges = old_marked_edges;
        partition.district_tree = old_district_tree;
        partition.cross_edges = old_cross_edges;
        partition.forest_edges = old_forest_edges;
        return 0;
    }
}
