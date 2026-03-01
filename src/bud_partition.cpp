/********************************************************
 * BUD MCMC Redistricting Sampler
 * BUD Partition Implementation
 ********************************************************/

#include "bud_partition.h"
#include "wilson.h"
#include "random.h"
#include <algorithm>
#include <numeric>
#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

// ============================================================
// Cuttable tree closure operations
// ============================================================

NodeClosure init_node_closure(double val, double lower, double upper) {
    NodeClosure closure;
    closure[0] = {{val, val}};
    return closure;
}

std::vector<PopInterval> clean_closure(std::vector<PopInterval>& closure, double epsilon) {
    if (closure.empty()) return closure;
    std::sort(closure.begin(), closure.end());
    std::vector<PopInterval> cleaned;
    double t1 = closure[0].lo, t2 = closure[0].hi;
    for (size_t i = 1; i < closure.size(); i++) {
        double nt1 = closure[i].lo, nt2 = closure[i].hi;
        if (nt1 <= t2 + epsilon) {
            t2 = std::max(t2, nt2);
        } else {
            cleaned.push_back({t1, t2});
            t1 = nt1;
            t2 = nt2;
        }
    }
    cleaned.push_back({t1, t2});
    return cleaned;
}

void merge_closures(NodeClosure& c1, const NodeClosure& c2, double lower, double upper) {
    if (c2.empty()) return;
    if (c1.empty()) {
        c1 = c2;
        return;
    }

    NodeClosure merged;
    long total_ops = 0;
    for (auto& [d1, ivs1] : c1) {
        for (auto& [d2, ivs2] : c2) {
            int dnew = d1 + d2;
            auto& dest = merged[dnew];
            for (auto& ec1 : ivs1) {
                for (auto& ec2 : ivs2) {
                    dest.push_back({ec1.lo + ec2.lo, ec1.hi + ec2.hi});
                    if (++total_ops > 1000000) {
                        throw std::runtime_error("merge_closures: too many operations");
                    }
                }
            }
            std::sort(dest.begin(), dest.end());
            // Remove intervals entirely above upper bound
            while (!dest.empty() && dest.back().lo > upper) {
                dest.pop_back();
            }
        }
    }

    // Clean overlapping intervals
    double epsilon = upper - lower;
    for (auto it = merged.begin(); it != merged.end(); ) {
        if (it->second.empty()) {
            it = merged.erase(it);
        } else {
            it->second = clean_closure(it->second, epsilon);
            ++it;
        }
    }

    c1 = merged;
}

void add_cuts_to_closures(NodeClosure& closure, double lower, double upper) {
    std::vector<int> flagged;
    for (auto& [dists, ivs] : closure) {
        if (ivs.back().hi >= lower) {
            flagged.push_back(dists + 1);
        }
    }
    double epsilon = upper - lower;
    for (int d : flagged) {
        if (closure.find(d) == closure.end()) {
            closure[d] = {{0, 0}};
        } else {
            closure[d].insert(closure[d].begin(), {0, 0});
            closure[d] = clean_closure(closure[d], epsilon);
        }
    }
}

bool cuttable_huh(const NodeClosure& root_closure, int num_dists) {
    auto it = root_closure.find(num_dists);
    if (it == root_closure.end()) return false;
    return !it->second.empty() && it->second[0].lo == 0;
}

// ============================================================
// LCT traversal helpers for cuttable info computation
// ============================================================

/*
 * Recursive helper to compute cuttable info for an LCT subtree.
 * Returns the closure for this subtree.
 */
static int g_cuttable_call_count = 0;

static NodeClosure get_cuttable_tree_info_helper(
    LinkCutTree& lct,
    LCTNode* node,
    const std::vector<int>& pop_values,
    double lower, double upper,
    CuttableInfo& info,
    bool reversed,
    NodeClosure prev_closures
) {
    if (node == nullptr) return prev_closures;

    if (++g_cuttable_call_count % 1000 == 0) {
        Rcpp::checkUserInterrupt();
    }

    reversed = reversed ^ node->reversed;
    int lc_idx = reversed ? 1 : 0;
    int rc_idx = 1 - lc_idx;

    // Init closure for this node
    NodeClosure node_closure = init_node_closure(pop_values[node->vertex], lower, upper);

    // Process right child (deeper in preferred path)
    NodeClosure child_closure = get_cuttable_tree_info_helper(
        lct, node->children[rc_idx], pop_values, lower, upper, info, reversed, prev_closures);
    merge_closures(node_closure, child_closure, lower, upper);

    // Process path children
    for (LCTNode* pc : node->path_children) {
        NodeClosure pc_closure = get_cuttable_tree_info_helper(
            lct, pc, pop_values, lower, upper, info, false, NodeClosure());
        merge_closures(node_closure, pc_closure, lower, upper);
    }

    add_cuts_to_closures(node_closure, lower, upper);
    info.node_closures[node->vertex] = node_closure;

    // Continue along left child (shallower in preferred path)
    return get_cuttable_tree_info_helper(
        lct, node->children[lc_idx], pop_values, lower, upper, info, reversed, node_closure);
}

void set_cuttable_tree_info(LinkCutTree& lct, int root_vertex,
                            const std::vector<int>& pop_values,
                            double lower, double upper,
                            CuttableInfo& info) {
    lct.evert(root_vertex);
    LCTNode* root = lct.node(root_vertex);

    bool rev = root->reversed;
    int lc_idx = rev ? 1 : 0;
    int rc_idx = 1 - lc_idx;

    NodeClosure root_closure = init_node_closure(pop_values[root_vertex], lower, upper);

    NodeClosure child_closure = get_cuttable_tree_info_helper(
        lct, root->children[rc_idx], pop_values, lower, upper, info, rev, NodeClosure());
    merge_closures(root_closure, child_closure, lower, upper);

    for (LCTNode* pc : root->path_children) {
        NodeClosure pc_closure = get_cuttable_tree_info_helper(
            lct, pc, pop_values, lower, upper, info, false, NodeClosure());
        merge_closures(root_closure, pc_closure, lower, upper);
    }

    add_cuts_to_closures(root_closure, lower, upper);
    info.node_closures[root_vertex] = root_closure;
    info.root = root_vertex;
}

void set_cuttable_path_info(LinkCutTree& lct,
                            const std::vector<int>& path,
                            const std::vector<int>& pop_values,
                            double lower, double upper,
                            CuttableInfo& info) {
    // path[0] is deepest, path.back() is root
    // Process from leaves toward root
    for (size_t i = 0; i < path.size(); i++) {
        int v = path[i];
        LCTNode* node = lct.node(v);
        NodeClosure node_closure = init_node_closure(pop_values[v], lower, upper);

        // Get all LCT neighbors
        lct.expose(v);
        // Collect neighbor vertices
        std::vector<int> neighbors;
        // Left child in splay tree
        auto get_rightmost = [](LCTNode* n, bool rev) -> LCTNode* {
            if (!n) return nullptr;
            while (true) {
                rev ^= n->reversed;
                int rc = rev ? 0 : 1;
                if (n->children[rc] == nullptr) return n;
                n = n->children[rc];
            }
        };
        auto get_leftmost = [](LCTNode* n, bool rev) -> LCTNode* {
            if (!n) return nullptr;
            while (true) {
                rev ^= n->reversed;
                int lc = rev ? 1 : 0;
                if (n->children[lc] == nullptr) return n;
                n = n->children[lc];
            }
        };

        bool rev = node->reversed;
        int lc_idx = rev ? 1 : 0;
        int rc_idx = 1 - lc_idx;
        if (node->children[lc_idx]) {
            LCTNode* rm = get_rightmost(node->children[lc_idx], rev);
            if (rm) neighbors.push_back(rm->vertex);
        } else if (node->path_parent) {
            neighbors.push_back(node->path_parent->vertex);
        }
        if (node->children[rc_idx]) {
            LCTNode* lm = get_leftmost(node->children[rc_idx], rev);
            if (lm) neighbors.push_back(lm->vertex);
        }
        for (LCTNode* pc : node->path_children) {
            LCTNode* lm = get_leftmost(pc, false);
            if (lm) neighbors.push_back(lm->vertex);
        }

        // Determine parent (next node on path)
        int parent = -1;
        if (i + 1 < path.size()) parent = path[i + 1];

        // Merge closures from non-parent neighbors
        for (int nb : neighbors) {
            if (nb != parent) {
                merge_closures(node_closure, info.node_closures[nb], lower, upper);
            }
        }

        add_cuts_to_closures(node_closure, lower, upper);
        info.node_closures[v] = node_closure;
    }

    info.root = path.back();
}

// ============================================================
// BUDPartition implementation
// ============================================================

BUDPartition::BUDPartition(int n_vertices, int n_districts)
    : LCTPartition(n_vertices, n_districts),
      district_tree(n_districts),
      cuttable_info(n_vertices),
      tmp_cycle_lct(n_vertices),
      pop_values(n_vertices, 0),
      pop_lower(0), pop_upper(0) {
}

int BUDPartition::init_from_plan(const Graph& g,
                                  const arma::uvec& plan,
                                  const arma::uvec& population,
                                  const arma::uvec& county_assignments,
                                  double lower, double upper) {
    // Initialize base members (don't call parent init_from_plan,
    // we replicate its logic to capture forest edges)
    graph = &g;
    pop = &population;
    counties = &county_assignments;

    int V = n_vertices;

    lct = LinkCutTree(V);
    forest_edges.clear();
    tmp_cycle_lct = LinkCutTree(V);
    district_tree.assign(n_districts, std::set<int>());
    pop_lower = lower;
    pop_upper = upper;

    for (int i = 0; i < V; i++) {
        pop_values[i] = (int)population(i);
    }

    // For each district, sample spanning tree via Wilson's and record edges
    Multigraph cg = county_graph(g, county_assignments);

    for (int d = 0; d < n_districts; d++) {
        std::vector<bool> ignore(V);
        int n_in_district = 0;
        for (int i = 0; i < V; i++) {
            if ((int)plan(i) == d + 1) {
                ignore[i] = false;
                n_in_district++;
            } else {
                ignore[i] = true;
            }
        }
        if (n_in_district == 0) {
            Rcpp::stop("District %d has no vertices", d + 1);
        }

        Tree tree = init_tree(V);
        int root;
        std::vector<bool> visited(V);
        int result = sample_sub_ust(g, tree, V, root, visited, ignore,
                                     population, lower, upper, county_assignments, cg);
        if (result != 0) return 1;

        district_roots[d] = root;

        // Load tree into LCT AND record forest edges
        std::queue<int> queue;
        queue.push(root);
        node_to_district[root] = d;
        while (!queue.empty()) {
            int v = queue.front();
            queue.pop();
            for (int child : tree[v]) {
                lct.evert(child);
                lct.link(child, v);
                node_to_district[child] = d;
                forest_edges.push_back({child, v});
                queue.push(child);
            }
        }

        district_pop[d] = 0;
        for (int i = 0; i < V; i++) {
            if ((int)plan(i) == d + 1) {
                district_pop[d] += population(i);
            }
        }
    }

    find_cross_district_edges();

    // Sample district tree and marked edges
    sample_district_tree();
    sample_marked_edges();

    return 0;
}

void BUDPartition::sample_district_tree() {
    // Build weighted district graph
    int nd = n_districts;
    // district adjacency with weights (sum of edge weights between districts)
    std::vector<std::vector<std::pair<int, double>>> dist_adj(nd);

    for (auto& [dp, edges] : cross_edges) {
        if (edges.empty()) continue;
        double weight = 0;
        for (auto& e : edges) {
            weight += e.weight;
        }
        int d1 = dp.first;
        int d2 = dp.second;
        dist_adj[d1].push_back({d2, weight});
        dist_adj[d2].push_back({d1, weight});
    }

    // Wilson's algorithm on district graph
    std::vector<bool> in_tree(nd, false);
    std::vector<int> next(nd, -1);

    // Start from district 0
    in_tree[0] = true;

    for (int i = 1; i < nd; i++) {
        if (in_tree[i]) continue;

        // Random walk from i until hitting tree
        int u = i;
        while (!in_tree[u]) {
            if (dist_adj[u].empty()) {
                // Disconnected district graph - shouldn't happen
                return;
            }
            // Sample neighbor proportional to weight
            double total_w = 0;
            for (auto& [nb, w] : dist_adj[u]) {
                total_w += w;
            }
            double r = r_unif() * total_w;
            double cum = 0;
            int chosen = dist_adj[u][0].first;
            for (auto& [nb, w] : dist_adj[u]) {
                cum += w;
                if (r <= cum) {
                    chosen = nb;
                    break;
                }
            }
            next[u] = chosen;
            u = chosen;
        }

        // Trace back and add to tree
        u = i;
        while (!in_tree[u]) {
            in_tree[u] = true;
            int v = next[u];
            // Add edge u-v to district tree
            district_tree[u].insert(v);
            district_tree[v].insert(u);
            u = v;
        }
    }
}

void BUDPartition::sample_marked_edges() {
    marked_edges.clear();

    // For each edge in district tree, sample a node-level edge
    std::set<DistrictPair> visited;
    for (int d = 0; d < n_districts; d++) {
        for (int nb : district_tree[d]) {
            DistrictPair dp = {std::min(d, nb), std::max(d, nb)};
            if (visited.count(dp)) continue;
            visited.insert(dp);

            // Get cross-district edges for this pair
            auto it = cross_edges.find(dp);
            if (it == cross_edges.end() || it->second.empty()) continue;

            const EdgeSet& edges = it->second;
            // Sample proportional to weight
            double total_w = 0;
            for (auto& e : edges) {
                total_w += e.weight;
            }
            double r = r_unif() * total_w;
            double cum = 0;
            std::pair<int, int> sampled = {edges.begin()->u, edges.begin()->v};
            for (auto& e : edges) {
                cum += e.weight;
                if (r <= cum) {
                    sampled = {e.u, e.v};
                    break;
                }
            }
            marked_edges[dp] = sampled;
        }
    }
}

std::vector<DistrictPair> BUDPartition::get_district_path(int d1, int d2) const {
    // BFS on district tree
    std::vector<int> parent(n_districts, -1);
    std::vector<bool> visited(n_districts, false);
    std::queue<int> q;
    q.push(d1);
    visited[d1] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        if (u == d2) break;
        for (int nb : district_tree[u]) {
            if (!visited[nb]) {
                visited[nb] = true;
                parent[nb] = u;
                q.push(nb);
            }
        }
    }

    // Reconstruct path as district pairs
    std::vector<DistrictPair> path;
    if (!visited[d2]) return path; // d2 not reachable
    int cur = d2;
    while (cur != d1) {
        int p = parent[cur];
        if (p < 0) return std::vector<DistrictPair>(); // safety
        path.push_back({std::min(p, cur), std::max(p, cur)});
        cur = p;
    }
    std::reverse(path.begin(), path.end());
    return path;
}

void BUDPartition::rebuild_district_tree() {
    for (int d = 0; d < n_districts; d++) {
        district_tree[d].clear();
    }
    for (auto& [dp, edge] : marked_edges) {
        district_tree[dp.first].insert(dp.second);
        district_tree[dp.second].insert(dp.first);
    }
}

void BUDPartition::extract_forest_edges() {
    // Extract forest edges by traversing the tree from each district root.
    // Uses DFS through graph edges, checking LCT connectivity.
    forest_edges.clear();

    // For each district, do BFS from root.
    // Since within a district, the spanning tree is a tree on a connected subgraph,
    // we can identify tree edges by checking: after evert(root), for each node u,
    // find its tree-parent by the LCT expose path.
    // Simpler approach: BFS from root, for each unvisited same-district neighbor v of u,
    // try cutting v from the LCT - if u and v become disconnected, (u,v) was a tree edge.
    // But this is destructive.
    //
    // Most practical approach: BFS using the graph edges. In a spanning tree of the
    // district's induced subgraph, each node except root has exactly one parent edge.
    // We can find this by: for each node u, evert(root), then the LCT path root->u
    // gives us the tree path; the last edge is (parent_of_u, u).
    //
    // But this is O(V * log V). Let's use a different approach:
    // We know the total number of tree edges is V - n_districts.
    // BFS from root, for each unvisited same-district neighbor, check if it's
    // in the same LCT component as the root (always yes in a spanning tree).
    // So we must use a different criterion.
    //
    // Correct approach: evert(root), then for each node u (except root),
    // find its parent by exposing u and getting the previous node on the path.
    // This is complex with splay trees.
    //
    // Simplest correct approach: for each node u, check all same-district graph
    // neighbors v > u. For each, temporarily cut and check. O(E log V).
    // Too expensive for general use, so let's only call this when needed.

    std::vector<bool> visited(n_vertices, false);
    for (int d = 0; d < n_districts; d++) {
        int root = district_roots[d];
        if (root < 0 || visited[root]) continue;

        // Evert root to make it the actual root of the LCT tree for this district
        lct.evert(root);

        // For each node in this district, find its parent via the LCT path from node to root
        // The parent is the second node on the path from node to root after evert(root)
        std::vector<int> district_nodes;
        for (int i = 0; i < n_vertices; i++) {
            if (node_to_district[i] == d) {
                district_nodes.push_back(i);
                visited[i] = true;
            }
        }

        for (int u : district_nodes) {
            if (u == root) continue;
            // After evert(root), expose(u) gives us the path from root to u
            // The parent of u is the node just before u on this path
            lct.expose(u);
            // Get the path from root to u
            std::vector<int> path = lct.find_path(u);
            if (path.size() >= 2) {
                // Parent of u is path[path.size()-2]
                int parent = path[path.size() - 2];
                forest_edges.push_back({u, parent});
            }
        }
    }
}

void BUDPartition::rebuild_lct_from_forest() {
    lct = LinkCutTree(n_vertices);
    for (auto& [u, v] : forest_edges) {
        lct.evert(u);
        lct.link(u, v);
    }
}

void BUDPartition::update_cross_district_edges(const Graph& g,
                                                const std::vector<int>& affected_districts) {
    std::set<int> affected_set(affected_districts.begin(), affected_districts.end());

    // Remove all cross-edge entries involving affected districts
    std::vector<DistrictPair> to_remove;
    for (auto& [key, edges] : cross_edges) {
        if (affected_set.count(key.first) || affected_set.count(key.second)) {
            to_remove.push_back(key);
        }
    }
    for (auto& key : to_remove) {
        cross_edges.erase(key);
    }

    // Rescan edges for affected districts
    for (int u = 0; u < n_vertices; u++) {
        int d_u = node_to_district[u];
        if (!affected_set.count(d_u)) continue;
        for (int v : g[u]) {
            if (v > u) {
                int d_v = node_to_district[v];
                if (d_u != d_v) {
                    DistrictPair key(std::min(d_u, d_v), std::max(d_u, d_v));
                    cross_edges[key].insert(CWEdge(u, v));
                }
            }
        }
    }
    // Also scan edges where u is not affected but v is
    for (int u = 0; u < n_vertices; u++) {
        int d_u = node_to_district[u];
        if (affected_set.count(d_u)) continue;
        for (int v : g[u]) {
            if (v > u) {
                int d_v = node_to_district[v];
                if (d_v != d_u && affected_set.count(d_v)) {
                    DistrictPair key(std::min(d_u, d_v), std::max(d_u, d_v));
                    cross_edges[key].insert(CWEdge(u, v));
                }
            }
        }
    }
}

void BUDPartition::print_bud_state(int verbosity) const {
    LCTPartition::print_state(verbosity);
    if (verbosity >= 2) {
        Rcout << "\nDistrict tree edges: " << marked_edges.size() << "\n";
        for (auto& [dp, edge] : marked_edges) {
            Rcout << "  Districts (" << dp.first << "," << dp.second
                  << ") -> nodes (" << edge.first << "," << edge.second << ")\n";
        }
    }
}
