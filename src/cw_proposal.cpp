/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Cycle Walk Proposal Implementation
 ********************************************************/

#include "cw_proposal.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <map>

// Forward declaration for recursive helper
static int topological_sort_helper(
    std::map<int, int>& cut_pop,
    LCTNode* node,
    LCTNode* source,
    const LCTPartition& partition,
    bool reversed,
    int mass);

/*
 * Compute subtree populations for all vertices reachable from root.
 * After evert(root), traverses the LCT and for each vertex v, computes
 * the total population of the subtree rooted at v.
 *
 * This is equivalent to the Julia topological_sort function.
 */
static std::map<int, int> compute_subtree_pops(
    LCTPartition& partition,
    int root_vertex
) {
    std::map<int, int> cut_pop;
    LinkCutTree& lct = partition.lct;
    const arma::uvec& pop = *(partition.pop);

    // Evert to make root_vertex the root
    lct.evert(root_vertex);
    LCTNode* root = lct.node(root_vertex);

    // Push any pending reversal flags
    // The root after evert should have reversed=true, need to handle this
    bool rev = root->reversed;
    int lc = rev ? 1 : 0;
    int rc = 1 - lc;

    // Traverse the tree
    int total = 0;

    // Process right child (in tree order)
    if (root->children[rc] != nullptr) {
        total += topological_sort_helper(cut_pop, root->children[rc], root,
                                          partition, rev, 0);
    }

    // Process path children
    for (LCTNode* child : root->path_children) {
        total += topological_sort_helper(cut_pop, child, root, partition, false, 0);
    }

    // Set root's population
    int root_pop = pop(root_vertex);
    cut_pop[root_vertex] = root_pop + total;

    return cut_pop;
}

/*
 * Recursive helper for topological_sort.
 * Returns the total population of the subtree rooted at node.
 */
static int topological_sort_helper(
    std::map<int, int>& cut_pop,
    LCTNode* node,
    LCTNode* source,
    const LCTPartition& partition,
    bool reversed,
    int mass
) {
    if (node == nullptr) return 0;

    const arma::uvec& pop = *(partition.pop);
    int remainder = 0;

    // Handle reversal flag
    reversed = reversed ^ node->reversed;
    int lc = reversed ? 1 : 0;
    int rc = 1 - lc;

    // Process right subtree first (these are "below" us in tree order)
    if (node->children[rc] != nullptr) {
        remainder += topological_sort_helper(cut_pop, node->children[rc], node,
                                              partition, reversed, mass);
    }

    // Process path children
    for (LCTNode* child : node->path_children) {
        remainder += topological_sort_helper(cut_pop, child, node, partition, false, 0);
    }

    // This node's value
    int node_pop = pop(node->vertex);
    cut_pop[node->vertex] = remainder + node_pop + mass;

    // Process left subtree (if it's not where we came from)
    if (node->children[lc] != nullptr && node->children[lc] != source) {
        remainder += topological_sort_helper(cut_pop, node->children[lc], node,
                                              partition, reversed, cut_pop[node->vertex]);
    }

    return remainder + node_pop;
}

bool get_random_adjacent_districts(const LCTPartition& partition,
                                    int& d1, int& d2) {
    auto pairs = partition.get_adjacent_district_pairs();
    if (pairs.empty()) return false;

    int idx = r_int((int)pairs.size());
    d1 = pairs[idx].first;
    d2 = pairs[idx].second;
    return true;
}

bool get_random_edge_pair(const LCTPartition& partition,
                          int d1, int d2,
                          CWEdge& e1, CWEdge& e2) {
    const EdgeSet& edges = partition.get_cross_edges(d1, d2);
    int n_edges = (int)edges.size();

    // Need at least 2 edges to form a cycle
    if (n_edges < 2) return false;

    // Sample two distinct edges
    int idx1 = r_int(n_edges);
    int idx2 = r_int(n_edges - 1);
    if (idx2 >= idx1) idx2++;

    auto it = edges.begin();
    std::advance(it, idx1);
    e1 = *it;

    it = edges.begin();
    std::advance(it, idx2);
    e2 = *it;

    return true;
}

bool get_cycle_paths(LCTPartition& partition,
                     const CWEdge& e1, const CWEdge& e2,
                     std::vector<int>& path1, std::vector<int>& path2) {
    // e1 connects (u1, v1), e2 connects (u2, v2)
    // We need to determine which endpoints are in which district
    int u1 = e1.u, v1 = e1.v;
    int u2 = e2.u, v2 = e2.v;

    int d1 = partition.get_district(u1);
    int d2 = partition.get_district(v1);

    // Ensure u1, u2 are in d1 and v1, v2 are in d2
    if (partition.get_district(u1) != d1) std::swap(u1, v1);
    if (partition.get_district(u2) != d1) std::swap(u2, v2);

    // Now u1, u2 are in district d1 (or same district)
    // and v1, v2 are in district d2

    LinkCutTree& lct = partition.lct;

    // Path in district d1: from u1 to u2
    // First evert u1 to make it the root
    lct.evert(u1);
    path1 = lct.find_path(u2);

    // Path in district d2: from v1 to v2
    lct.evert(v1);
    path2 = lct.find_path(v2);

    return !path1.empty() && !path2.empty();
}

/*
 * Compute collapsed cycle weights - the "marginal" population that would
 * be transferred to the other district if we cut at each position.
 *
 * This matches the Julia get_collapsed_cycle_weights function.
 */
std::vector<int> get_collapsed_cycle_weights(
    LCTPartition& partition,
    const std::vector<int>& path1,
    const std::vector<int>& path2
) {
    int cycle_len = (int)(path1.size() + path2.size());
    std::vector<int> collapsed_weights(cycle_len);

    // Get subtree populations rooted at u1 (first vertex of path1)
    int u1 = path1[0];
    std::map<int, int> u_cut_pop = compute_subtree_pops(partition, u1);

    // Get subtree populations rooted at v1 (first vertex of path2)
    int v1 = path2[0];
    std::map<int, int> v_cut_pop = compute_subtree_pops(partition, v1);

    // Build collapsed weights for path1 (reversed order in cycle)
    // The cycle goes: path1_reversed, then path2_reversed
    // So position 0 is path1[end], position path1.size()-1 is path1[0]
    for (size_t ii = 0; ii < path1.size(); ii++) {
        int vertex = path1[path1.size() - 1 - ii];  // Reversed
        collapsed_weights[ii] = u_cut_pop[vertex];
        if (ii > 0) {
            int next_vertex = path1[path1.size() - ii];  // Previous in reversed order
            collapsed_weights[ii] -= u_cut_pop[next_vertex];
        }
    }

    // Build collapsed weights for path2 (reversed order in cycle)
    for (size_t ii = 0; ii < path2.size(); ii++) {
        int pos = path1.size() + ii;
        int vertex = path2[path2.size() - 1 - ii];  // Reversed
        collapsed_weights[pos] = v_cut_pop[vertex];
        if (ii > 0) {
            int next_vertex = path2[path2.size() - ii];  // Previous in reversed order
            collapsed_weights[pos] -= v_cut_pop[next_vertex];
        }
    }

    return collapsed_weights;
}

std::vector<std::pair<int, int>> find_valid_cut_pairs(
    const std::vector<int>& cycle_pops,
    int initial_cut,
    int total_pop,
    double lower, double upper) {

    std::vector<std::pair<int, int>> valid_pairs;
    int n = (int)cycle_pops.size();

    // Compute prefix sums for efficient range queries
    std::vector<int> prefix(n + 1, 0);
    for (int i = 0; i < n; i++) {
        prefix[i + 1] = prefix[i] + cycle_pops[i];
    }

    // The cycle consists of:
    // - path1: tree edges in district d1 (positions 1 to initial_cut-1)
    // - boundary edge e2 at position initial_cut (NOT a tree edge!)
    // - path2: tree edges in district d2 (positions initial_cut+1 to n-1)
    // - boundary edge e1 at position n/0 (NOT a tree edge!)
    //
    // We can ONLY cut tree edges, not boundary edges.
    // Boundary edge e2 is at position initial_cut
    // Boundary edge e1 is at position n (which wraps to 0)

    // Try all pairs of cut positions
    // A cut at position i means we cut before vertex i
    // So cutting at positions (i, j) where i < j gives us:
    //   - One part: vertices [i, j)
    //   - Other part: vertices [0, i) + [j, n)
    for (int i = 1; i < n; i++) {
        // Skip boundary edge e2 at position initial_cut
        if (i == initial_cut) continue;

        for (int j = i + 1; j <= n; j++) {
            // Skip boundary edge e1 at position n
            if (j == n) continue;

            // Skip boundary edge e2 at position initial_cut
            if (j == initial_cut) continue;

            // Skip the initial configuration (no change to districts)
            // The initial config has cuts at positions 1 and initial_cut+1
            // which corresponds to keeping the original district boundaries
            if (i == 1 && j == initial_cut + 1) continue;

            int pop1 = prefix[j] - prefix[i];
            int pop2 = total_pop - pop1;

            if (pop1 >= lower && pop1 <= upper &&
                pop2 >= lower && pop2 <= upper) {
                valid_pairs.push_back({i, j});
            }
        }
    }

    return valid_pairs;
}

void apply_update(LCTPartition& partition,
                  const CycleWalkUpdate& update) {
    if (!update.valid) return;

    LinkCutTree& lct = partition.lct;

    // Apply cuts
    for (const auto& cut : update.cuts) {
        lct.evert(cut.first);
        lct.cut(cut.second);
    }

    // Apply links
    for (const auto& link : update.links) {
        lct.evert(link.first);
        lct.link(link.first, link.second);
    }

    // Update district assignments and roots
    int d1 = update.changed_districts.first;
    int d2 = update.changed_districts.second;

    // Find new roots for the changed districts
    // The new roots are the roots of the trees after cuts and links
    if (!update.cuts.empty()) {
        int new_root1 = lct.find_root(update.cuts[0].first);
        int new_root2 = lct.find_root(update.cuts[0].second);
        partition.district_roots[d1] = new_root1;
        partition.district_roots[d2] = new_root2;
    }

    // Reassign districts by traversing from new roots
    // This is a simplified version - full implementation would use BFS
    for (int v = 0; v < partition.n_vertices; v++) {
        int root = lct.find_root(v);
        if (root == partition.district_roots[d1]) {
            partition.node_to_district[v] = d1;
        } else if (root == partition.district_roots[d2]) {
            partition.node_to_district[v] = d2;
        }
    }

    // Recompute populations
    partition.district_pop[d1] = 0;
    partition.district_pop[d2] = 0;
    const arma::uvec& pop = *(partition.pop);
    for (int v = 0; v < partition.n_vertices; v++) {
        if (partition.node_to_district[v] == d1) {
            partition.district_pop[d1] += pop(v);
        } else if (partition.node_to_district[v] == d2) {
            partition.district_pop[d2] += pop(v);
        }
    }

    // Recompute cross-district edges
    partition.cross_edges.clear();
    const Graph& g = *(partition.graph);
    for (int u = 0; u < partition.n_vertices; u++) {
        int d_u = partition.node_to_district[u];
        for (int v : g[u]) {
            if (v > u) {
                int d_v = partition.node_to_district[v];
                if (d_u != d_v) {
                    DistrictPair key(std::min(d_u, d_v), std::max(d_u, d_v));
                    partition.cross_edges[key].insert(CWEdge(u, v));
                }
            }
        }
    }
}

int cycle_walk(LCTPartition& partition,
               double lower, double upper,
               double target,
               Rcpp::List constraints,
               double& accept_ratio) {
    accept_ratio = 0.0;

    // Step 1: Pick random adjacent districts
    int d1, d2;
    if (!get_random_adjacent_districts(partition, d1, d2)) {
        return -1;  // No adjacent districts
    }

    // Step 2: Pick two random boundary edges
    CWEdge e1(0, 0), e2(0, 0);
    if (!get_random_edge_pair(partition, d1, d2, e1, e2)) {
        return -2;  // Need at least 2 boundary edges
    }

    // Step 3: Get the cycle paths
    std::vector<int> path1, path2;
    if (!get_cycle_paths(partition, e1, e2, path1, path2)) {
        return -3;  // Couldn't get cycle paths
    }

    // Step 4: Get collapsed cycle weights (subtree populations)
    std::vector<int> cycle_pops = get_collapsed_cycle_weights(partition, path1, path2);

    // Total population of the two districts
    int total_pop = partition.district_pop[d1] + partition.district_pop[d2];

    // Initial cut is at the boundary between path1 and path2
    int initial_cut = (int)path1.size();

    // Step 5: Find valid cut pairs
    std::vector<std::pair<int, int>> valid_pairs =
        find_valid_cut_pairs(cycle_pops, initial_cut, total_pop, lower, upper);

    if (valid_pairs.empty()) {
        // No valid cuts found - restore roots and return
        partition.lct.evert(partition.district_roots[d1]);
        partition.lct.evert(partition.district_roots[d2]);
        return -4;  // No valid cut pairs found
    }

    // Step 6: Sample a cut pair (uniform for now, should be weighted)
    int sample_idx = r_int((int)valid_pairs.size());
    auto [cut1, cut2] = valid_pairs[sample_idx];

    // Step 7: Compute MH acceptance ratio
    // For now, use a simple ratio based on number of valid pairs
    // Full implementation would account for boundary edge counts
    int old_boundary = (int)partition.get_cross_edges(d1, d2).size();

    // Store old state for potential revert
    std::vector<int> old_node_to_district = partition.node_to_district;
    std::vector<int> old_district_pop = partition.district_pop;
    std::vector<int> old_district_roots = partition.district_roots;
    CrossEdgeMap old_cross_edges = partition.cross_edges;

    // Calculate constraint penalty for OLD state (before update)
    // Districts that changed: d1 and d2
    std::vector<int> changed_districts = {d1, d2};
    Rcpp::NumericVector psi_vec;
    Rcpp::CharacterVector constr_names;
    if (constraints.size() > 0) {
        constr_names = constraints.names();
        psi_vec = Rcpp::NumericVector(constr_names.size());
        psi_vec.names() = constr_names;
    }

    // Get old plan as uvec for constraint evaluation
    arma::uvec old_plan(partition.n_vertices);
    for (int i = 0; i < partition.n_vertices; i++) {
        old_plan(i) = old_node_to_district[i] + 1;  // 1-indexed
    }
    arma::subview_col<arma::uword> old_plan_view = old_plan.subvec(0, partition.n_vertices - 1);

    double old_constraint_penalty = 0.0;
    if (constraints.size() > 0) {
        old_constraint_penalty = calc_gibbs_tgt(
            old_plan_view, partition.n_districts, partition.n_vertices,
            changed_districts, psi_vec, *(partition.pop), target,
            *(partition.graph), constraints
        );
    }

    // Create the update
    CycleWalkUpdate update;
    update.changed_districts = {d1, d2};
    update.valid = true;

    // Determine cuts and links based on selected cut positions
    // This is simplified - the full algorithm needs to track which edges
    // in the paths correspond to the cut positions

    // Get vertices at cut positions
    int cycle_len = (int)cycle_pops.size();
    auto get_cycle_vertex = [&](int pos) -> int {
        if (pos < (int)path1.size()) {
            return path1[pos];
        } else {
            return path2[path2.size() - 1 - (pos - path1.size())];
        }
    };

    // The cuts are at positions cut1 and cut2 in the cycle
    // We cut edges (cycle[cut1-1], cycle[cut1]) and (cycle[cut2-1], cycle[cut2 % cycle_len])
    int v_cut1_from = get_cycle_vertex(cut1 - 1);
    int v_cut1_to = get_cycle_vertex(cut1 % cycle_len);
    int v_cut2_from = get_cycle_vertex(cut2 - 1);
    int v_cut2_to = get_cycle_vertex(cut2 % cycle_len);

    update.cuts.push_back({v_cut1_from, v_cut1_to});
    update.cuts.push_back({v_cut2_from, v_cut2_to});

    // Links are the original boundary edges
    update.links.push_back({e1.u, e1.v});
    update.links.push_back({e2.u, e2.v});

    // Apply the update
    apply_update(partition, update);

    // Calculate constraint penalty for NEW state (after update)
    arma::uvec new_plan = partition.get_plan();
    arma::subview_col<arma::uword> new_plan_view = new_plan.subvec(0, partition.n_vertices - 1);

    double new_constraint_penalty = 0.0;
    if (constraints.size() > 0) {
        // Reset psi_vec for new calculation
        psi_vec = Rcpp::NumericVector(constr_names.size());
        psi_vec.names() = constr_names;
        new_constraint_penalty = calc_gibbs_tgt(
            new_plan_view, partition.n_districts, partition.n_vertices,
            changed_districts, psi_vec, *(partition.pop), target,
            *(partition.graph), constraints
        );
    }

    // Compute new boundary edge count
    int new_boundary = (int)partition.get_cross_edges(d1, d2).size();

    // Compute log-MH ratio
    double log_mh_ratio = 0.0;

    // Boundary edge ratio (proposal ratio)
    if (new_boundary > 0 && old_boundary > 0) {
        log_mh_ratio += std::log((double)(old_boundary * (old_boundary - 1))) -
                        std::log((double)(new_boundary * (new_boundary - 1)));
    }

    // Constraint ratio (target ratio)
    // Convention: subtract new, add old (because we want pi(new)/pi(old))
    // and constraints are log-penalties (higher = worse)
    log_mh_ratio += old_constraint_penalty - new_constraint_penalty;

    // Convert to acceptance probability
    accept_ratio = std::min(1.0, std::exp(log_mh_ratio));

    // MH accept/reject
    if (r_unif() < accept_ratio) {
        return 1;  // Accepted
    } else {
        // Reject - revert the update
        // First undo LCT changes: cut the links, relink the cuts
        LinkCutTree& lct = partition.lct;
        for (const auto& link : update.links) {
            lct.evert(link.first);
            lct.cut(link.second);
        }
        for (const auto& cut : update.cuts) {
            lct.evert(cut.first);
            lct.link(cut.first, cut.second);
        }

        // Restore partition state
        partition.node_to_district = old_node_to_district;
        partition.district_pop = old_district_pop;
        partition.district_roots = old_district_roots;
        partition.cross_edges = old_cross_edges;

        return 0;  // Rejected
    }
}
