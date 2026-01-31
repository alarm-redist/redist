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

    // Path in district d1: from u2 to u1 (root)
    // First evert u1 to make it the root
    lct.evert(u1);
    path1 = lct.find_path(u2);
    // find_path returns root-to-u, but Julia wants u-to-root, so reverse
    std::reverse(path1.begin(), path1.end());

    // Path in district d2: from v2 to v1 (root)
    lct.evert(v1);
    path2 = lct.find_path(v2);
    // Same: reverse to get v2-to-root ordering like Julia
    std::reverse(path2.begin(), path2.end());

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

    // After reversal in get_cycle_paths:
    // path1 = [u2, ..., u1] where u1 is the root (at the END after reversal)
    // path2 = [v2, ..., v1] where v1 is the root (at the END after reversal)
    // Julia uses uPath[1] = u1 (BEFORE reversal) = path1.back() (AFTER reversal)
    
    // Get subtree populations rooted at u1 (LAST vertex of path1 after reversal)
    int u1 = path1.back();
    std::map<int, int> u_cut_pop = compute_subtree_pops(partition, u1);

    // Get subtree populations rooted at v1 (LAST vertex of path2 after reversal)
    int v1 = path2.back();
    std::map<int, int> v_cut_pop = compute_subtree_pops(partition, v1);

    // Build collapsed weights for path1
    // After reversal: path1 = [u2, ..., u1] (same as Julia's uPath_rev)
    // Julia iterates: uPath_rev[1] → collapsed[1], uPath_rev[2] → collapsed[2], etc.
    // In C++ (0-indexed): path1[0] → collapsed[0], path1[1] → collapsed[1], etc.
    for (size_t ii = 0; ii < path1.size(); ii++) {
        int vertex = path1[ii];  // Direct mapping, same order as Julia's uPath_rev
        collapsed_weights[ii] = u_cut_pop[vertex];
        if (ii > 0) {
            int prev_vertex = path1[ii - 1];  // Previous vertex in path (closer to u2)
            collapsed_weights[ii] -= u_cut_pop[prev_vertex];
        }
    }

    // Build collapsed weights for path2
    // After reversal: path2 = [v2, ..., v1]
    // Julia fills from end: vPath[end] → collapsed[end], vPath[end-1] → collapsed[end-1], etc.
    // Julia's vPath = [v1, ..., v2], so vPath[end] = v2
    // Our path2 = [v2, ..., v1], so we need path2[0]=v2 → collapsed[end], etc.
    for (size_t ii = 0; ii < path2.size(); ii++) {
        // Fill from end of cycle backwards
        int pos = cycle_len - 1 - ii;
        
        // path2[ii] gives us vertices from v2 towards v1
        // We want collapsed[end] = v2, collapsed[end-1] = next towards v1, etc.
        int vertex = path2[ii];
        collapsed_weights[pos] = v_cut_pop[vertex];
        
        if (ii > 0) {
            // Previous vertex in path2 is closer to v2
            int prev_vertex = path2[ii - 1];
            collapsed_weights[pos] -= v_cut_pop[prev_vertex];
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
    // prefix[i] = sum of cycle_pops[0..i-1]
    std::vector<int> prefix(n + 1, 0);
    for (int i = 0; i < n; i++) {
        prefix[i + 1] = prefix[i] + cycle_pops[i];
    }

    // Match Julia's find_cuttable_edge_pairs exactly:
    // Julia uses 1-indexed: cut1 = 1:path_length, cut2 = cut1:path_length-1
    // For pair (cut1, cut2), pop1 = sum(cycle_weights[cut1:cut2])
    //
    // In 0-indexed C++, pair (cut1, cut2) represents:
    // pop1 = sum of cycle_pops[cut1-1..cut2-1] = prefix[cut2] - prefix[cut1-1]
    //
    // Note: boundary edges at positions 1 and initial_cut+1 in Julia's indexing
    // are NOT excluded - they're handled by cancellation in get_cuts_and_links.
    // Only the identity pair (1, initial_cut) is removed.

    for (int cut1 = 1; cut1 <= n; cut1++) {
        for (int cut2 = cut1; cut2 <= n - 1; cut2++) {
            // Skip the identity configuration (no change to districts)
            // Julia: delete!(possible_pairs, (1, initial_cut_index))
            if (cut1 == 1 && cut2 == initial_cut) continue;

            // pop1 = sum of cycle_pops at positions cut1-1 to cut2-1 (0-indexed)
            int pop1 = prefix[cut2] - prefix[cut1 - 1];
            int pop2 = total_pop - pop1;

            if (pop1 >= lower && pop1 <= upper &&
                pop2 >= lower && pop2 <= upper) {
                valid_pairs.push_back({cut1, cut2});
            }
        }
    }

    return valid_pairs;
}

/*
 * Find the position of a link vertex in the cycle paths.
 * Returns 1-indexed position matching Julia.
 *
 * Corresponds to Julia's get_link_path_ind()
 */
static int get_link_path_ind(
    int link_vertex,
    const std::vector<int>& path1,
    const std::vector<int>& path2
) {
    int path1_len = (int)path1.size();
    int path2_len = (int)path2.size();

    // Check if link_vertex is at path1 end
    if (path1[path1_len - 1] == link_vertex) {
        return 1;
    }
    // Check if link_vertex is at path1 start
    else if (path1[0] == link_vertex) {
        return path1_len;
    }
    // Check if link_vertex is at path2 end
    else if (path2[path2_len - 1] == link_vertex) {
        return path1_len + path2_len;
    }
    // Check if link_vertex is at path2 start
    else if (path2[0] == link_vertex) {
        return path1_len + 1;
    }

    // Should never get here
    throw std::runtime_error("Couldn't find link vertex in paths");
}

/*
 * Determine if we should swap root assignments after the swap.
 *
 * This implements Julia's swap_assignment_check() logic.
 * The function determines which district should get which new root
 * based on population overlaps.
 */
static bool swap_assignment_check(
    int path_ind,
    int cut1,
    int cut2,
    const std::vector<int>& path1,
    const std::vector<int>& path2,
    const std::vector<int>& cycle_weights
) {
    int path1_len = (int)path1.size();
    int cycle_len = (int)cycle_weights.size();

    // Compute overlap1: population in the interval [cut1, cut2]
    // that comes from path1 (the "u" district)
    int overlap1 = 0;
    int tot_pop = 0;
    for (int w : cycle_weights) tot_pop += w;

    // Case 1: If cut1 <= path1_len, add weights in [cut1, min(path1_len, cut2)]
    if (cut1 <= path1_len) {
        for (int i = cut1; i <= std::min(path1_len, cut2); i++) {
            overlap1 += cycle_weights[i - 1];  // Convert to 0-indexed
        }
    }
    // Case 2: If cut1 > path1_len+1, add weights from [path1_len+1, cut1-1]
    else if (cut1 > path1_len + 1) {
        for (int i = path1_len + 1; i < cut1; i++) {
            overlap1 += cycle_weights[i - 1];
        }
    }

    // Case 3: If cut2 < cycle_len, add weights from [max(path1_len+1, cut2+1), end]
    if (cut2 < cycle_len) {
        for (int i = std::max(path1_len + 1, cut2 + 1); i <= cycle_len; i++) {
            overlap1 += cycle_weights[i - 1];
        }
    }

    // Determine which path gets more of the interval
    bool uPathToInterval = (2 * overlap1 > tot_pop);

    // Check positions
    bool l11_in_interval = (cut1 <= path_ind && path_ind <= cut2);
    bool l11_in_uPath = (path_ind <= path1_len);

    // XOR logic from Julia: (l11_in_uPath XOR l11_in_interval) XOR !uPathToInterval
    // This is equivalent to: !((l11_in_uPath XOR l11_in_interval) XOR uPathToInterval)
    // Or: (l11_in_uPath != l11_in_interval) == uPathToInterval
    return (l11_in_uPath != l11_in_interval) == uPathToInterval;
}

void apply_update(LCTPartition& partition,
                  const CycleWalkUpdate& update) {
    if (!update.valid) return;

    LinkCutTree& lct = partition.lct;

    // Apply cuts (Julia does NOT evert before cutting)
    for (const auto& cut : update.cuts) {
        // Julia: cut!(partition.lct.nodes[cut[2]])
        lct.cut(cut.second);
    }

    // Apply links (Julia: evert!(link[1]); link!(link[1], link[2]))
    for (const auto& link : update.links) {
        lct.evert(link.first);
        lct.link(link.first, link.second);
    }

    // Update district assignments and roots
    int d1 = update.changed_districts.first;
    int d2 = update.changed_districts.second;

    // Find new roots for the changed districts
    if (!update.cuts.empty()) {
        int new_root1 = lct.find_root(update.cuts[0].first);
        int new_root2 = lct.find_root(update.cuts[0].second);

        // Determine correct assignment using swap_link11 logic
        // This matches Julia's complex XOR logic for root assignment
        if (!update.links.empty()) {
            int link11_vertex = update.links[0].first;
            int link11_dist_cur = partition.node_to_district[link11_vertex];
            int r11_new_root = lct.find_root(link11_vertex);

            // Determine which new_root corresponds to link11's new tree
            int r11_root_ind_new = (r11_new_root != new_root1) ? 2 : 1;
            // Determine which district link11 currently belongs to
            int r11_root_ind_cur = (link11_dist_cur != d1) ? 2 : 1;

            // XOR logic: (r11_root_ind_new == r11_root_ind_cur) XOR !swap_link11
            // This is equivalent to: (equal AND swap) OR (!equal AND !swap)
            bool should_swap = (r11_root_ind_new == r11_root_ind_cur) == update.swap_link11;

            if (should_swap) {
                std::swap(new_root1, new_root2);
            }
        }

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
               CycleWalkDiagnostics& diagnostics) {
    // Initialize diagnostics
    diagnostics = CycleWalkDiagnostics();

    // Step 1: Pick random adjacent districts
    int d1, d2;
    if (!get_random_adjacent_districts(partition, d1, d2)) {
        diagnostics.status = -1;
        return -1;  // No adjacent districts
    }

    // Step 2: Pick two random boundary edges
    CWEdge e1(0, 0), e2(0, 0);
    if (!get_random_edge_pair(partition, d1, d2, e1, e2)) {
        diagnostics.status = -2;
        return -2;  // Need at least 2 boundary edges
    }

    // Step 3: Get the cycle paths
    std::vector<int> path1, path2;
    if (!get_cycle_paths(partition, e1, e2, path1, path2)) {
        diagnostics.status = -3;
        return -3;  // Couldn't get cycle paths
    }

    // Step 4: Get collapsed cycle weights (subtree populations)
    std::vector<int> cycle_pops = get_collapsed_cycle_weights(partition, path1, path2);

    // Track cycle length
    diagnostics.cycle_length = (int)cycle_pops.size();

    // Total population of the two districts
    int total_pop = partition.district_pop[d1] + partition.district_pop[d2];

    // Verify cycle integrity
    int cycle_pop_sum = 0;
    for (int p : cycle_pops) cycle_pop_sum += p;
    if (cycle_pop_sum != total_pop) {
        Rcpp::Rcout << "[ERROR] Cycle pop sum mismatch: sum=" << cycle_pop_sum
                    << ", total=" << total_pop << "\n";
    }

    // Initial cut is at the boundary between path1 and path2
    int initial_cut = (int)path1.size();

    // Step 5: Find valid cut pairs
    std::vector<std::pair<int, int>> valid_pairs =
        find_valid_cut_pairs(cycle_pops, initial_cut, total_pop, lower, upper);

    // Track number of valid cuts
    diagnostics.n_valid_cuts = (int)valid_pairs.size();

    if (valid_pairs.empty()) {
        // No valid cuts found - restore roots and return
        diagnostics.status = -4;
        partition.lct.evert(partition.district_roots[d1]);
        partition.lct.evert(partition.district_roots[d2]);
        return -4;  // No valid cut pairs found
    }

    // Save number of valid pairs for MH ratio
    int n_valid_pairs_fwd = (int)valid_pairs.size();

    // Step 6: Sample a cut pair weighted by edge weights
    // Build cumulative distribution weighted by 1/(w1*w2)
    auto get_edge_at_position = [&](int edge_ind) -> std::pair<int, int> {
        int path1_len = (int)path1.size();
        int path2_len = (int)path2.size();

        if (edge_ind == 1) {
            return {path1[0], path2[0]};
        } else if (edge_ind <= path1_len) {
            return {path1[edge_ind - 1], path1[edge_ind - 2]};
        } else if (edge_ind == path1_len + 1) {
            return {path1[path1_len - 1], path2[path2_len - 1]};
        } else {
            int ind = edge_ind - path1_len - 1;
            return {path2[path2_len - ind], path2[path2_len - ind - 1]};
        }
    };

    // Compute cumulative edge weight products for sampling
    std::vector<double> cum_edge_weight_product(n_valid_pairs_fwd);
    double cumsum = 0.0;

    for (int i = 0; i < n_valid_pairs_fwd; i++) {
        auto [c1, c2] = valid_pairs[i];

        // Get edges for this cut pair
        auto [e1_u, e1_v] = get_edge_at_position(c1);
        auto [e2_u, e2_v] = get_edge_at_position(c2 + 1);

        // Get edge weights from partition
        double w1 = partition.get_edge_weight(e1_u, e1_v);
        double w2 = partition.get_edge_weight(e2_u, e2_v);

        // Add 1/(w1*w2) to cumulative sum
        cumsum += 1.0 / (w1 * w2);
        cum_edge_weight_product[i] = cumsum;
    }

    // Sample proportional to 1/(w1*w2)
    double rand_samp = r_unif() * cum_edge_weight_product[n_valid_pairs_fwd - 1];
    int sample_idx = 0;
    while (sample_idx < n_valid_pairs_fwd - 1 &&
           rand_samp > cum_edge_weight_product[sample_idx]) {
        sample_idx++;
    }

    auto [cut1, cut2] = valid_pairs[sample_idx];

    // Get the actual edges for selected cuts (for edge weight ratio)
    auto [selected_e1_u, selected_e1_v] = get_edge_at_position(cut1);
    auto [selected_e2_u, selected_e2_v] = get_edge_at_position(cut2 + 1);
    double w1_cuts = partition.get_edge_weight(selected_e1_u, selected_e1_v);
    double w2_cuts = partition.get_edge_weight(selected_e2_u, selected_e2_v);
    double w1w2_cuts_inv = 1.0 / (w1_cuts * w2_cuts);

    // Get edge weights for the boundary edges (links)
    double w1_links = partition.get_edge_weight(e1.u, e1.v);
    double w2_links = partition.get_edge_weight(e2.u, e2.v);
    double w1w2_links_inv = 1.0 / (w1_links * w2_links);

    double sum_edge_weight_products = cum_edge_weight_product[n_valid_pairs_fwd - 1];

    // Step 7: Compute MH acceptance ratio
    int old_boundary = (int)partition.get_cross_edges(d1, d2).size();

    // Store old state for potential revert
    std::vector<int> old_node_to_district = partition.node_to_district;
    std::vector<int> old_district_pop = partition.district_pop;
    std::vector<int> old_district_roots = partition.district_roots;
    CrossEdgeMap old_cross_edges = partition.cross_edges;

    // Calculate constraint penalty for OLD state (before update)
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

    // Get edges at selected cut positions
    // Julia: e1 = get_node_indices_from_paths(edge_inds[1], ...)
    //        e2 = get_node_indices_from_paths(edge_inds[2]+1, ...)
    auto [fe1_u, fe1_v] = get_edge_at_position(cut1);
    auto [fe2_u, fe2_v] = get_edge_at_position(cut2 + 1);

    // Julia's get_cuts_and_links - handle boundary edge cancellation
    // Initial boundary edges (links candidates)
    // IMPORTANT: Use (u1,v1) and (u2,v2) ordering to match Julia's behavior
    // where link[0][0] is the endpoint in district d1
    // Extract the d1 endpoints from the boundary edges
    int e1_d1, e1_d2, e2_d1, e2_d2;
    if (partition.get_district(e1.u) == d1) {
        e1_d1 = e1.u; e1_d2 = e1.v;
    } else {
        e1_d1 = e1.v; e1_d2 = e1.u;
    }
    if (partition.get_district(e2.u) == d1) {
        e2_d1 = e2.u; e2_d2 = e2.v;
    } else {
        e2_d1 = e2.v; e2_d2 = e2.u;
    }
    
    std::vector<std::pair<int, int>> links = {{e1_d1, e1_d2}, {e2_d1, e2_d2}};
    // Selected cut edges
    std::vector<std::pair<int, int>> cuts = {{fe1_u, fe1_v}, {fe2_u, fe2_v}};

    // Helper to check if two edges are the same (order-independent)
    auto edges_equal = [](std::pair<int, int> a, std::pair<int, int> b) {
        return (a.first == b.first && a.second == b.second) ||
               (a.first == b.second && a.second == b.first);
    };

    // If a boundary edge is in cuts, it cancels (remove from both)
    if (edges_equal(links[0], cuts[0]) || edges_equal(links[0], cuts[1])) {
        // Remove link[0] and the matching cut
        if (edges_equal(links[0], cuts[0])) {
            cuts.erase(cuts.begin());
        } else {
            cuts.erase(cuts.begin() + 1);
        }
        links.erase(links.begin());
    } else if (edges_equal(links[1], cuts[0]) || edges_equal(links[1], cuts[1])) {
        // Remove link[1] and the matching cut
        if (edges_equal(links[1], cuts[0])) {
            cuts.erase(cuts.begin());
        } else {
            cuts.erase(cuts.begin() + 1);
        }
        links.erase(links.begin() + 1);
    }

    update.cuts = cuts;
    update.links = links;

    // Compute swap_link11 to determine correct root assignment
    if (!links.empty()) {
        int link11_vertex = links[0].first;
        int path_ind_l11 = get_link_path_ind(link11_vertex, path1, path2);
        update.swap_link11 = swap_assignment_check(
            path_ind_l11, cut1, cut2, path1, path2, cycle_pops
        );
    } else {
        update.swap_link11 = false;
    }

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

    // Count adjacent district pairs before and after
    // We count keys in cross_edges that involve d1 or d2
    auto count_adj_dists_involving = [](const CrossEdgeMap& cross_edges, int d1, int d2) {
        int count = 0;
        for (const auto& [key, edges] : cross_edges) {
            if (key.first == d1 || key.first == d2 ||
                key.second == d1 || key.second == d2) {
                if (!edges.empty()) count++;
            }
        }
        return count;
    };

    int old_adj_dists_d1d2 = count_adj_dists_involving(old_cross_edges, d1, d2);
    int new_adj_dists_d1d2 = count_adj_dists_involving(partition.cross_edges, d1, d2);
    int delta_adj_dists = new_adj_dists_d1d2 - old_adj_dists_d1d2;

    // Total number of adjacent district pairs (for proposal probability)
    int old_adj_dists_total = 0;
    for (const auto& [key, edges] : old_cross_edges) {
        if (!edges.empty()) old_adj_dists_total++;
    }

    // Compute log-MH ratio
    double log_mh_ratio = 0.0;
    double log_adj_ratio = 0.0;
    double log_edge_ratio = 0.0;
    double log_weight_ratio = 0.0;
    double log_constraint_ratio = 0.0;

    // Adjacent district ratio (accounts for change in number of adjacent pairs)
    if (old_adj_dists_total > 0 && old_adj_dists_total + delta_adj_dists > 0) {
        log_adj_ratio = std::log((double)old_adj_dists_total) -
                        std::log((double)(old_adj_dists_total + delta_adj_dists));
        log_mh_ratio += log_adj_ratio;
    }

    // Boundary edge ratio (proposal ratio)
    if (new_boundary > 0 && old_boundary > 0) {
        log_edge_ratio = std::log((double)(old_boundary * (old_boundary - 1))) -
                         std::log((double)(new_boundary * (new_boundary - 1)));
        log_mh_ratio += log_edge_ratio;
    }

    // Edge weight ratio for weighted graphs
    double new_sum_weights = sum_edge_weight_products + w1w2_links_inv - w1w2_cuts_inv;
    if (sum_edge_weight_products > 0 && new_sum_weights > 0) {
        log_weight_ratio = std::log(sum_edge_weight_products) - std::log(new_sum_weights);
        log_mh_ratio += log_weight_ratio;
    }

    // Constraint ratio (target ratio)
    log_constraint_ratio = old_constraint_penalty - new_constraint_penalty;
    log_mh_ratio += log_constraint_ratio;

    // Convert to acceptance probability
    double accept_ratio = std::min(1.0, std::exp(log_mh_ratio));

    // Track acceptance probability in diagnostics
    diagnostics.accept_prob = accept_ratio;

    // MH accept/reject
    double rand_val = r_unif();
    bool do_accept = rand_val < accept_ratio;
    
    if (do_accept) {
        diagnostics.status = 1;  // Accepted
        return 1;
    } else {
        // Reject - revert the update
        diagnostics.status = 0;  // Rejected
        // First undo LCT changes: cut the links, relink the cuts
        LinkCutTree& lct = partition.lct;
        for (const auto& link : update.links) {
            // Julia: evert!(link[2]), cut!(link[1])
            lct.evert(link.second);
            lct.cut(link.first);
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

        // Revert LCT roots
        partition.lct.evert(old_district_roots[d1]);
        partition.lct.evert(old_district_roots[d2]);

        return 0;  // Rejected
    }
}
