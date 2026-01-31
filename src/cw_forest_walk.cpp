/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * Internal Forest Walk Implementation
 ********************************************************/

#include "cw_forest_walk.h"
#include <algorithm>

bool get_random_internal_edge(LCTPartition& partition,
                               int& u, int& v,
                               int max_attempts) {
    const Graph& g = *(partition.graph);
    int V = partition.n_vertices;
    
    // Build list of all edges (as pairs where u < v)
    // This is O(E) but only done once per call
    std::vector<std::pair<int, int>> all_edges;
    for (int i = 0; i < V; i++) {
        for (int j : g[i]) {
            if (j > i) {  // Only count each edge once
                all_edges.push_back({i, j});
            }
        }
    }
    
    if (all_edges.empty()) return false;

    // Try to find an internal edge (both endpoints in same district)
    // Pick uniformly at random from all edges, like Julia does
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        // Pick a random edge uniformly
        int edge_idx = r_int((int)all_edges.size());
        u = all_edges[edge_idx].first;
        v = all_edges[edge_idx].second;

        // Check if same district (internal edge) by checking same root
        int root_u = partition.lct.find_root(u);
        int root_v = partition.lct.find_root(v);
        
        if (root_u == root_v) {
            return true;
        }
    }
    return false;
}

int internal_forest_walk(LCTPartition& partition, int max_attempts) {
    // Step 1: Pick a random internal edge
    int u, v;
    if (!get_random_internal_edge(partition, u, v, max_attempts)) {
        return 1;  // No internal edge found
    }

    LinkCutTree& lct = partition.lct;
    LCTNode* node_u = lct.node(u);
    LCTNode* node_v = lct.node(v);

    // Get the original root before we modify the tree
    int original_root = lct.find_root(u);

    // Step 2: Evert u to make it the root, then find path to v
    lct.evert(u);
    std::vector<int> path = lct.find_path(v);

    // Path goes from root (u) to v
    // If path length is 2, the edge (u,v) is already in the tree
    if (path.size() <= 2) {
        // The edge is already in the spanning tree, nothing to do
        // But we need to repair the partition (root may have changed)
        int new_root = lct.find_root(v);
        int district = partition.node_to_district[v];
        partition.district_roots[district] = new_root;
        return 0;
    }

    // Step 3: Compute cumulative weights for edges in the path
    // For unweighted graphs, each edge has weight 1, so 1/weight = 1
    // pathWeights[i] = cumulative weight up to edge (path[i-1], path[i])
    int path_len = (int)path.size();
    std::vector<double> cumulative_weights(path_len);
    cumulative_weights[0] = 0.0;

    double cum_weight = 0.0;
    for (int i = 1; i < path_len; i++) {
        // Weight of edge (path[i-1], path[i])
        double edge_weight = partition.get_edge_weight(path[i - 1], path[i]);
        cum_weight += 1.0 / edge_weight;
        cumulative_weights[i] = cum_weight;
    }

    // Add weight for the proposed new edge (u, v)
    double new_edge_weight = partition.get_edge_weight(u, v);
    double total_weight = cum_weight + 1.0 / new_edge_weight;

    // Step 4: Sample a random position
    double rand_sample = r_unif() * total_weight;

    // If we land past the last path edge, we're selecting the new edge
    // In that case, do nothing (the new edge is already "virtually" there)
    if (rand_sample > cumulative_weights[path_len - 1]) {
        // Selected the new edge - nothing to cut, just repair partition
        int new_root = lct.find_root(v);
        int district = partition.node_to_district[v];
        partition.district_roots[district] = new_root;
        return 0;
    }

    // Step 5: Find which edge to cut
    // We cut edge (path[edge_idx], path[edge_idx+1]) where
    // cumulative_weights[edge_idx] < rand_sample <= cumulative_weights[edge_idx+1]
    int edge_idx = -1;
    for (int i = 0; i < path_len - 1; i++) {
        if (rand_sample > cumulative_weights[i] &&
            rand_sample <= cumulative_weights[i + 1]) {
            edge_idx = i;
            break;
        }
    }

    if (edge_idx < 0) {
        // Should not happen
        return 1;
    }

    // Step 6: Cut the selected edge and link the new edge
    // We cut edge (path[edge_idx], path[edge_idx+1])
    // After evert(u), path[0] = u is the root
    // So path[edge_idx+1] is the node to cut (it's a child of path[edge_idx])

    int cut_child = path[edge_idx + 1];
    lct.cut(cut_child);

    // Find the new root BEFORE linking (like Julia does)
    // After cutting, v is in a different component from u
    // We need v's root to update the partition
    int new_root = lct.find_root(v);

    // Link u to v (u is already the root after evert)
    lct.link(u, v);

    // Step 7: Update partition - the district root may have changed
    // Use the new_root we found before linking
    int district = partition.node_to_district[v];
    partition.district_roots[district] = new_root;

    return 0;
}
