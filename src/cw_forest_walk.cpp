#include "cw_forest_walk.h"

#include <algorithm>
#include <vector>

bool get_random_internal_edge(LCTGraphPlan &plan, MapParams const &map_params,
                              RNGState &rng_state, int &u, int &v, int max_attempts) {
    Graph const &g = map_params.g;
    int V = map_params.V;

    std::vector<std::pair<int, int>> all_edges;
    for (int i = 0; i < V; i++) {
        for (int j : g[i]) {
            if (j > i)
                all_edges.push_back({i, j});
        }
    }
    if (all_edges.empty())
        return false;

    LinkCutTree &lct = *plan.lct;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        int edge_idx = rng_state.r_int(static_cast<int>(all_edges.size()));
        u = all_edges[edge_idx].first;
        v = all_edges[edge_idx].second;
        if (lct.find_root(u) == lct.find_root(v))
            return true;
    }
    return false;
}

int internal_forest_walk(LCTGraphPlan &plan, MapParams const &map_params, RNGState &rng_state,
                         int max_attempts) {
    int u, v;
    if (!get_random_internal_edge(plan, map_params, rng_state, u, v, max_attempts))
        return 1;

    LinkCutTree &lct = *plan.lct;

    lct.evert(u);
    std::vector<int> path = lct.find_path(v);

    if (path.size() <= 2) {
        // Edge already in the tree.
        int new_root = lct.find_root(v);
        int region = plan.region_ids[v];
        plan.district_roots[region] = new_root;
        return 0;
    }

    int path_len = static_cast<int>(path.size());
    std::vector<double> cumulative_weights(path_len);
    cumulative_weights[0] = 0.0;
    double cum_weight = 0.0;
    for (int i = 1; i < path_len; i++) {
        double w = plan.get_edge_weight(path[i - 1], path[i]);
        cum_weight += 1.0 / w;
        cumulative_weights[i] = cum_weight;
    }
    double new_edge_weight = plan.get_edge_weight(u, v);
    double total_weight = cum_weight + 1.0 / new_edge_weight;

    double rand_sample = rng_state.r_unif() * total_weight;

    if (rand_sample > cumulative_weights[path_len - 1]) {
        int new_root = lct.find_root(v);
        int region = plan.region_ids[v];
        plan.district_roots[region] = new_root;
        return 0;
    }

    int edge_idx = -1;
    for (int i = 0; i < path_len - 1; i++) {
        if (rand_sample > cumulative_weights[i] && rand_sample <= cumulative_weights[i + 1]) {
            edge_idx = i;
            break;
        }
    }
    if (edge_idx < 0)
        return 1;

    int cut_child = path[edge_idx + 1];
    lct.cut(cut_child);
    int new_root = lct.find_root(v);
    lct.link(u, v);

    int region = plan.region_ids[v];
    plan.district_roots[region] = new_root;
    return 0;
}
