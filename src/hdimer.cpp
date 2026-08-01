#include "hdimer.h"

#include "hdimer_matching.h"
#include "hungarian.h"
#include "mcmc_gibbs.h"
#include "random.h"
#include "smc_base.h"
#include "tree_op.h"

#include <algorithm>
#include <chrono>
#include <cli/progress.h>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_map>
#include <utility>

using namespace arma;
using namespace Rcpp;

namespace {

constexpr double CUT_SUPPORT = 1.0;

using Clock = std::chrono::steady_clock;
using ForestAdjacency = std::vector<std::vector<std::pair<int, int>>>;

struct MultiEdge {
    int from;
    int to;
    int u;
    int v;
};

struct Fragment {
    std::vector<int> vertices;
    std::vector<std::vector<int>> units;
    int old_district = -1;
    double population = 0.0;
    double cut_normalizer = 0.0;
};

struct CandidateEdge {
    int u;
    int v;
    int fragment_u;
    int fragment_v;
    double log_weight;
};

struct BaseEdge {
    int u;
    int v;
    double cut_bias;
};

struct PairCache {
    int boundary_edges = 0;
    bool population_valid = false;
    bool hierarchy_valid = false;
    double log_target = 0.0;
    std::vector<int> shared_units;
};

struct HDimerWorkspace {
    HDimerWorkspace(int n_vertices, int n_districts)
        : full_adjacency(n_vertices), half_adjacency(n_vertices), parent(n_vertices),
          parent_edge(n_vertices), roots(n_districts), district_sizes(n_districts),
          subtree(n_vertices), district_order(n_districts), fragment_of(n_vertices),
          reroot_sum(n_vertices), district_vertices(n_districts), cut_edge_marks(n_vertices) {
        fragments.reserve(2 * n_districts);
    }

    ForestAdjacency full_adjacency;
    ForestAdjacency half_adjacency;
    std::vector<int> parent;
    std::vector<int> parent_edge;
    std::vector<int> roots;
    std::vector<int> district_sizes;
    std::vector<double> subtree;
    std::vector<std::vector<int>> district_order;
    std::vector<int> fragment_of;
    std::vector<Fragment> fragments;
    std::vector<double> reroot_sum;
    std::vector<std::vector<int>> district_vertices;
    std::vector<int> cut_edge_marks;
    int cut_edge_mark = 0;
};

double elapsed_seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

[[noreturn]] void stop_dimer_numerical_failure(const char *phase) {
    Rcpp::stop("Numerical failure while %s; no chain was returned.", phase);
}

std::uint64_t edge_key(int u, int v) {
    std::uint32_t left = static_cast<std::uint32_t>(std::min(u, v));
    std::uint32_t right = static_cast<std::uint32_t>(std::max(u, v));
    return (static_cast<std::uint64_t>(left) << 32) | right;
}

double hierarchy_edge_bias(int u, int v, const umat &hierarchy, const vec &cut_bias) {
    double bias = 1.0;
    for (uword level = 0; level < hierarchy.n_cols; level++) {
        if (hierarchy(u, level) != hierarchy(v, level)) {
            bias += cut_bias(level);
        }
    }
    if (!std::isfinite(bias)) {
        stop_dimer_numerical_failure("computing hierarchical cut weights");
    }
    return bias;
}

double evaluate_constraint_penalty(const uvec &plan, int n_districts, double target,
                                   const uvec &population, const Graph &graph,
                                   const List &constraints) {
    if (constraints.size() == 0) {
        return 0.0;
    }

    umat plan_matrix(plan.n_elem, 1);
    plan_matrix.col(0) = plan;
    std::vector<int> districts(n_districts);
    std::iota(districts.begin(), districts.end(), 1);

    CharacterVector constraint_names = CharacterVector::create(
        "pop_dev", "splits", "multisplits", "total_splits", "segregation", "grp_pow",
        "grp_hinge", "grp_inv_hinge", "compet", "status_quo", "incumbency", "polsby",
        "fry_hold", "log_st", "edges_removed", "qps", "custom");
    NumericVector psi(constraint_names.size());
    psi.names() = constraint_names;
    auto plan_view = plan_matrix.col(0);
    double penalty = calc_gibbs_tgt(plan_view, n_districts, static_cast<int>(plan.n_elem),
                                    districts, psi, population, target, graph, constraints);
    double global_scale = 0.0;
    if (n_districts > 0) {
        global_scale = static_cast<double>(n_districts - 1) / n_districts;
    }
    if (constraints.containsElementNamed("log_st")) {
        List instances = constraints["log_st"];
        for (int i = 0; i < instances.size(); i++) {
            List instance = instances[i];
            double strength = instance["strength"];
            if (strength != 0.0) {
                penalty += strength * global_scale *
                           eval_log_st(plan_view, graph, as<uvec>(instance["admin"]),
                                       n_districts);
            }
        }
    }
    if (constraints.containsElementNamed("edges_removed")) {
        List instances = constraints["edges_removed"];
        for (int i = 0; i < instances.size(); i++) {
            List instance = instances[i];
            double strength = instance["strength"];
            if (strength != 0.0) {
                penalty += strength * global_scale * eval_er(plan_view, graph, n_districts);
            }
        }
    }
    if (!std::isfinite(penalty)) {
        Rcpp::stop("The redist_constr penalty was not finite on a sampled plan.");
    }
    return penalty;
}

bool sample_multigraph_tree(int n_nodes, const std::vector<MultiEdge> &edges,
                            std::vector<HDimerForestEdge> &selected) {
    if (n_nodes <= 1) {
        return true;
    }

    std::vector<std::vector<int>> incident(n_nodes);
    for (int edge_id = 0; edge_id < static_cast<int>(edges.size()); edge_id++) {
        const MultiEdge &edge = edges[edge_id];
        incident[edge.from].push_back(edge_id);
        incident[edge.to].push_back(edge_id);
    }

    for (const auto &neighbors : incident) {
        if (neighbors.empty()) {
            return false;
        }
    }

    std::vector<bool> connected(n_nodes, false);
    std::vector<int> stack = {0};
    connected[0] = true;
    for (std::size_t cursor = 0; cursor < stack.size(); cursor++) {
        int current = stack[cursor];
        for (int edge_id : incident[current]) {
            const MultiEdge &edge = edges[edge_id];
            int neighbor = edge.from == current ? edge.to : edge.from;
            if (!connected[neighbor]) {
                connected[neighbor] = true;
                stack.push_back(neighbor);
            }
        }
    }
    if (static_cast<int>(stack.size()) != n_nodes) {
        return false;
    }

    std::vector<bool> visited(n_nodes, false);
    std::vector<int> next_node(n_nodes, -1);
    std::vector<int> next_edge(n_nodes, -1);
    visited[0] = true;

    for (int start = 1; start < n_nodes; start++) {
        if (visited[start]) {
            continue;
        }

        int current = start;
        std::uint64_t walk_steps = 0;
        while (!visited[current]) {
            const auto &neighbors = incident[current];
            int edge_id = neighbors[r_int(static_cast<uint32_t>(neighbors.size()))];
            const MultiEdge &edge = edges[edge_id];
            int proposal = edge.from == current ? edge.to : edge.from;
            next_node[current] = proposal;
            next_edge[current] = edge_id;
            current = proposal;
            walk_steps++;
            if (walk_steps % 1000000 == 0) {
                Rcpp::checkUserInterrupt();
            }
        }

        current = start;
        while (!visited[current]) {
            visited[current] = true;
            int edge_id = next_edge[current];
            if (edge_id < 0) {
                return false;
            }
            const MultiEdge &edge = edges[edge_id];
            selected.push_back({edge.u, edge.v});
            current = next_node[current];
        }
    }

    return true;
}

bool sample_vertex_tree(const Graph &graph, const std::vector<int> &vertices,
                        const std::vector<bool> &active,
                        std::vector<HDimerForestEdge> &selected) {
    if (vertices.size() <= 1) {
        return true;
    }

    std::unordered_map<int, int> local;
    local.reserve(vertices.size());
    for (int i = 0; i < static_cast<int>(vertices.size()); i++) {
        local[vertices[i]] = i;
    }

    std::vector<MultiEdge> edges;
    for (int u : vertices) {
        for (int v : graph[u]) {
            if (v > u && active[v] && local.find(v) != local.end()) {
                edges.push_back({local[u], local[v], u, v});
            }
        }
    }

    std::size_t before = selected.size();
    if (!sample_multigraph_tree(static_cast<int>(vertices.size()), edges, selected)) {
        return false;
    }
    return selected.size() - before == vertices.size() - 1;
}

bool connect_hierarchy_level(const Graph &graph, const std::vector<int> &vertices,
                             const std::vector<bool> &active, const umat &hierarchy,
                             int parent_level, int child_level,
                             std::vector<HDimerForestEdge> &selected) {
    std::map<unsigned int, std::vector<int>> parent_vertices;
    for (int vertex : vertices) {
        unsigned int parent = 0;
        if (parent_level >= 0) {
            parent = hierarchy(vertex, parent_level);
        }
        parent_vertices[parent].push_back(vertex);
    }

    for (const auto &entry : parent_vertices) {
        const std::vector<int> &inside = entry.second;
        std::map<unsigned int, int> child_local;
        for (int vertex : inside) {
            unsigned int child = hierarchy(vertex, child_level);
            if (child_local.find(child) == child_local.end()) {
                child_local[child] = static_cast<int>(child_local.size());
            }
        }

        int n_children = static_cast<int>(child_local.size());
        if (n_children <= 1) {
            continue;
        }

        std::vector<MultiEdge> edges;
        for (int u : inside) {
            int from = child_local[hierarchy(u, child_level)];
            for (int v : graph[u]) {
                if (v <= u || !active[v]) {
                    continue;
                }
                if (parent_level >= 0 &&
                    hierarchy(v, parent_level) != hierarchy(u, parent_level)) {
                    continue;
                }
                unsigned int child_v = hierarchy(v, child_level);
                if (child_v == hierarchy(u, child_level)) {
                    continue;
                }
                auto found = child_local.find(child_v);
                if (found != child_local.end()) {
                    edges.push_back({from, found->second, u, v});
                }
            }
        }

        std::size_t before = selected.size();
        if (!sample_multigraph_tree(n_children, edges, selected)) {
            return false;
        }
        if (selected.size() - before != static_cast<std::size_t>(n_children - 1)) {
            return false;
        }
    }

    return true;
}

bool sample_hierarchical_tree(const Graph &graph, const std::vector<int> &vertices,
                              const umat &hierarchy, std::vector<HDimerForestEdge> &selected) {
    int n_vertices = static_cast<int>(graph.size());
    std::vector<bool> active(n_vertices, false);
    for (int vertex : vertices) {
        active[vertex] = true;
    }

    if (hierarchy.n_cols == 0) {
        return sample_vertex_tree(graph, vertices, active, selected);
    }

    // A hierarchy-compatible tree has a unique product decomposition: trees
    // inside the finest units, then multigraph trees joining children inside
    // each parent. Uniform Wilson draws for every factor therefore give a
    // uniform draw from all hierarchy-compatible trees on the district.
    int finest = static_cast<int>(hierarchy.n_cols) - 1;
    std::map<unsigned int, std::vector<int>> finest_vertices;
    for (int vertex : vertices) {
        finest_vertices[hierarchy(vertex, finest)].push_back(vertex);
    }
    for (const auto &entry : finest_vertices) {
        if (!sample_vertex_tree(graph, entry.second, active, selected)) {
            return false;
        }
    }

    for (int parent = finest - 1; parent >= 0; parent--) {
        if (!connect_hierarchy_level(graph, vertices, active, hierarchy, parent, parent + 1,
                                     selected)) {
            return false;
        }
    }
    return connect_hierarchy_level(graph, vertices, active, hierarchy, -1, 0, selected);
}

void fill_forest_adjacency(int n_vertices, const std::vector<HDimerForestEdge> &forest,
                           ForestAdjacency &adjacency) {
    if (static_cast<int>(adjacency.size()) != n_vertices) {
        adjacency.resize(n_vertices);
    }
    for (auto &neighbors : adjacency) {
        neighbors.clear();
    }
    for (int edge_id = 0; edge_id < static_cast<int>(forest.size()); edge_id++) {
        const HDimerForestEdge &edge = forest[edge_id];
        adjacency[edge.u].push_back({edge.v, edge_id});
        adjacency[edge.v].push_back({edge.u, edge_id});
    }
}

void fill_forest_adjacency_excluding(int n_vertices,
                                     const std::vector<HDimerForestEdge> &forest,
                                     const std::vector<int> &cut_edge_marks, int cut_edge_mark,
                                     ForestAdjacency &adjacency) {
    if (static_cast<int>(adjacency.size()) != n_vertices) {
        adjacency.resize(n_vertices);
    }
    for (auto &neighbors : adjacency) {
        neighbors.clear();
    }
    for (int edge_id = 0; edge_id < static_cast<int>(forest.size()); edge_id++) {
        if (cut_edge_marks[edge_id] == cut_edge_mark) {
            continue;
        }
        const HDimerForestEdge &edge = forest[edge_id];
        adjacency[edge.u].push_back({edge.v, edge_id});
        adjacency[edge.v].push_back({edge.u, edge_id});
    }
}

bool fill_district_vertices(const uvec &plan, int n_districts,
                            std::vector<std::vector<int>> &district_vertices) {
    for (auto &vertices : district_vertices) {
        vertices.clear();
    }
    for (int vertex = 0; vertex < static_cast<int>(plan.n_elem); vertex++) {
        int district = static_cast<int>(plan(vertex)) - 1;
        if (district < 0 || district >= n_districts) {
            return false;
        }
        district_vertices[district].push_back(vertex);
    }
    return true;
}

bool choose_cut_edges(const std::vector<HDimerForestEdge> &forest,
                      const ForestAdjacency &adjacency, const uvec &plan, int n_districts,
                      const vec &population, const umat &hierarchy, const vec &cut_bias,
                      const std::vector<bool> &selected_districts, HDimerWorkspace &workspace,
                      std::vector<int> &cut_ids) {
    std::fill(workspace.roots.begin(), workspace.roots.end(), -1);
    std::fill(workspace.district_sizes.begin(), workspace.district_sizes.end(), 0);
    std::fill(workspace.parent.begin(), workspace.parent.end(), -2);
    std::fill(workspace.parent_edge.begin(), workspace.parent_edge.end(), -1);
    std::fill(workspace.subtree.begin(), workspace.subtree.end(), 0.0);
    for (auto &order : workspace.district_order) {
        order.clear();
    }
    std::vector<int> &roots = workspace.roots;
    std::vector<int> &district_sizes = workspace.district_sizes;
    std::vector<int> &parent = workspace.parent;
    std::vector<int> &parent_edge = workspace.parent_edge;
    std::vector<double> &subtree = workspace.subtree;
    for (int district = 0; district < n_districts; district++) {
        if (!selected_districts[district]) {
            continue;
        }
        const std::vector<int> &vertices = workspace.district_vertices[district];
        if (!vertices.empty()) {
            roots[district] = vertices.front();
            district_sizes[district] = static_cast<int>(vertices.size());
        }
    }

    cut_ids.clear();
    cut_ids.reserve(std::count(selected_districts.begin(), selected_districts.end(), true));
    for (int district = 0; district < n_districts; district++) {
        if (!selected_districts[district]) {
            continue;
        }
        if (district_sizes[district] <= 1 || roots[district] < 0) {
            return false;
        }

        int root = roots[district];
        std::vector<int> &district_order = workspace.district_order[district];
        district_order.reserve(district_sizes[district]);
        parent[root] = -1;
        district_order.push_back(root);
        for (std::size_t cursor = 0; cursor < district_order.size(); cursor++) {
            int u = district_order[cursor];
            for (const auto &[v, edge_id] : adjacency[u]) {
                if (parent[v] != -2 || static_cast<int>(plan(v)) - 1 != district) {
                    continue;
                }
                parent[v] = u;
                parent_edge[v] = edge_id;
                district_order.push_back(v);
            }
        }
        if (static_cast<int>(district_order.size()) != district_sizes[district]) {
            return false;
        }

        for (auto it = district_order.rbegin(); it != district_order.rend(); ++it) {
            int u = *it;
            subtree[u] += population(u);
            if (parent[u] >= 0) {
                subtree[parent[u]] += subtree[u];
            }
        }

        double total = subtree[root];
        if (!std::isfinite(total)) {
            stop_dimer_numerical_failure("computing district populations for a tree cut");
        }
        double weight_sum = 0.0;
        for (int u : district_order) {
            if (parent[u] < 0) {
                continue;
            }
            int edge_id = parent_edge[u];
            const HDimerForestEdge &edge = forest[edge_id];
            double bias = hierarchy_edge_bias(edge.u, edge.v, hierarchy, cut_bias);
            double mass = subtree[u] * (total - subtree[u]) + CUT_SUPPORT;
            double weight = bias * mass;
            if (!std::isfinite(weight) || !(weight > 0.0)) {
                stop_dimer_numerical_failure("computing tree-cut probabilities");
            }
            weight_sum += weight;
        }
        if (!std::isfinite(weight_sum) || !(weight_sum > 0.0)) {
            stop_dimer_numerical_failure("normalizing tree-cut probabilities");
        }
        double draw = r_unif() * weight_sum;
        double cumulative = 0.0;
        int selected_edge = -1;
        for (int u : district_order) {
            if (parent[u] < 0) {
                continue;
            }
            int edge_id = parent_edge[u];
            const HDimerForestEdge &edge = forest[edge_id];
            double bias = hierarchy_edge_bias(edge.u, edge.v, hierarchy, cut_bias);
            double mass = subtree[u] * (total - subtree[u]) + CUT_SUPPORT;
            double weight = bias * mass;
            if (!std::isfinite(weight) || !(weight > 0.0)) {
                stop_dimer_numerical_failure("sampling a tree cut");
            }
            cumulative += weight;
            if (!std::isfinite(cumulative)) {
                stop_dimer_numerical_failure("sampling a tree cut");
            }
            selected_edge = edge_id;
            if (draw <= cumulative) {
                break;
            }
        }
        if (selected_edge < 0) {
            return false;
        }
        cut_ids.push_back(selected_edge);
    }

    return true;
}

} // namespace

bool sample_hierarchical_forest(const Graph &graph, const uvec &plan, int n_districts,
                                const umat &hierarchy, std::vector<HDimerForestEdge> &forest) {
    std::vector<std::vector<int>> district_vertices(n_districts);
    for (int vertex = 0; vertex < static_cast<int>(plan.n_elem); vertex++) {
        int district = static_cast<int>(plan(vertex)) - 1;
        if (district < 0 || district >= n_districts) {
            return false;
        }
        district_vertices[district].push_back(vertex);
    }

    forest.clear();
    forest.reserve(plan.n_elem - n_districts);
    for (const auto &vertices : district_vertices) {
        if (vertices.empty() || !sample_hierarchical_tree(graph, vertices, hierarchy, forest)) {
            return false;
        }
    }
    return forest.size() == plan.n_elem - n_districts;
}

namespace {

bool sample_hierarchical_forest_subset(const Graph &graph, const uvec &plan, int n_districts,
                                       const umat &hierarchy,
                                       const std::vector<bool> &selected_districts,
                                       const std::vector<std::vector<int>> &district_vertices,
                                       std::vector<HDimerForestEdge> &forest) {
    std::vector<HDimerForestEdge> refreshed;
    refreshed.reserve(plan.n_elem - n_districts);
    for (const HDimerForestEdge &edge : forest) {
        int left_district = static_cast<int>(plan(edge.u)) - 1;
        int right_district = static_cast<int>(plan(edge.v)) - 1;
        if (left_district < 0 || left_district >= n_districts ||
            right_district != left_district) {
            return false;
        }
        if (!selected_districts[left_district]) {
            refreshed.push_back(edge);
        }
    }
    for (int district = 0; district < n_districts; district++) {
        if (!selected_districts[district]) {
            continue;
        }
        const std::vector<int> &vertices = district_vertices[district];
        if (vertices.empty() ||
            !sample_hierarchical_tree(graph, vertices, hierarchy, refreshed)) {
            return false;
        }
    }
    if (refreshed.size() != plan.n_elem - n_districts) {
        return false;
    }
    forest = std::move(refreshed);
    return true;
}

bool build_fragments(const std::vector<HDimerForestEdge> &forest,
                     const ForestAdjacency &adjacency, const uvec &plan, const vec &population,
                     const umat &hierarchy, const vec &cut_bias,
                     const std::vector<bool> &selected_districts, HDimerWorkspace &workspace) {
    int n_vertices = static_cast<int>(population.n_elem);
    std::vector<int> &fragment_of = workspace.fragment_of;
    std::vector<Fragment> &fragments = workspace.fragments;
    std::vector<double> &reroot_sum = workspace.reroot_sum;
    std::vector<int> &parent = workspace.parent;
    std::vector<int> &parent_edge = workspace.parent_edge;
    std::vector<double> &subtree = workspace.subtree;
    std::fill(fragment_of.begin(), fragment_of.end(), -1);
    fragments.clear();
    std::fill(reroot_sum.begin(), reroot_sum.end(), 0.0);
    std::fill(parent.begin(), parent.end(), -2);
    std::fill(parent_edge.begin(), parent_edge.end(), -1);
    std::fill(subtree.begin(), subtree.end(), 0.0);

    for (int start = 0; start < n_vertices; start++) {
        if (fragment_of[start] >= 0) {
            continue;
        }
        int start_district = static_cast<int>(plan(start)) - 1;
        if (start_district < 0 ||
            start_district >= static_cast<int>(selected_districts.size()) ||
            !selected_districts[start_district]) {
            continue;
        }
        int fragment_id = static_cast<int>(fragments.size());
        Fragment fragment;
        fragment.units.resize(hierarchy.n_cols);
        fragment.old_district = static_cast<int>(plan(start)) - 1;

        parent[start] = -1;
        fragment_of[start] = fragment_id;
        fragment.vertices.push_back(start);
        for (std::size_t cursor = 0; cursor < fragment.vertices.size(); cursor++) {
            int u = fragment.vertices[cursor];
            if (static_cast<int>(plan(u)) - 1 != fragment.old_district) {
                return false;
            }
            fragment.population += population(u);
            if (!std::isfinite(fragment.population)) {
                stop_dimer_numerical_failure("computing fragment populations");
            }
            for (const auto &[v, edge_id] : adjacency[u]) {
                if (fragment_of[v] < 0) {
                    fragment_of[v] = fragment_id;
                    parent[v] = u;
                    parent_edge[v] = edge_id;
                    fragment.vertices.push_back(v);
                }
            }
        }

        for (uword level = 0; level < hierarchy.n_cols; level++) {
            std::vector<int> units;
            units.reserve(fragment.vertices.size());
            for (int vertex : fragment.vertices) {
                units.push_back(static_cast<int>(hierarchy(vertex, level)));
            }
            std::sort(units.begin(), units.end());
            units.erase(std::unique(units.begin(), units.end()), units.end());
            fragment.units[level] = std::move(units);
        }

        int root = fragment.vertices.front();
        for (auto it = fragment.vertices.rbegin(); it != fragment.vertices.rend(); ++it) {
            int u = *it;
            subtree[u] += population(u);
            if (parent[u] >= 0) {
                subtree[parent[u]] += subtree[u];
            }
        }

        double root_sum = 0.0;
        for (int u : fragment.vertices) {
            if (parent[u] < 0) {
                continue;
            }
            const HDimerForestEdge &edge = forest[parent_edge[u]];
            double bias = hierarchy_edge_bias(edge.u, edge.v, hierarchy, cut_bias);
            double cut_weight =
                bias * (subtree[u] * (fragment.population - subtree[u]) + CUT_SUPPORT);
            double reroot_weight = bias * subtree[u];
            if (!std::isfinite(cut_weight) || !std::isfinite(reroot_weight)) {
                stop_dimer_numerical_failure("constructing fragment cut weights");
            }
            fragment.cut_normalizer += cut_weight;
            root_sum += reroot_weight;
            if (!std::isfinite(fragment.cut_normalizer) || !std::isfinite(root_sum)) {
                stop_dimer_numerical_failure("accumulating fragment cut weights");
            }
        }
        reroot_sum[root] = root_sum;
        for (int u : fragment.vertices) {
            if (parent[u] < 0) {
                continue;
            }
            const HDimerForestEdge &edge = forest[parent_edge[u]];
            double bias = hierarchy_edge_bias(edge.u, edge.v, hierarchy, cut_bias);
            reroot_sum[u] =
                reroot_sum[parent[u]] + bias * (fragment.population - 2.0 * subtree[u]);
            if (!std::isfinite(reroot_sum[u])) {
                stop_dimer_numerical_failure("rerooting fragment cut weights");
            }
        }

        fragments.push_back(std::move(fragment));
    }

    return true;
}

bool initialize_pair_cache(const Fragment &left, const Fragment &right, const umat &hierarchy,
                           const vec &split_penalty, double boundary_penalty, double lower,
                           double upper, PairCache &cache) {
    double combined_population = left.population + right.population;
    if (!std::isfinite(combined_population)) {
        stop_dimer_numerical_failure("combining fragment populations");
    }
    cache.population_valid = combined_population >= lower && combined_population <= upper;
    cache.hierarchy_valid = true;
    cache.log_target = boundary_penalty * cache.boundary_edges;
    if (!std::isfinite(cache.log_target)) {
        stop_dimer_numerical_failure("computing connector target weights");
    }
    cache.shared_units.assign(hierarchy.n_cols, -1);
    for (uword level = 0; level < hierarchy.n_cols; level++) {
        const std::vector<int> &left_units = left.units[level];
        const std::vector<int> &right_units = right.units[level];
        std::size_t left_at = 0;
        std::size_t right_at = 0;
        int shared_count = 0;
        int shared_unit = -1;
        while (left_at < left_units.size() && right_at < right_units.size()) {
            if (left_units[left_at] < right_units[right_at]) {
                left_at++;
            } else if (right_units[right_at] < left_units[left_at]) {
                right_at++;
            } else {
                shared_count++;
                shared_unit = left_units[left_at];
                left_at++;
                right_at++;
                if (shared_count > 1) {
                    cache.hierarchy_valid = false;
                    return false;
                }
            }
        }
        cache.shared_units[level] = shared_unit;
        int union_count =
            static_cast<int>(left_units.size() + right_units.size()) - shared_count;
        cache.log_target -= split_penalty(level) * union_count;
        if (!std::isfinite(cache.log_target)) {
            stop_dimer_numerical_failure("computing connector target weights");
        }
    }
    return cache.population_valid && cache.hierarchy_valid;
}

bool hierarchy_join_valid(const PairCache &cache, int u, int v, const umat &hierarchy) {
    for (uword level = 0; level < hierarchy.n_cols; level++) {
        int shared_unit = cache.shared_units[level];
        if (shared_unit >= 0 && (static_cast<int>(hierarchy(u, level)) != shared_unit ||
                                 static_cast<int>(hierarchy(v, level)) != shared_unit)) {
            return false;
        }
    }
    return true;
}

std::vector<CandidateEdge> build_candidate_edges(
    const std::vector<BaseEdge> &base_edges, const std::vector<int> &fragment_of,
    const std::vector<Fragment> &fragments, const std::vector<double> &reroot_sum,
    const umat &hierarchy, const vec &split_penalty, double boundary_penalty, double lower,
    double upper, const std::unordered_map<std::uint64_t, int> &old_cut_slots,
    std::vector<int> &old_candidate_ids) {
    int n_fragments = static_cast<int>(fragments.size());
    std::size_t table_size = static_cast<std::size_t>(n_fragments) * n_fragments;
    std::vector<int> pair_cache_index(table_size, -1);
    std::vector<std::pair<int, int>> touched_pairs;
    std::vector<PairCache> pair_caches;
    for (const BaseEdge &edge : base_edges) {
        int left = fragment_of[edge.u];
        int right = fragment_of[edge.v];
        if (left < 0 || right < 0 || left == right) {
            continue;
        }
        if (left > right) {
            std::swap(left, right);
        }
        std::size_t pair_id = static_cast<std::size_t>(left) * n_fragments + right;
        int cache_id = pair_cache_index[pair_id];
        if (cache_id < 0) {
            cache_id = static_cast<int>(pair_caches.size());
            pair_cache_index[pair_id] = cache_id;
            pair_caches.emplace_back();
            touched_pairs.push_back({left, right});
        }
        pair_caches[cache_id].boundary_edges++;
    }

    for (int cache_id = 0; cache_id < static_cast<int>(pair_caches.size()); cache_id++) {
        const auto &[left, right] = touched_pairs[cache_id];
        initialize_pair_cache(fragments[left], fragments[right], hierarchy, split_penalty,
                              boundary_penalty, lower, upper, pair_caches[cache_id]);
    }

    std::vector<CandidateEdge> candidates;
    candidates.reserve(base_edges.size());
    old_candidate_ids.assign(old_cut_slots.size(), -1);
    for (const BaseEdge &edge : base_edges) {
        int left_id = fragment_of[edge.u];
        int right_id = fragment_of[edge.v];
        if (left_id < 0 || right_id < 0 || left_id == right_id) {
            continue;
        }
        int u = edge.u;
        int v = edge.v;
        if (left_id > right_id) {
            std::swap(left_id, right_id);
            std::swap(u, v);
        }
        std::size_t pair_id = static_cast<std::size_t>(left_id) * n_fragments + right_id;
        const PairCache &cache = pair_caches[pair_cache_index[pair_id]];
        if (!cache.population_valid || !cache.hierarchy_valid ||
            !hierarchy_join_valid(cache, u, v, hierarchy)) {
            continue;
        }

        const Fragment &left = fragments[left_id];
        const Fragment &right = fragments[right_id];
        double numerator = edge.cut_bias * (left.population * right.population + CUT_SUPPORT);
        if (!std::isfinite(numerator) || !(numerator > 0.0)) {
            stop_dimer_numerical_failure("computing reverse cut probabilities");
        }
        // Reconstruct the cut-law normalizer for the proposed joined tree
        // in O(1). The reroot sums account for the other fragment's
        // population on every internal edge.
        double denominator = left.cut_normalizer + right.population * reroot_sum[u] +
                             right.cut_normalizer + left.population * reroot_sum[v] + numerator;
        if (!(denominator > 0.0) || !std::isfinite(denominator)) {
            stop_dimer_numerical_failure("normalizing reverse cut probabilities");
        }

        // This is target(plan, tree) times the reverse auxiliary-cut
        // probability. A product of these edge weights is therefore the
        // exact conditional law on global fragment matchings.
        double log_weight = std::log(numerator) - std::log(denominator) + cache.log_target;
        if (!std::isfinite(log_weight)) {
            stop_dimer_numerical_failure("constructing connector matching weights");
        }
        int candidate_id = static_cast<int>(candidates.size());
        candidates.push_back({u, v, left_id, right_id, log_weight});
        auto old_cut = old_cut_slots.find(edge_key(u, v));
        if (old_cut != old_cut_slots.end()) {
            old_candidate_ids[old_cut->second] = candidate_id;
        }
    }

    return candidates;
}

bool relabel_selected_fragments(const std::vector<Fragment> &fragments,
                                const std::vector<int> &fragment_of,
                                std::vector<std::pair<int, int>> component_pairs,
                                const std::vector<int> &selected_labels, const uvec &old_plan,
                                int n_districts, uvec &new_plan) {
    int n_selected = static_cast<int>(selected_labels.size());
    if (static_cast<int>(component_pairs.size()) != n_selected) {
        return false;
    }
    std::vector<int> selected_slot(n_districts, -1);
    for (int slot = 0; slot < n_selected; slot++) {
        int label = selected_labels[slot];
        if (label < 0 || label >= n_districts || selected_slot[label] >= 0) {
            return false;
        }
        selected_slot[label] = slot;
    }
    std::sort(component_pairs.begin(), component_pairs.end(),
              [](const auto &left, const auto &right) {
                  return std::min(left.first, left.second) <
                         std::min(right.first, right.second);
              });

    std::vector<std::vector<double>> cost(n_selected, std::vector<double>(n_selected, 0.0));
    std::vector<std::vector<double>> overlap(n_selected, std::vector<double>(n_selected, 0.0));
    std::vector<int> component_of_fragment(fragments.size(), -1);
    double max_overlap = 0.0;
    for (int component = 0; component < n_selected; component++) {
        const auto &[left, right] = component_pairs[component];
        if (left < 0 || right < 0 || left >= static_cast<int>(fragments.size()) ||
            right >= static_cast<int>(fragments.size()) || left == right ||
            component_of_fragment[left] >= 0 || component_of_fragment[right] >= 0) {
            return false;
        }
        component_of_fragment[left] = component;
        component_of_fragment[right] = component;
        for (int fragment_id : {left, right}) {
            const Fragment &fragment = fragments[fragment_id];
            int slot = fragment.old_district < 0 || fragment.old_district >= n_districts
                           ? -1
                           : selected_slot[fragment.old_district];
            if (slot < 0) {
                return false;
            }
            overlap[component][slot] += fragment.population;
            max_overlap = std::max(max_overlap, overlap[component][slot]);
        }
    }
    for (int component = 0; component < n_selected; component++) {
        for (int slot = 0; slot < n_selected; slot++) {
            cost[component][slot] = max_overlap - overlap[component][slot];
        }
    }

    HungarianAlgorithm solver;
    std::vector<int> assignment;
    solver.Solve(cost, assignment);
    if (static_cast<int>(assignment.size()) != n_selected) {
        return false;
    }

    new_plan = old_plan;
    for (int vertex = 0; vertex < static_cast<int>(old_plan.n_elem); vertex++) {
        int fragment_id = fragment_of[vertex];
        if (fragment_id < 0) {
            continue;
        }
        int component = component_of_fragment[fragment_id];
        if (component < 0 || assignment[component] < 0 || assignment[component] >= n_selected) {
            return false;
        }
        new_plan(vertex) = selected_labels[assignment[component]] + 1;
    }
    return true;
}

std::vector<int> sample_update_districts(int n_districts, int ell) {
    std::vector<int> districts(n_districts);
    std::iota(districts.begin(), districts.end(), 0);
    if (ell == n_districts) {
        return districts;
    }
    for (int slot = 0; slot < ell; slot++) {
        int other = slot + r_int(static_cast<uint32_t>(n_districts - slot));
        std::swap(districts[slot], districts[other]);
    }
    districts.resize(ell);
    std::sort(districts.begin(), districts.end());
    return districts;
}

void selected_base_edges(const std::vector<std::vector<BaseEdge>> &edges_by_lower_vertex,
                         const uvec &plan, const std::vector<bool> &selected_districts,
                         const std::vector<int> &selected_labels,
                         const std::vector<std::vector<int>> &district_vertices,
                         std::vector<BaseEdge> &output) {
    output.clear();
    for (int district : selected_labels) {
        for (int vertex : district_vertices[district]) {
            for (const BaseEdge &edge : edges_by_lower_vertex[vertex]) {
                int other_district = static_cast<int>(plan(edge.v)) - 1;
                if (other_district >= 0 &&
                    other_district < static_cast<int>(selected_districts.size()) &&
                    selected_districts[other_district]) {
                    output.push_back(edge);
                }
            }
        }
    }
}

HDimerStepDiagnostics
run_dimer_step(const std::vector<BaseEdge> &base_edges,
               const std::vector<std::vector<BaseEdge>> &edges_by_lower_vertex,
               const vec &population, const umat &hierarchy, const vec &split_penalty,
               const vec &cut_bias, double boundary_penalty, double lower, double upper,
               int n_districts, uvec &plan, double &constraint_penalty,
               const uvec &constraint_population, double constraint_target, const Graph &graph,
               const List &constraints, std::vector<HDimerForestEdge> &forest,
               const std::vector<int> &selected_labels, HDimerWorkspace &workspace) {
    HDimerStepDiagnostics diagnostics;
    std::vector<bool> selected_districts(n_districts, false);
    for (int label : selected_labels) {
        if (label < 0 || label >= n_districts || selected_districts[label]) {
            diagnostics.no_matching = true;
            return diagnostics;
        }
        selected_districts[label] = true;
    }
    std::vector<BaseEdge> selected_edges;
    const std::vector<BaseEdge> *candidate_base_edges = &base_edges;
    if (static_cast<int>(selected_labels.size()) < n_districts) {
        selected_base_edges(edges_by_lower_vertex, plan, selected_districts, selected_labels,
                            workspace.district_vertices, selected_edges);
        candidate_base_edges = &selected_edges;
    }
    Clock::time_point phase_start = Clock::now();
    fill_forest_adjacency(static_cast<int>(population.n_elem), forest,
                          workspace.full_adjacency);
    std::vector<int> cut_ids;
    if (!choose_cut_edges(forest, workspace.full_adjacency, plan, n_districts, population,
                          hierarchy, cut_bias, selected_districts, workspace, cut_ids)) {
        diagnostics.cut_seconds = elapsed_seconds(phase_start);
        diagnostics.no_matching = true;
        return diagnostics;
    }

    diagnostics.cut_seconds = elapsed_seconds(phase_start);

    phase_start = Clock::now();
    if (++workspace.cut_edge_mark == std::numeric_limits<int>::max()) {
        std::fill(workspace.cut_edge_marks.begin(), workspace.cut_edge_marks.end(), 0);
        workspace.cut_edge_mark = 1;
    }
    for (int edge_id : cut_ids) {
        workspace.cut_edge_marks[edge_id] = workspace.cut_edge_mark;
    }
    fill_forest_adjacency_excluding(static_cast<int>(population.n_elem), forest,
                                    workspace.cut_edge_marks, workspace.cut_edge_mark,
                                    workspace.half_adjacency);
    if (!build_fragments(forest, workspace.half_adjacency, plan, population, hierarchy,
                         cut_bias, selected_districts, workspace) ||
        static_cast<int>(workspace.fragments.size()) !=
            2 * static_cast<int>(selected_labels.size())) {
        diagnostics.fragment_seconds = elapsed_seconds(phase_start);
        diagnostics.no_matching = true;
        return diagnostics;
    }
    diagnostics.fragment_seconds = elapsed_seconds(phase_start);
    const std::vector<int> &fragment_of = workspace.fragment_of;
    const std::vector<Fragment> &fragments = workspace.fragments;
    const std::vector<double> &reroot_sum = workspace.reroot_sum;

    phase_start = Clock::now();
    std::unordered_map<std::uint64_t, int> old_cut_slots;
    old_cut_slots.reserve(2 * cut_ids.size());
    for (int slot = 0; slot < static_cast<int>(cut_ids.size()); slot++) {
        const HDimerForestEdge &edge = forest[cut_ids[slot]];
        old_cut_slots[edge_key(edge.u, edge.v)] = slot;
    }
    std::vector<int> old_candidate_ids;
    std::vector<CandidateEdge> candidates = build_candidate_edges(
        *candidate_base_edges, fragment_of, fragments, reroot_sum, hierarchy, split_penalty,
        boundary_penalty, lower, upper, old_cut_slots, old_candidate_ids);
    std::vector<HDimerMatchingEdge> matching_edges;
    matching_edges.reserve(candidates.size());
    for (int candidate_id = 0; candidate_id < static_cast<int>(candidates.size());
         candidate_id++) {
        const CandidateEdge &candidate = candidates[candidate_id];
        matching_edges.push_back(
            {candidate.fragment_u, candidate.fragment_v, candidate.log_weight, candidate_id});
    }

    int n_fragments = static_cast<int>(fragments.size());
    std::vector<bool> old_fragment_pairs(static_cast<std::size_t>(n_fragments) * n_fragments,
                                         false);
    std::vector<int> old_fragment_degree(n_fragments, 0);
    double old_log_weight = 0.0;
    for (int slot = 0; slot < static_cast<int>(cut_ids.size()); slot++) {
        int candidate_id = old_candidate_ids[slot];
        if (candidate_id < 0 || candidate_id >= static_cast<int>(candidates.size())) {
            diagnostics.candidate_seconds = elapsed_seconds(phase_start);
            diagnostics.no_matching = true;
            return diagnostics;
        }
        old_log_weight += candidates[candidate_id].log_weight;
        if (!std::isfinite(old_log_weight)) {
            stop_dimer_numerical_failure("accumulating matching weights");
        }
        int left = candidates[candidate_id].fragment_u;
        int right = candidates[candidate_id].fragment_v;
        if (left > right) {
            std::swap(left, right);
        }
        old_fragment_pairs[static_cast<std::size_t>(left) * n_fragments + right] = true;
        old_fragment_degree[left]++;
        old_fragment_degree[right]++;
    }
    if (std::any_of(old_fragment_degree.begin(), old_fragment_degree.end(),
                    [](int degree) { return degree != 1; })) {
        diagnostics.candidate_seconds = elapsed_seconds(phase_start);
        diagnostics.no_matching = true;
        return diagnostics;
    }
    diagnostics.candidate_seconds = elapsed_seconds(phase_start);

    phase_start = Clock::now();
    HDimerMatchingResult matching = sample_weighted_perfect_matching(
        static_cast<int>(fragments.size()), matching_edges, true);
    diagnostics.matching_seconds = elapsed_seconds(phase_start);
    diagnostics.planar = matching.status != HDimerMatchingStatus::nonplanar;
    diagnostics.numerical_failure = matching.status == HDimerMatchingStatus::numerical_failure;
    diagnostics.used_subset_dp = matching.used_subset_dp;
    if (matching.status == HDimerMatchingStatus::numerical_failure) {
        Rcpp::stop(
            "Numerical failure while sampling a perfect matching; no chain was returned.");
    }
    if (matching.status == HDimerMatchingStatus::invalid_input ||
        matching.status == HDimerMatchingStatus::no_perfect_matching) {
        Rcpp::stop("The matching backend rejected a candidate graph with a known perfect "
                   "matching.");
    }
    if (matching.status == HDimerMatchingStatus::nonplanar) {
        return diagnostics;
    }
    diagnostics.old_matching_probability =
        std::exp(std::min(0.0, old_log_weight - matching.log_partition));

    phase_start = Clock::now();
    std::vector<HDimerForestEdge> proposed_forest;
    proposed_forest.reserve(forest.size());
    for (int edge_id = 0; edge_id < static_cast<int>(forest.size()); edge_id++) {
        if (workspace.cut_edge_marks[edge_id] != workspace.cut_edge_mark) {
            proposed_forest.push_back(forest[edge_id]);
        }
    }
    std::vector<std::pair<int, int>> component_pairs;
    component_pairs.reserve(selected_labels.size());
    for (int selected_id : matching.selected_edge_indices) {
        if (selected_id < 0 || selected_id >= static_cast<int>(candidates.size())) {
            Rcpp::stop("The matching backend returned an invalid connector index.");
        }
        const CandidateEdge &candidate = candidates[selected_id];
        proposed_forest.push_back({candidate.u, candidate.v});
        int left = std::min(candidate.fragment_u, candidate.fragment_v);
        int right = std::max(candidate.fragment_u, candidate.fragment_v);
        component_pairs.push_back({left, right});
        if (!old_fragment_pairs[static_cast<std::size_t>(left) * n_fragments + right]) {
            diagnostics.matching_edges_changed++;
        }
    }
    if (static_cast<int>(proposed_forest.size()) !=
            static_cast<int>(population.n_elem) - n_districts ||
        static_cast<int>(component_pairs.size()) != static_cast<int>(selected_labels.size())) {
        Rcpp::stop("The matching backend returned an incomplete perfect matching.");
    }

    uvec proposed_plan;
    if (diagnostics.matching_edges_changed == 0) {
        proposed_plan = plan;
    } else if (!relabel_selected_fragments(fragments, fragment_of, component_pairs,
                                           selected_labels, plan, n_districts, proposed_plan)) {
        Rcpp::stop("Failed to relabel a sampled hierarchical split-dimer plan.");
    }

    if (constraints.size() > 0) {
        double proposed_constraint_penalty = evaluate_constraint_penalty(
            proposed_plan, n_districts, constraint_target, constraint_population, graph,
            constraints);
        double log_accept = constraint_penalty - proposed_constraint_penalty;
        if (log_accept < 0.0 && std::log(r_unif()) > log_accept) {
            diagnostics.constraint_rejected = true;
            diagnostics.relabel_seconds = elapsed_seconds(phase_start);
            return diagnostics;
        }
        constraint_penalty = proposed_constraint_penalty;
    }

    for (int vertex = 0; vertex < static_cast<int>(plan.n_elem); vertex++) {
        if (proposed_plan(vertex) != plan(vertex)) {
            diagnostics.vertices_changed++;
        }
    }
    diagnostics.changed = diagnostics.vertices_changed > 0;
    plan = std::move(proposed_plan);
    forest = std::move(proposed_forest);
    diagnostics.relabel_seconds = elapsed_seconds(phase_start);
    return diagnostics;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List hier_dimer_plans(int N, Rcpp::List adj, const arma::uvec &init, const arma::vec &pop,
                            const arma::umat &hierarchy, const arma::vec &split_penalty,
                            const arma::vec &cut_bias, double boundary_penalty, double target,
                            double lower, double upper, Rcpp::List constraints, int thin,
                            int refresh, int ell, int verbosity) {
    seed_rng(static_cast<int>(Rcpp::sample(INT_MAX, 1)[0]));
    Graph graph = list_to_graph(adj);
    int n_vertices = static_cast<int>(graph.size());
    if (N < 1 || thin < 1 || refresh < 1) {
        Rcpp::stop("N, thin, and refresh must be positive integers.");
    }
    if (init.n_elem == 0 || init.n_elem != static_cast<uword>(n_vertices) ||
        pop.n_elem != static_cast<uword>(n_vertices)) {
        Rcpp::stop("Initial plan and population must have one value per map unit.");
    }
    if (!pop.is_finite() || arma::any(pop < 0.0)) {
        Rcpp::stop("Population must be finite and non-negative.");
    }
    int n_districts = static_cast<int>(init.max());
    arma::uvec district_labels = arma::unique(init);
    if (init.min() != 1 || district_labels.n_elem != static_cast<uword>(n_districts)) {
        Rcpp::stop("Initial district labels must be consecutive positive integers.");
    }
    if (ell < 2 || ell > n_districts) {
        Rcpp::stop("l must be between two and the number of districts.");
    }
    if (hierarchy.n_rows != static_cast<uword>(n_vertices)) {
        Rcpp::stop("Hierarchy must have one row per map unit.");
    }
    if (split_penalty.n_elem != hierarchy.n_cols || cut_bias.n_elem != hierarchy.n_cols) {
        Rcpp::stop("Split penalties and cut biases must match the hierarchy levels.");
    }
    if (!split_penalty.is_finite() || !cut_bias.is_finite() || arma::any(split_penalty < 0.0) ||
        arma::any(cut_bias < 0.0) || !std::isfinite(boundary_penalty) ||
        boundary_penalty < 0.0 || !std::isfinite(target) || target <= 0.0 ||
        !std::isfinite(lower) || !std::isfinite(upper) || lower < 0.0 || upper < 0.0 ||
        lower > upper) {
        Rcpp::stop("Weights and population bounds must be finite and non-negative.");
    }

    uvec constraint_population;
    if (constraints.size() > 0) {
        constraint_population.set_size(pop.n_elem);
        for (uword i = 0; i < pop.n_elem; i++) {
            double rounded = std::round(pop(i));
            if (!std::isfinite(pop(i)) || pop(i) < 0.0 ||
                rounded > static_cast<double>(std::numeric_limits<uword>::max()) ||
                std::fabs(pop(i) - rounded) > 1e-8) {
                Rcpp::stop("redist_constr support requires integer-valued populations.");
            }
            constraint_population(i) = static_cast<uword>(rounded);
        }
    }

    std::vector<BaseEdge> base_edges;
    std::vector<std::vector<BaseEdge>> edges_by_lower_vertex(n_vertices);
    for (int u = 0; u < n_vertices; u++) {
        for (int v : graph[u]) {
            if (v > u) {
                BaseEdge edge = {u, v, hierarchy_edge_bias(u, v, hierarchy, cut_bias)};
                base_edges.push_back(edge);
                edges_by_lower_vertex[u].push_back(edge);
            }
        }
    }

    uvec plan = init;
    double constraint_penalty = evaluate_constraint_penalty(
        plan, n_districts, target, constraint_population, graph, constraints);
    std::vector<HDimerForestEdge> forest;
    HDimerWorkspace workspace(n_vertices, n_districts);
    Clock::time_point forest_start = Clock::now();
    if (!sample_hierarchical_forest(graph, plan, n_districts, hierarchy, forest)) {
        Rcpp::stop("Initial plan is outside the hierarchical spanning-tree support.");
    }
    double forest_seconds = elapsed_seconds(forest_start);
    int forest_draws = 1;
    int tree_district_draws = n_districts;

    int n_saved = N / thin;
    umat plans(n_vertices, n_saved + 1, fill::zeros);
    plans.col(0) = plan;
    IntegerVector changed(n_saved + 1);

    int attempts = 0;
    int changed_count = 0;
    int identity_count = 0;
    int nonplanar_count = 0;
    int numerical_failure_count = 0;
    int no_matching_count = 0;
    int constraint_proposals = 0;
    int constraint_rejection_count = 0;
    int subset_dp_count = 0;
    double changed_edges_sum = 0.0;
    double changed_vertices_sum = 0.0;
    double old_probability_sum = 0.0;
    int old_probability_n = 0;
    double cut_seconds = 0.0;
    double fragment_seconds = 0.0;
    double candidate_seconds = 0.0;
    double matching_seconds = 0.0;
    double relabel_seconds = 0.0;

    RObject bar = R_NilValue;
    if (verbosity >= 1) {
        Rcout << "HIERARCHICAL SPLIT-DIMER MCMC\n";
        Rcout << "Sampling " << N << " transitions on " << n_vertices << " units and "
              << n_districts << " districts.\n";
        bar = cli_progress_bar(N, cli_config(false));
    }

    int saved = 0;
    for (int iteration = 1; iteration <= N; iteration++) {
        if (!fill_district_vertices(plan, n_districts, workspace.district_vertices)) {
            Rcpp::stop("A sampled plan has invalid district labels.");
        }
        std::vector<int> selected_labels = sample_update_districts(n_districts, ell);
        std::vector<bool> selected_districts(n_districts, false);
        for (int label : selected_labels) {
            selected_districts[label] = true;
        }
        if (iteration > 1 && (iteration - 1) % refresh == 0) {
            forest_start = Clock::now();
            if (!sample_hierarchical_forest_subset(graph, plan, n_districts, hierarchy,
                                                   selected_districts,
                                                   workspace.district_vertices, forest)) {
                Rcpp::stop("A sampled plan left the hierarchical spanning-tree support.");
            }
            forest_seconds += elapsed_seconds(forest_start);
            forest_draws++;
            tree_district_draws += ell;
        }

        HDimerStepDiagnostics step =
            run_dimer_step(base_edges, edges_by_lower_vertex, pop, hierarchy, split_penalty,
                           cut_bias, boundary_penalty, lower, upper, n_districts, plan,
                           constraint_penalty, constraint_population, target, graph, constraints,
                           forest, selected_labels, workspace);
        attempts++;
        if (step.changed) {
            changed_count++;
        } else {
            identity_count++;
        }
        if (!step.planar) {
            nonplanar_count++;
        }
        if (step.numerical_failure) {
            numerical_failure_count++;
        }
        if (step.no_matching) {
            no_matching_count++;
        }
        if (constraints.size() > 0 && !step.no_matching && step.planar) {
            constraint_proposals++;
        }
        if (step.constraint_rejected) {
            constraint_rejection_count++;
        }
        if (step.used_subset_dp) {
            subset_dp_count++;
        }
        changed_edges_sum += step.matching_edges_changed;
        changed_vertices_sum += step.vertices_changed;
        cut_seconds += step.cut_seconds;
        fragment_seconds += step.fragment_seconds;
        candidate_seconds += step.candidate_seconds;
        matching_seconds += step.matching_seconds;
        relabel_seconds += step.relabel_seconds;
        if (std::isfinite(step.old_matching_probability)) {
            old_probability_sum += step.old_matching_probability;
            old_probability_n++;
        }

        if (iteration % thin == 0) {
            changed[saved + 1] = arma::any(plan != plans.col(saved)) ? 1 : 0;
            plans.col(saved + 1) = plan;
            saved++;
        }
        if (verbosity >= 1 && CLI_SHOULD_TICK) {
            cli_progress_set(bar, iteration);
        }
        if (iteration % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }
    if (verbosity >= 1) {
        cli_progress_done(bar);
    }

    double old_probability_mean = NA_REAL;
    if (old_probability_n > 0) {
        old_probability_mean = old_probability_sum / old_probability_n;
    }
    double constraint_acceptance = NA_REAL;
    if (constraints.size() == 0) {
        constraint_acceptance = 1.0;
    } else if (constraint_proposals > 0) {
        constraint_acceptance =
            static_cast<double>(constraint_proposals - constraint_rejection_count) /
            constraint_proposals;
    }
    List diagnostics = List::create(
        Named("proposal_attempts") = attempts, Named("proposal_changed") = changed_count,
        Named("proposal_identity") = identity_count,
        Named("proposal_nonplanar") = nonplanar_count,
        Named("proposal_numerical_failure") = numerical_failure_count,
        Named("proposal_no_matching") = no_matching_count,
        Named("constraint_proposals") = constraint_proposals,
        Named("constraint_rejections") = constraint_rejection_count,
        Named("constraint_acceptance") = constraint_acceptance,
        Named("proposal_subset_dp") = subset_dp_count,
        Named("mean_vertices_changed") = changed_vertices_sum / attempts,
        Named("mean_matching_edges_changed") = changed_edges_sum / attempts,
        Named("mean_old_matching_prob") = old_probability_mean, Named("l") = ell,
        Named("tree_draws") = forest_draws, Named("tree_district_draws") = tree_district_draws,
        Named("seconds_tree_sampling") = forest_seconds, Named("seconds_cut") = cut_seconds,
        Named("seconds_fragments") = fragment_seconds,
        Named("seconds_candidates") = candidate_seconds,
        Named("seconds_matching") = matching_seconds,
        Named("seconds_relabel") = relabel_seconds);

    return List::create(Named("plans") = plans, Named("transition_changed") = changed,
                        Named("diagnostics") = diagnostics,
                        Named("algorithm") = "hierarchical_split_dimer");
}
