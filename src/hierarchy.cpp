#include "hierarchy.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace {

class DisjointSet {
  public:
    explicit DisjointSet(int size) : parent_(size), rank_(size, 0) {
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    bool unite(int first, int second) {
        int first_root = find(first);
        int second_root = find(second);
        if (first_root == second_root) return false;

        if (rank_[first_root] < rank_[second_root]) {
            std::swap(first_root, second_root);
        }
        parent_[second_root] = first_root;
        if (rank_[first_root] == rank_[second_root]) rank_[first_root]++;
        return true;
    }

  private:
    int find(int value) {
        if (parent_[value] != value) parent_[value] = find(parent_[value]);
        return parent_[value];
    }

    std::vector<int> parent_;
    std::vector<int> rank_;
};

bool same_unit(const uvec &level, int first, int second) {
    return level(first) == level(second);
}

std::vector<int> active_vertices_for_plan(const uvec &plan) {
    std::vector<int> vertices;
    vertices.reserve(plan.n_elem);
    for (uword i = 0; i < plan.n_elem; i++) {
        if (plan(i) > 0) vertices.push_back(i);
    }
    return vertices;
}

bool connected_by_label(const Graph &g,
                        const std::vector<int> &vertices,
                        const std::vector<int> &labels,
                        const uvec &level) {
    if (vertices.empty()) return true;

    std::vector<bool> active(g.size(), false);
    std::vector<bool> visited(g.size(), false);
    for (int vertex : vertices) active[vertex] = true;

    std::vector<int> stack;
    stack.reserve(vertices.size());
    std::set<std::pair<int, uword>> completed;
    for (int start : vertices) {
        if (visited[start]) continue;

        int label = labels[start];
        uword unit = level(start);
        if (!completed.insert({label, unit}).second) return false;
        visited[start] = true;
        stack.push_back(start);
        while (!stack.empty()) {
            int vertex = stack.back();
            stack.pop_back();
            for (int neighbor : g[vertex]) {
                if (!active[neighbor] || visited[neighbor]) continue;
                if (labels[neighbor] != label || level(neighbor) != unit) continue;
                visited[neighbor] = true;
                stack.push_back(neighbor);
            }
        }
    }
    return true;
}

bool level_is_forest(const Graph &g,
                     const uvec &plan,
                     const uvec &level) {
    int n_districts = 0;
    for (uword i = 0; i < plan.n_elem; i++) {
        n_districts = std::max(n_districts, static_cast<int>(plan(i)));
    }
    if (n_districts <= 1) return true;

    std::map<uword, std::set<int>> districts_by_unit;
    for (uword i = 0; i < plan.n_elem; i++) {
        if (plan(i) > 0) districts_by_unit[level(i)].insert(plan(i));
    }

    std::map<uword, std::vector<std::pair<int, int>>> edges_by_unit;
    for (int vertex = 0; vertex < static_cast<int>(g.size()); vertex++) {
        if (plan(vertex) == 0) continue;
        for (int neighbor : g[vertex]) {
            if (neighbor <= vertex || plan(neighbor) == 0) continue;
            if (!same_unit(level, vertex, neighbor)) continue;
            int first = plan(vertex);
            int second = plan(neighbor);
            if (first == second) continue;
            if (first > second) std::swap(first, second);
            edges_by_unit[level(vertex)].push_back({first, second});
        }
    }

    int split_count = 0;
    int edge_count = 0;
    DisjointSet forest(n_districts + 1);
    for (const auto &[unit, districts] : districts_by_unit) {
        split_count += static_cast<int>(districts.size()) - 1;
        auto edge_it = edges_by_unit.find(unit);
        if (edge_it == edges_by_unit.end()) continue;
        for (const auto &[first, second] : edge_it->second) {
            edge_count++;
            if (!forest.unite(first, second)) return false;
        }
    }

    return split_count <= n_districts - 1 && edge_count == split_count;
}

double log_tree_count(const Graph &g, const std::vector<int> &vertices) {
    int size = static_cast<int>(vertices.size());
    if (size <= 1) return 0;

    std::vector<int> positions(g.size(), -1);
    for (int i = 0; i < size; i++) positions[vertices[i]] = i;

    mat laplacian = zeros<mat>(size - 1, size - 1);
    for (int i = 1; i < size; i++) {
        int vertex = vertices[i];
        int position = i - 1;
        for (int neighbor : g[vertex]) {
            int neighbor_position = positions[neighbor];
            if (neighbor_position < 0) continue;
            laplacian(position, position) += 1;
            if (neighbor_position > 0) {
                laplacian(position, neighbor_position - 1) -= 1;
            }
        }
    }

    double log_determinant = 0;
    double sign = 0;
    arma::log_det(log_determinant, sign, laplacian);
    if (sign <= 0 || !std::isfinite(log_determinant)) {
        return -std::numeric_limits<double>::infinity();
    }
    return log_determinant;
}

double log_quotient_tree_count(const Graph &g,
                               const uvec &plan,
                               int district,
                               const uvec &parent_level,
                               uword parent_unit,
                               const uvec &child_level,
                               bool use_parent) {
    std::map<uword, int> child_positions;
    std::vector<int> vertices;
    for (uword i = 0; i < plan.n_elem; i++) {
        if (plan(i) != static_cast<uword>(district)) continue;
        if (use_parent && parent_level(i) != parent_unit) {
            continue;
        }
        vertices.push_back(i);
        if (child_positions.find(child_level(i)) == child_positions.end()) {
            int position = child_positions.size();
            child_positions.emplace(child_level(i), position);
        }
    }

    int size = static_cast<int>(child_positions.size());
    if (size <= 1) return 0;

    mat laplacian = zeros<mat>(size - 1, size - 1);
    for (int vertex : vertices) {
        int position = child_positions.at(child_level(vertex));
        for (int neighbor : g[vertex]) {
            if (plan(neighbor) != static_cast<uword>(district)) continue;
            if (use_parent && parent_level(neighbor) != parent_unit) {
                continue;
            }
            int neighbor_position = child_positions.at(child_level(neighbor));
            if (neighbor_position == position) continue;
            if (position > 0) {
                laplacian(position - 1, position - 1) += 1;
                if (neighbor_position > 0) {
                    laplacian(position - 1, neighbor_position - 1) -= 1;
                }
            }
        }
    }

    double log_determinant = 0;
    double sign = 0;
    arma::log_det(log_determinant, sign, laplacian);
    if (sign <= 0 || !std::isfinite(log_determinant)) {
        return -std::numeric_limits<double>::infinity();
    }
    return log_determinant;
}

} // namespace

HierarchySpec hierarchy_spec_from_control(const List &control,
                                          const uvec &fallback_counties) {
    HierarchySpec hierarchy;
    bool has_hierarchy_mode = false;
    if (control.containsElementNamed("hierarchy_mode")) {
        has_hierarchy_mode = true;
        int mode = Rcpp::as<int>(control["hierarchy_mode"]);
        if (mode == 1) hierarchy.mode = HierarchyMode::speedup;
        if (mode == 2) hierarchy.mode = HierarchyMode::strict;
    }

    bool has_multiple_fallback_units = false;
    if (fallback_counties.n_elem > 0) {
        std::set<uword> fallback_units(
            fallback_counties.begin(), fallback_counties.end()
        );
        has_multiple_fallback_units = fallback_units.size() > 1;
    }
    if (!has_hierarchy_mode && has_multiple_fallback_units) {
        hierarchy.mode = HierarchyMode::speedup;
    }

    if (control.containsElementNamed("hierarchy_levels")) {
        IntegerMatrix level_matrix = control["hierarchy_levels"];
        int V = level_matrix.nrow();
        for (int j = 0; j < level_matrix.ncol(); j++) {
            uvec level(V);
            for (int i = 0; i < V; i++) level(i) = level_matrix(i, j);
            hierarchy.levels.push_back(level);
        }
    } else if (hierarchy.mode != HierarchyMode::none && fallback_counties.n_elem > 0) {
        hierarchy.levels.push_back(fallback_counties);
    }

    for (int j = static_cast<int>(hierarchy.levels.size()) - 1; j >= 0; j--) {
        hierarchy.sampler_levels.push_back(hierarchy.levels[j]);
    }
    return hierarchy;
}

bool hierarchy_region_connected(const Graph &g,
                                const std::vector<int> &region_vertices,
                                const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    std::vector<int> labels(g.size(), 1);
    for (const uvec &level : hierarchy.levels) {
        if (!connected_by_label(g, region_vertices, labels, level)) return false;
    }
    return true;
}

bool hierarchy_plan_connected(const Graph &g,
                              const uvec &plan,
                              const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    std::vector<int> vertices = active_vertices_for_plan(plan);
    std::vector<int> labels(plan.n_elem, 0);
    for (uword i = 0; i < plan.n_elem; i++) labels[i] = plan(i);
    for (const uvec &level : hierarchy.levels) {
        if (!connected_by_label(g, vertices, labels, level)) return false;
    }
    return true;
}

bool hierarchy_plan_valid(const Graph &g,
                          const uvec &plan,
                          const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    if (!hierarchy_plan_connected(g, plan, hierarchy)) return false;
    for (const uvec &level : hierarchy.levels) {
        if (!level_is_forest(g, plan, level)) return false;
    }
    return true;
}

int admissible_boundary_count(const Graph &g,
                              const uvec &plan,
                              int first_label,
                              const std::vector<int> &other_labels,
                              const std::vector<int> &region_vertices,
                              const HierarchySpec &hierarchy) {
    if (region_vertices.empty()) return 0;

    std::vector<bool> is_other(plan.n_elem, false);
    for (int label : other_labels) {
        for (uword i = 0; i < plan.n_elem; i++) {
            if (static_cast<int>(plan(i)) == label) is_other[i] = true;
        }
    }

    std::vector<bool> in_region(plan.n_elem, false);
    for (int vertex : region_vertices) in_region[vertex] = true;

    int deepest_shared_level = -1;
    uword deepest_shared_unit = 0;
    for (int level_index = 0;
         level_index < static_cast<int>(hierarchy.levels.size());
         level_index++) {
        std::set<uword> first_units;
        std::set<uword> other_units;
        for (int vertex : region_vertices) {
            if (static_cast<int>(plan(vertex)) == first_label) {
                first_units.insert(hierarchy.levels[level_index](vertex));
            } else if (is_other[vertex]) {
                other_units.insert(hierarchy.levels[level_index](vertex));
            }
        }

        std::vector<uword> shared_units;
        std::set_intersection(first_units.begin(), first_units.end(),
                              other_units.begin(), other_units.end(),
                              std::back_inserter(shared_units));
        if (shared_units.size() > 1) return 0;
        if (shared_units.size() == 1) {
            deepest_shared_level = level_index;
            deepest_shared_unit = shared_units[0];
        }
    }

    int count = 0;
    for (int vertex : region_vertices) {
        if (static_cast<int>(plan(vertex)) != first_label) continue;
        for (int neighbor : g[vertex]) {
            if (!in_region[neighbor] || !is_other[neighbor]) continue;
            if (deepest_shared_level >= 0 &&
                (hierarchy.levels[deepest_shared_level](vertex) !=
                 deepest_shared_unit ||
                 hierarchy.levels[deepest_shared_level](neighbor) !=
                 deepest_shared_unit)) {
                continue;
            }
            count++;
        }
    }
    return count;
}

double log_hierarchical_tree_count(const Graph &g,
                                   const uvec &plan,
                                   int district,
                                   const HierarchySpec &hierarchy) {
    std::vector<int> district_vertices;
    for (uword i = 0; i < plan.n_elem; i++) {
        if (static_cast<int>(plan(i)) == district) district_vertices.push_back(i);
    }

    if (!hierarchy.enabled()) return log_tree_count(g, district_vertices);

    double result = 0;
    int n_levels = static_cast<int>(hierarchy.levels.size());
    const uvec &finest_level = hierarchy.levels.back();

    std::map<uword, std::vector<int>> finest_vertices;
    for (int vertex : district_vertices) {
        finest_vertices[finest_level(vertex)].push_back(vertex);
    }
    for (const auto &[unit, vertices] : finest_vertices) {
        (void) unit;
        result += log_tree_count(g, vertices);
    }

    for (int child_level_index = 0;
         child_level_index < n_levels;
         child_level_index++) {
        const uvec &child_level = hierarchy.levels[child_level_index];
        if (child_level_index == 0) {
            result += log_quotient_tree_count(g, plan, district,
                                              child_level, 0,
                                              child_level, false);
            continue;
        }

        const uvec &parent_level = hierarchy.levels[child_level_index - 1];
        std::set<uword> parent_units;
        for (int vertex : district_vertices) {
            parent_units.insert(parent_level(vertex));
        }
        for (uword parent_unit : parent_units) {
            result += log_quotient_tree_count(
                g, plan, district, parent_level,
                parent_unit, child_level, true);
        }
    }
    return result;
}
