#include "hierarchy.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <numeric>
#include <set>
#include <utility>

namespace {

constexpr int max_precomputed_unit_size = 12;

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

std::vector<int> active_vertices_for_plan(const uvec &plan) {
    std::vector<int> vertices;
    vertices.reserve(plan.n_elem);
    for (uword i = 0; i < plan.n_elem; i++) {
        if (plan(i) > 0) vertices.push_back(i);
    }
    return vertices;
}

bool connected_region_by_cache(const Graph &g,
                               const std::vector<int> &region_vertices,
                               const std::vector<char> &region_mark,
                               const HierarchyLevelCache &cache) {
    std::vector<char> visited(g.size(), false);
    std::vector<int> stack;
    stack.reserve(region_vertices.size());

    for (int unit_index = 0;
         unit_index < static_cast<int>(cache.units.size()); unit_index++) {
        const HierarchyUnitCache &unit = cache.units[unit_index];
        if (!unit.connected_masks.empty()) {
            std::uint64_t mask = 0;
            for (int vertex_index = 0;
                 vertex_index < static_cast<int>(unit.vertices.size());
                 vertex_index++) {
                if (region_mark[unit.vertices[vertex_index]]) {
                    mask |= std::uint64_t{1} << vertex_index;
                }
            }
            if (mask != 0 &&
                (unit.connected_masks[mask / 64] &
                 (std::uint64_t{1} << (mask % 64))) == 0) {
                return false;
            }
            continue;
        }

        bool found_component = false;
        for (int start : unit.vertices) {
            if (!region_mark[start] || visited[start]) continue;
            if (found_component) return false;
            found_component = true;
            visited[start] = true;
            stack.push_back(start);
            while (!stack.empty()) {
                int vertex = stack.back();
                stack.pop_back();
                for (int neighbor : g[vertex]) {
                    if (!region_mark[neighbor] || visited[neighbor] ||
                        cache.labels[neighbor] != unit_index + 1) {
                        continue;
                    }
                    visited[neighbor] = true;
                    stack.push_back(neighbor);
                }
            }
        }
    }
    return true;
}

template <typename Labels>
bool connected_by_label(const Graph &g,
                        const std::vector<int> &vertices,
                        const Labels &labels,
                        const HierarchyLevelCache &cache) {
    if (vertices.empty()) return true;

    std::vector<char> active(g.size(), false);
    for (int vertex : vertices) active[vertex] = true;

    std::vector<char> visited(g.size(), false);
    std::vector<int> stack;
    stack.reserve(vertices.size());

    for (int unit_index = 0;
         unit_index < static_cast<int>(cache.units.size()); unit_index++) {
        const HierarchyUnitCache &unit = cache.units[unit_index];
        if (!unit.connected_masks.empty()) {
            std::array<int, max_precomputed_unit_size> label_values;
            std::array<std::uint64_t, max_precomputed_unit_size> label_masks{};
            int n_labels = 0;
            for (int vertex_index = 0;
                 vertex_index < static_cast<int>(unit.vertices.size());
                 vertex_index++) {
                int vertex = unit.vertices[vertex_index];
                if (!active[vertex]) continue;
                int label = labels[vertex];
                int label_index = 0;
                while (label_index < n_labels &&
                       label_values[label_index] != label) {
                    label_index++;
                }
                if (label_index == n_labels) {
                    label_values[n_labels] = label;
                    n_labels++;
                }
                label_masks[label_index] |= std::uint64_t{1} << vertex_index;
            }

            for (int label_index = 0; label_index < n_labels; label_index++) {
                std::uint64_t mask = label_masks[label_index];
                if ((unit.connected_masks[mask / 64] &
                     (std::uint64_t{1} << (mask % 64))) == 0) {
                    return false;
                }
            }
            continue;
        }

        std::vector<int> completed_labels;
        for (int start : unit.vertices) {
            if (!active[start] || visited[start]) continue;

            int label = labels[start];
            if (std::find(completed_labels.begin(), completed_labels.end(),
                          label) != completed_labels.end()) {
                return false;
            }
            completed_labels.push_back(label);
            visited[start] = true;
            stack.push_back(start);
            while (!stack.empty()) {
                int vertex = stack.back();
                stack.pop_back();
                for (int neighbor : g[vertex]) {
                    if (!active[neighbor] || visited[neighbor]) continue;
                    if (cache.labels[neighbor] != unit_index + 1 ||
                        labels[neighbor] != label) continue;
                    visited[neighbor] = true;
                    stack.push_back(neighbor);
                }
            }
        }
    }
    return true;
}

bool level_is_forest(const Graph &g,
                     const uvec &plan,
                     const HierarchyLevelCache &cache) {
    int n_districts = 0;
    for (uword i = 0; i < plan.n_elem; i++) {
        n_districts = std::max(n_districts, static_cast<int>(plan(i)));
    }
    if (n_districts <= 1) return true;

    std::vector<std::set<int>> districts_by_unit(cache.units.size());
    for (uword i = 0; i < plan.n_elem; i++) {
        if (plan(i) > 0) {
            districts_by_unit[cache.labels[i] - 1].insert(plan(i));
        }
    }

    std::vector<std::vector<std::pair<int, int>>> edges_by_unit(
        cache.units.size());
    for (int vertex = 0; vertex < static_cast<int>(g.size()); vertex++) {
        if (plan(vertex) == 0) continue;
        for (int neighbor : g[vertex]) {
            if (neighbor <= vertex || plan(neighbor) == 0) continue;
            if (cache.labels[vertex] != cache.labels[neighbor]) continue;
            int first = plan(vertex);
            int second = plan(neighbor);
            if (first == second) continue;
            if (first > second) std::swap(first, second);
            edges_by_unit[cache.labels[vertex] - 1].push_back({first, second});
        }
    }

    int split_count = 0;
    int edge_count = 0;
    DisjointSet forest(n_districts + 1);
    for (int unit = 0;
         unit < static_cast<int>(districts_by_unit.size()); unit++) {
        if (districts_by_unit[unit].empty()) continue;
        split_count += static_cast<int>(districts_by_unit[unit].size()) - 1;
        for (const auto &[first, second] : edges_by_unit[unit]) {
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

HierarchyLevelCache build_level_cache(const Graph &g, const uvec &level) {
    HierarchyLevelCache cache;
    std::map<uword, int> relabel;
    cache.labels.resize(level.n_elem);
    for (uword vertex = 0; vertex < level.n_elem; vertex++) {
        auto label = relabel.emplace(level(vertex), relabel.size() + 1);
        cache.labels[vertex] = label.first->second;
    }

    cache.units.resize(relabel.size());
    cache.local_positions.assign(level.n_elem, -1);
    for (uword vertex = 0; vertex < level.n_elem; vertex++) {
        int unit = cache.labels[vertex] - 1;
        cache.local_positions[vertex] = cache.units[unit].vertices.size();
        cache.units[unit].vertices.push_back(vertex);
    }

    for (int unit_index = 0;
         unit_index < static_cast<int>(cache.units.size()); unit_index++) {
        HierarchyUnitCache &unit = cache.units[unit_index];
        int size = unit.vertices.size();
        if (size > max_precomputed_unit_size) continue;

        std::vector<std::uint64_t> neighbors(size, 0);
        for (int vertex_index = 0; vertex_index < size; vertex_index++) {
            int vertex = unit.vertices[vertex_index];
            for (int neighbor : g[vertex]) {
                if (cache.labels[neighbor] != unit_index + 1) continue;
                int neighbor_index = cache.local_positions[neighbor];
                neighbors[vertex_index] |= (std::uint64_t{1} << neighbor_index);
            }
        }

        int n_masks = 1 << size;
        unit.connected_masks.assign((n_masks + 63) / 64, 0);
        for (int mask = 1; mask < n_masks; mask++) {
            int first = 0;
            while ((mask & (1 << first)) == 0) first++;

            std::uint64_t seen = std::uint64_t{1} << first;
            std::uint64_t frontier = seen;
            while (frontier != 0) {
                int vertex = 0;
                while ((frontier & (std::uint64_t{1} << vertex)) == 0) {
                    vertex++;
                }
                frontier &= ~(std::uint64_t{1} << vertex);
                std::uint64_t added = neighbors[vertex] & mask & ~seen;
                seen |= added;
                frontier |= added;
            }

            if (seen == static_cast<std::uint64_t>(mask)) {
                unit.connected_masks[mask / 64] |=
                    std::uint64_t{1} << (mask % 64);
            }
        }
    }
    return cache;
}

} // namespace

HierarchySpec hierarchy_spec_from_control(const List &control,
                                          const uvec &fallback_counties,
                                          const Graph &g) {
    HierarchySpec hierarchy;
    if (control.containsElementNamed("hierarchy_mode")) {
        int mode = Rcpp::as<int>(control["hierarchy_mode"]);
        if (mode == 1) hierarchy.mode = HierarchyMode::connected;
        if (mode == 2) hierarchy.mode = HierarchyMode::strict;
        if (mode == 3) hierarchy.mode = HierarchyMode::heuristic;
        if (mode < 0 || mode > 3) {
            Rcpp::stop("Invalid hierarchy mode.");
        }
    }

    if (control.containsElementNamed("hierarchy_levels")) {
        IntegerMatrix level_matrix = control["hierarchy_levels"];
        int V = level_matrix.nrow();
        for (int j = 0; j < level_matrix.ncol(); j++) {
            uvec level(V);
            for (int i = 0; i < V; i++) level(i) = level_matrix(i, j);
            hierarchy.levels.push_back(level);
        }
    } else if (hierarchy.mode != HierarchyMode::none &&
               fallback_counties.n_elem > 0) {
        hierarchy.levels.push_back(fallback_counties);
    }

    if (hierarchy.enabled()) {
        for (const uvec &level : hierarchy.levels) {
            hierarchy.level_caches.push_back(build_level_cache(g, level));
        }
    }

    for (int j = static_cast<int>(hierarchy.levels.size()) - 1; j >= 0; j--) {
        hierarchy.sampler_levels.push_back(hierarchy.levels[j]);
    }

    int n_levels = hierarchy.sampler_levels.size();
    hierarchy.sampler_labels.resize(n_levels);
    hierarchy.sampler_group_counts.resize(n_levels);
    hierarchy.sampler_parents.resize(n_levels);
    for (int level_index = 0; level_index < n_levels; level_index++) {
        const uvec &level = hierarchy.sampler_levels[level_index];
        std::map<uword, int> relabel;
        std::vector<int> labels(level.n_elem);
        for (uword vertex = 0; vertex < level.n_elem; vertex++) {
            auto label = relabel.emplace(level(vertex), relabel.size() + 1);
            labels[vertex] = label.first->second;
        }
        hierarchy.sampler_labels[level_index] = std::move(labels);
        hierarchy.sampler_group_counts[level_index] = relabel.size();
    }

    for (int level_index = 0; level_index < n_levels; level_index++) {
        int n_groups = hierarchy.sampler_group_counts[level_index];
        hierarchy.sampler_parents[level_index].assign(n_groups, -1);
        if (level_index == n_levels - 1) continue;

        const std::vector<int> &child_labels =
            hierarchy.sampler_labels[level_index];
        const std::vector<int> &parent_labels =
            hierarchy.sampler_labels[level_index + 1];
        for (uword vertex = 0; vertex < child_labels.size(); vertex++) {
            int child = child_labels[vertex] - 1;
            hierarchy.sampler_parents[level_index][child] =
                parent_labels[vertex] - 1;
        }
    }

    if (n_levels > 0) {
        const std::vector<int> &finest_labels = hierarchy.sampler_labels[0];
        int V = finest_labels.size();
        hierarchy.sampler_finest_off.resize(V + 1);
        for (int vertex = 0; vertex < V; vertex++) {
            hierarchy.sampler_finest_off[vertex] =
                hierarchy.sampler_finest_adj.size();
            for (int neighbor : g[vertex]) {
                if (finest_labels[neighbor] == finest_labels[vertex]) {
                    hierarchy.sampler_finest_adj.push_back(neighbor);
                }
            }
        }
        hierarchy.sampler_finest_off[V] =
            hierarchy.sampler_finest_adj.size();
    }

    return hierarchy;
}

bool hierarchy_region_connected(const Graph &g,
                                const std::vector<int> &region_vertices,
                                const std::vector<char> &region_mark,
                                const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    for (const HierarchyLevelCache &cache : hierarchy.level_caches) {
        if (!connected_region_by_cache(g, region_vertices, region_mark, cache)) {
            return false;
        }
    }
    return true;
}

bool hierarchy_plan_connected(const Graph &g,
                              const uvec &plan,
                              const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    std::vector<int> vertices = active_vertices_for_plan(plan);
    for (const HierarchyLevelCache &cache : hierarchy.level_caches) {
        if (!connected_by_label(g, vertices, plan, cache)) {
            return false;
        }
    }
    return true;
}

bool hierarchy_plan_valid(const Graph &g,
                          const uvec &plan,
                          const HierarchySpec &hierarchy) {
    if (!hierarchy.enabled()) return true;
    if (!hierarchy_plan_connected(g, plan, hierarchy)) return false;
    for (const HierarchyLevelCache &cache : hierarchy.level_caches) {
        if (!level_is_forest(g, plan, cache)) return false;
    }
    return true;
}

int admissible_boundary_count(const Graph &g,
                              const uvec &plan,
                              int first_label,
                              const std::vector<int> &other_labels,
                              const std::vector<int> &region_vertices,
                              const std::vector<char> &region_mark,
                              const HierarchySpec &hierarchy) {
    if (region_vertices.empty()) return 0;

    if (!hierarchy.enabled()) {
        int count = 0;
        for (int vertex : region_vertices) {
            if (static_cast<int>(plan(vertex)) != first_label) continue;
            for (int neighbor : g[vertex]) {
                if (!region_mark[neighbor]) continue;
                if (std::find(other_labels.begin(), other_labels.end(),
                              static_cast<int>(plan(neighbor))) ==
                    other_labels.end()) continue;
                count++;
            }
        }
        return count;
    }

    int deepest_shared_level = -1;
    int deepest_shared_unit = 0;
    std::vector<int> first_units;
    std::vector<int> other_units;
    first_units.reserve(region_vertices.size());
    other_units.reserve(region_vertices.size());
    for (int level_index = 0;
         level_index < static_cast<int>(hierarchy.levels.size());
         level_index++) {
        first_units.clear();
        other_units.clear();
        const std::vector<int> &labels =
            hierarchy.level_caches[level_index].labels;
        for (int vertex : region_vertices) {
            if (static_cast<int>(plan(vertex)) == first_label) {
                first_units.push_back(labels[vertex]);
            } else if (std::find(other_labels.begin(), other_labels.end(),
                                 static_cast<int>(plan(vertex))) !=
                       other_labels.end()) {
                other_units.push_back(labels[vertex]);
            }
        }

        std::sort(first_units.begin(), first_units.end());
        first_units.erase(std::unique(first_units.begin(), first_units.end()),
                          first_units.end());
        std::sort(other_units.begin(), other_units.end());
        other_units.erase(std::unique(other_units.begin(), other_units.end()),
                          other_units.end());
        int shared_count = 0;
        int shared_unit = 0;
        for (int unit : first_units) {
            if (!std::binary_search(other_units.begin(), other_units.end(),
                                    unit)) {
                continue;
            }
            shared_count++;
            shared_unit = unit;
            if (shared_count > 1) return 0;
        }
        if (shared_count == 1) {
            deepest_shared_level = level_index;
            deepest_shared_unit = shared_unit;
        }
    }

    int count = 0;
    for (int vertex : region_vertices) {
        if (static_cast<int>(plan(vertex)) != first_label) continue;
        for (int neighbor : g[vertex]) {
            if (!region_mark[neighbor]) continue;
            if (std::find(other_labels.begin(), other_labels.end(),
                          static_cast<int>(plan(neighbor))) ==
                other_labels.end()) continue;
            if (deepest_shared_level >= 0 &&
                (hierarchy.level_caches[deepest_shared_level].labels[vertex] !=
                 deepest_shared_unit ||
                 hierarchy.level_caches[deepest_shared_level].labels[neighbor] !=
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
