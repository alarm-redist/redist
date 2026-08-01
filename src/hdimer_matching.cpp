#include "hdimer_matching.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <tuple>

namespace {
using PlanarGraph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
                          boost::property<boost::edge_index_t, int>>;

constexpr int SUBSET_DP_MAX_VERTICES = 14;
constexpr double PROBABILITY_SUM_TOLERANCE = 1e-8;
constexpr double NEGATIVE_PROBABILITY_TOLERANCE = 1e-10;

struct AggregatedEdge {
    int u;
    int v;
    double log_weight;
    std::vector<HDimerMatchingEdge> originals;
};

struct FaceTraversalStep {
    int edge;
    int start;
};

struct FaceCollector : boost::planar_face_traversal_visitor {
    explicit FaceCollector(const PlanarGraph &graph) : graph(graph) {}

    void begin_face() { faces.emplace_back(); }

    void next_vertex(PlanarGraph::vertex_descriptor vertex) {
        last_vertex = static_cast<int>(vertex);
    }

    void next_edge(PlanarGraph::edge_descriptor edge) {
        faces.back().push_back({get(boost::edge_index, graph, edge), last_vertex});
    }

    const PlanarGraph &graph;
    int last_vertex = -1;
    std::vector<std::vector<FaceTraversalStep>> faces;
};

enum class OrientationStatus { success, nonplanar, numerical_failure };

double log_sum_exp(double left, double right) {
    double negative_infinity = -std::numeric_limits<double>::infinity();
    if (left == negative_infinity) {
        return right;
    }
    if (right == negative_infinity) {
        return left;
    }
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (left < right) {
        std::swap(left, right);
    }
    return left + std::log1p(std::exp(right - left));
}

HDimerMatchingResult matching_failure(HDimerMatchingStatus status,
                                      bool planarity_checked = false, bool planar = false,
                                      bool used_subset_dp = false) {
    HDimerMatchingResult result;
    result.status = status;
    result.planarity_checked = planarity_checked;
    result.planar = planar;
    result.used_subset_dp = used_subset_dp;
    return result;
}

const char *matching_status_label(HDimerMatchingStatus status) {
    switch (status) {
    case HDimerMatchingStatus::success:
        return "success";
    case HDimerMatchingStatus::invalid_input:
        return "invalid_input";
    case HDimerMatchingStatus::no_perfect_matching:
        return "no_perfect_matching";
    case HDimerMatchingStatus::nonplanar:
        return "nonplanar";
    case HDimerMatchingStatus::numerical_failure:
        return "numerical_failure";
    }
    return "unknown";
}

bool sample_parallel_payload(const AggregatedEdge &edge, int &payload) {
    double uniform = R::runif(0, 1);
    double cumulative = 0.0;
    bool selected = false;
    payload = edge.originals.back().payload;
    for (const HDimerMatchingEdge &original : edge.originals) {
        double probability = std::exp(original.log_weight - edge.log_weight);
        if (!std::isfinite(probability) || !(probability > 0.0)) {
            return false;
        }
        cumulative += probability;
        if (!std::isfinite(cumulative)) {
            return false;
        }
        if (!selected && uniform <= cumulative) {
            payload = original.payload;
            selected = true;
        }
    }
    return std::abs(cumulative - 1.0) <= PROBABILITY_SUM_TOLERANCE;
}

int lowest_set_vertex(std::size_t mask) { return static_cast<int>(__builtin_ctzll(mask)); }

double subset_log_partition(std::size_t mask, int n, const std::vector<AggregatedEdge> &edges,
                            const std::vector<int> &pair_index,
                            const std::vector<std::size_t> &neighbor_mask,
                            std::vector<double> &memo) {
    if (!std::isnan(memo[mask])) {
        return memo[mask];
    }
    int first = lowest_set_vertex(mask);
    std::size_t remainder = mask ^ (std::size_t(1) << first);
    double value = -std::numeric_limits<double>::infinity();
    std::size_t available = remainder & neighbor_mask[first];
    while (available) {
        int second = lowest_set_vertex(available);
        std::size_t second_bit = std::size_t(1) << second;
        available ^= second_bit;
        int edge = pair_index[static_cast<std::size_t>(first) * n + second];
        double rest = subset_log_partition(remainder ^ second_bit, n, edges, pair_index,
                                           neighbor_mask, memo);
        if (std::isfinite(rest)) {
            value = log_sum_exp(value, edges[edge].log_weight + rest);
        }
    }
    memo[mask] = value;
    return value;
}

HDimerMatchingResult sample_subset_dp(int n, const std::vector<AggregatedEdge> &edges,
                                      const std::vector<int> &pair_index,
                                      bool planarity_checked, bool planar) {
    std::size_t n_states = std::size_t(1) << n;
    std::vector<std::size_t> neighbor_mask(n, 0);
    for (const AggregatedEdge &edge : edges) {
        neighbor_mask[edge.u] |= std::size_t(1) << edge.v;
        neighbor_mask[edge.v] |= std::size_t(1) << edge.u;
    }
    std::vector<double> memo(n_states, std::numeric_limits<double>::quiet_NaN());
    memo[0] = 0.0;
    std::size_t full_mask = n_states - 1;
    double total = subset_log_partition(full_mask, n, edges, pair_index, neighbor_mask, memo);
    if (!std::isfinite(total)) {
        HDimerMatchingStatus status = total == -std::numeric_limits<double>::infinity()
                                          ? HDimerMatchingStatus::no_perfect_matching
                                          : HDimerMatchingStatus::numerical_failure;
        return matching_failure(status, planarity_checked, planar, true);
    }

    std::vector<int> selected;
    selected.reserve(n / 2);
    std::size_t mask = full_mask;
    double max_probability_error = 0.0;
    double min_probability = 1.0;
    while (mask) {
        int first = lowest_set_vertex(mask);
        std::size_t remainder = mask ^ (std::size_t(1) << first);
        std::vector<int> available;
        std::vector<std::size_t> next_masks;
        std::vector<double> probabilities;
        double sum = 0.0;
        std::size_t choices = remainder & neighbor_mask[first];
        while (choices) {
            int second = lowest_set_vertex(choices);
            std::size_t second_bit = std::size_t(1) << second;
            choices ^= second_bit;
            int edge = pair_index[static_cast<std::size_t>(first) * n + second];
            std::size_t next_mask = remainder ^ second_bit;
            double rest = memo[next_mask];
            if (!std::isfinite(rest)) {
                continue;
            }
            double probability = std::exp(edges[edge].log_weight + rest - memo[mask]);
            if (!std::isfinite(probability) || probability <= 0.0) {
                return matching_failure(HDimerMatchingStatus::numerical_failure,
                                        planarity_checked, planar, true);
            }
            available.push_back(edge);
            next_masks.push_back(next_mask);
            probabilities.push_back(probability);
            sum += probability;
        }
        double error = std::abs(sum - 1.0);
        max_probability_error = std::max(max_probability_error, error);
        if (available.empty() || !std::isfinite(sum) || sum <= 0.0 ||
            error > PROBABILITY_SUM_TOLERANCE) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, planarity_checked,
                                    planar, true);
        }
        for (double &probability : probabilities) {
            probability /= sum;
            min_probability = std::min(min_probability, probability);
        }
        double uniform = R::runif(0, 1), cumulative = 0.0;
        int selected_index = static_cast<int>(available.size()) - 1;
        for (int index = 0; index < static_cast<int>(available.size()); index++) {
            cumulative += probabilities[index];
            if (uniform <= cumulative) {
                selected_index = index;
                break;
            }
        }
        const AggregatedEdge &edge = edges[available[selected_index]];
        int payload = -1;
        if (!sample_parallel_payload(edge, payload)) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, planarity_checked,
                                    planar, true);
        }
        selected.push_back(payload);
        mask = next_masks[selected_index];
    }

    HDimerMatchingResult result;
    result.status = HDimerMatchingStatus::success;
    result.planarity_checked = planarity_checked;
    result.planar = planar;
    result.used_subset_dp = true;
    result.selected_edge_indices = std::move(selected);
    result.log_partition = total;
    result.max_probability_error = max_probability_error;
    result.min_probability = min_probability;
    return result;
}
// Construct a clockwise-odd Kasteleyn orientation. A primal spanning tree and
// its complementary dual spanning tree make every non-root face independently
// orientable in leaf-to-root order. The root-face parity then follows because
// a component with a perfect matching has an even number of vertices.
OrientationStatus build_clockwise_odd_orientation(const std::vector<int> &component_vertices,
                                                  const std::vector<int> &component_edges,
                                                  const std::vector<AggregatedEdge> &all_edges,
                                                  std::vector<bool> &orientation) {
    if (component_edges.empty()) {
        return component_vertices.size() <= 1 ? OrientationStatus::success
                                              : OrientationStatus::numerical_failure;
    }

    std::map<int, int> local_vertex;
    for (std::size_t index = 0; index < component_vertices.size(); index++) {
        local_vertex[component_vertices[index]] = static_cast<int>(index);
    }

    PlanarGraph graph(component_vertices.size());
    std::vector<PlanarGraph::edge_descriptor> graph_edges(component_edges.size());
    for (std::size_t index = 0; index < component_edges.size(); index++) {
        const AggregatedEdge &edge = all_edges[component_edges[index]];
        auto inserted = add_edge(local_vertex[edge.u], local_vertex[edge.v], graph);
        if (!inserted.second) {
            return OrientationStatus::numerical_failure;
        }
        graph_edges[index] = inserted.first;
        put(boost::edge_index, graph, inserted.first, static_cast<int>(index));
    }

    std::vector<std::vector<PlanarGraph::edge_descriptor>> embedding(component_vertices.size());
    auto embedding_map =
        boost::make_iterator_property_map(embedding.begin(), get(boost::vertex_index, graph));
    if (!boost::boyer_myrvold_planarity_test(boost::boyer_myrvold_params::graph = graph,
                                             boost::boyer_myrvold_params::embedding =
                                                 embedding_map)) {
        return OrientationStatus::nonplanar;
    }

    FaceCollector face_collector(graph);
    boost::planar_face_traversal(graph, embedding_map, face_collector);
    if (face_collector.faces.empty()) {
        return OrientationStatus::numerical_failure;
    }

    std::vector<bool> primal_tree_edge(component_edges.size(), false);
    std::vector<bool> vertex_seen(component_vertices.size(), false);
    std::queue<int> queue;
    vertex_seen[0] = true;
    queue.push(0);
    while (!queue.empty()) {
        int vertex = queue.front();
        queue.pop();
        auto adjacent = adjacent_vertices(vertex, graph);
        for (auto iterator = adjacent.first; iterator != adjacent.second; ++iterator) {
            int neighbor = static_cast<int>(*iterator);
            if (vertex_seen[neighbor]) {
                continue;
            }
            vertex_seen[neighbor] = true;
            queue.push(neighbor);
            auto edge = boost::edge(vertex, neighbor, graph);
            if (!edge.second) {
                return OrientationStatus::numerical_failure;
            }
            primal_tree_edge[get(boost::edge_index, graph, edge.first)] = true;
        }
    }
    if (std::any_of(vertex_seen.begin(), vertex_seen.end(), [](bool seen) { return !seen; })) {
        return OrientationStatus::numerical_failure;
    }

    std::vector<int> first_face(component_edges.size(), -1);
    std::vector<int> second_face(component_edges.size(), -1);
    for (std::size_t face = 0; face < face_collector.faces.size(); face++) {
        for (const FaceTraversalStep &step : face_collector.faces[face]) {
            if (first_face[step.edge] < 0) {
                first_face[step.edge] = static_cast<int>(face);
            } else if (first_face[step.edge] != static_cast<int>(face)) {
                second_face[step.edge] = static_cast<int>(face);
            }
        }
    }

    std::vector<std::vector<std::pair<int, int>>> dual_adjacency(face_collector.faces.size());
    for (std::size_t edge = 0; edge < component_edges.size(); edge++) {
        if (primal_tree_edge[edge]) {
            continue;
        }
        if (first_face[edge] < 0 || second_face[edge] < 0) {
            return OrientationStatus::numerical_failure;
        }
        dual_adjacency[first_face[edge]].push_back({second_face[edge], static_cast<int>(edge)});
        dual_adjacency[second_face[edge]].push_back({first_face[edge], static_cast<int>(edge)});
    }

    std::vector<int> dual_parent(face_collector.faces.size(), -2);
    std::vector<int> dual_parent_edge(face_collector.faces.size(), -1);
    std::vector<int> face_order;
    dual_parent[0] = -1;
    queue.push(0);
    while (!queue.empty()) {
        int face = queue.front();
        queue.pop();
        face_order.push_back(face);
        for (const auto &[neighbor, edge] : dual_adjacency[face]) {
            if (dual_parent[neighbor] != -2) {
                continue;
            }
            dual_parent[neighbor] = face;
            dual_parent_edge[neighbor] = edge;
            queue.push(neighbor);
        }
    }
    if (face_order.size() != face_collector.faces.size()) {
        return OrientationStatus::numerical_failure;
    }

    // +1 means the edge follows Boost's stored source-to-target direction.
    std::vector<int> direction(component_edges.size(), 0);
    for (std::size_t edge = 0; edge < component_edges.size(); edge++) {
        if (primal_tree_edge[edge]) {
            direction[edge] = 1;
        }
    }
    for (auto iterator = face_order.rbegin(); iterator != face_order.rend(); ++iterator) {
        int face = *iterator;
        int parent_edge = dual_parent_edge[face];
        if (parent_edge < 0) {
            continue;
        }

        int clockwise_count = 0;
        int parent_start = -1;
        for (const FaceTraversalStep &step : face_collector.faces[face]) {
            if (step.edge == parent_edge) {
                parent_start = step.start;
                continue;
            }
            int source = static_cast<int>(boost::source(graph_edges[step.edge], graph));
            clockwise_count += (direction[step.edge] > 0) == (step.start == source);
        }
        if (parent_start < 0) {
            return OrientationStatus::numerical_failure;
        }
        int parent_source = static_cast<int>(boost::source(graph_edges[parent_edge], graph));
        bool orient_with_face = clockwise_count % 2 == 0;
        direction[parent_edge] = orient_with_face == (parent_start == parent_source) ? 1 : -1;
    }

    for (std::size_t edge = 0; edge < component_edges.size(); edge++) {
        if (direction[edge] == 0) {
            return OrientationStatus::numerical_failure;
        }
        orientation[component_edges[edge]] = direction[edge] > 0;
    }
    return OrientationStatus::success;
}

// Eliminate two rows and columns at a time. The skew-symmetric Schur update
// preserves the Pfaffian of the remaining block, while symmetric row/column
// pivoting changes only its sign. We need only log(abs(Pfaffian)).
bool log_abs_pfaffian(arma::mat matrix, double &log_pfaffian) {
    if (matrix.n_rows != matrix.n_cols || matrix.n_rows % 2 != 0) {
        return false;
    }
    int n_vertices = static_cast<int>(matrix.n_rows);
    log_pfaffian = 0.0;
    for (int eliminated = 0; eliminated < n_vertices; eliminated += 2) {
        int pivot_row = -1;
        int pivot_column = -1;
        double largest_pivot = 0.0;
        for (int row = eliminated; row < n_vertices; row++) {
            for (int column = row + 1; column < n_vertices; column++) {
                double candidate = std::abs(matrix(row, column));
                if (candidate > largest_pivot) {
                    largest_pivot = candidate;
                    pivot_row = row;
                    pivot_column = column;
                }
            }
        }
        if (largest_pivot == 0.0 || !std::isfinite(largest_pivot)) {
            return false;
        }

        matrix.swap_rows(eliminated, pivot_row);
        matrix.swap_cols(eliminated, pivot_row);
        if (pivot_column == eliminated) {
            pivot_column = pivot_row;
        }
        matrix.swap_rows(eliminated + 1, pivot_column);
        matrix.swap_cols(eliminated + 1, pivot_column);
        double pivot = matrix(eliminated, eliminated + 1);
        log_pfaffian += std::log(std::abs(pivot));
        if (!std::isfinite(log_pfaffian)) {
            return false;
        }

        for (int row = eliminated + 2; row < n_vertices; row++) {
            for (int column = row + 1; column < n_vertices; column++) {
                double update = (matrix(row, eliminated) * matrix(eliminated + 1, column) -
                                 matrix(row, eliminated + 1) * matrix(eliminated, column)) /
                                pivot;
                if (!std::isfinite(update)) {
                    return false;
                }
                matrix(row, column) += update;
                if (!std::isfinite(matrix(row, column))) {
                    return false;
                }
                matrix(column, row) = -matrix(row, column);
            }
        }
    }
    return true;
}
} // namespace

HDimerMatchingResult
sample_weighted_perfect_matching(int n_vertices,
                                 const std::vector<HDimerMatchingEdge> &input_edges,
                                 bool has_perfect_matching_witness) {
    if (n_vertices < 0 || n_vertices % 2 != 0) {
        return matching_failure(HDimerMatchingStatus::invalid_input);
    }
    if (n_vertices == 0) {
        HDimerMatchingResult result;
        result.status = HDimerMatchingStatus::success;
        result.planarity_checked = true;
        result.planar = true;
        result.log_partition = 0.0;
        result.max_probability_error = 0.0;
        result.min_probability = 1.0;
        return result;
    }

    std::vector<int> pair_index(static_cast<std::size_t>(n_vertices) * n_vertices, -1);
    std::vector<AggregatedEdge> edges;
    for (const HDimerMatchingEdge &input_edge : input_edges) {
        if (input_edge.u < 0 || input_edge.v < 0 || input_edge.u >= n_vertices ||
            input_edge.v >= n_vertices || input_edge.u == input_edge.v ||
            !std::isfinite(input_edge.log_weight)) {
            return matching_failure(HDimerMatchingStatus::invalid_input);
        }
        int left = std::min(input_edge.u, input_edge.v);
        int right = std::max(input_edge.u, input_edge.v);
        int &index = pair_index[static_cast<std::size_t>(left) * n_vertices + right];
        if (index < 0) {
            index = static_cast<int>(edges.size());
            pair_index[static_cast<std::size_t>(right) * n_vertices + left] = index;
            edges.push_back({left, right, -std::numeric_limits<double>::infinity(), {}});
        }
        edges[index].log_weight = log_sum_exp(edges[index].log_weight, input_edge.log_weight);
        if (!std::isfinite(edges[index].log_weight)) {
            return matching_failure(HDimerMatchingStatus::numerical_failure);
        }
        edges[index].originals.push_back(input_edge);
    }
    if (edges.empty()) {
        return matching_failure(HDimerMatchingStatus::no_perfect_matching, true, true);
    }

    // Lexicographic insertion order substantially reduces Boost's embedding
    // cost on the small fragment graphs without changing any matching weight.
    std::sort(edges.begin(), edges.end(),
              [](const AggregatedEdge &left, const AggregatedEdge &right) {
                  return std::tie(left.u, left.v) < std::tie(right.u, right.v);
              });
    std::fill(pair_index.begin(), pair_index.end(), -1);
    for (int index = 0; index < static_cast<int>(edges.size()); index++) {
        pair_index[static_cast<std::size_t>(edges[index].u) * n_vertices + edges[index].v] =
            index;
        pair_index[static_cast<std::size_t>(edges[index].v) * n_vertices + edges[index].u] =
            index;
    }

    // The split-dimer caller supplies its old matching as a witness. On these
    // small graphs subset DP is exact and has no reason to test planarity.
    if (has_perfect_matching_witness && n_vertices <= SUBSET_DP_MAX_VERTICES) {
        return sample_subset_dp(n_vertices, edges, pair_index, false, false);
    }

    std::vector<std::vector<int>> adjacency(n_vertices);
    for (std::size_t edge = 0; edge < edges.size(); edge++) {
        adjacency[edges[edge].u].push_back(static_cast<int>(edge));
        adjacency[edges[edge].v].push_back(static_cast<int>(edge));
    }

    std::vector<bool> vertex_seen(n_vertices, false);
    std::vector<bool> edge_seen(edges.size(), false);
    std::vector<bool> orientation(edges.size(), false);
    std::vector<double> gauge(n_vertices, 0.0);
    for (int root = 0; root < n_vertices; root++) {
        if (vertex_seen[root]) {
            continue;
        }

        std::vector<int> component_vertices;
        std::vector<int> component_edges;
        std::queue<int> queue;
        vertex_seen[root] = true;
        queue.push(root);
        double component_max_log_weight = -std::numeric_limits<double>::infinity();
        while (!queue.empty()) {
            int vertex = queue.front();
            queue.pop();
            component_vertices.push_back(vertex);
            for (int edge : adjacency[vertex]) {
                if (!edge_seen[edge]) {
                    edge_seen[edge] = true;
                    component_edges.push_back(edge);
                    component_max_log_weight =
                        std::max(component_max_log_weight, edges[edge].log_weight);
                }
                int neighbor = edges[edge].u == vertex ? edges[edge].v : edges[edge].u;
                if (!vertex_seen[neighbor]) {
                    vertex_seen[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }

        OrientationStatus orientation_status = build_clockwise_odd_orientation(
            component_vertices, component_edges, edges, orientation);
        if (orientation_status == OrientationStatus::nonplanar) {
            return matching_failure(HDimerMatchingStatus::nonplanar, true, false);
        }
        if (orientation_status == OrientationStatus::numerical_failure) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
        }
        if (std::isfinite(component_max_log_weight)) {
            for (int vertex : component_vertices) {
                gauge[vertex] = 0.5 * component_max_log_weight;
            }
        }
    }

    if (!has_perfect_matching_witness) {
        PlanarGraph full_graph(n_vertices);
        for (const AggregatedEdge &edge : edges) {
            add_edge(edge.u, edge.v, full_graph);
        }
        std::vector<PlanarGraph::vertex_descriptor> mate(n_vertices);
        auto mate_map = boost::make_iterator_property_map(mate.begin(),
                                                          get(boost::vertex_index, full_graph));
        boost::edmonds_maximum_cardinality_matching(full_graph, mate_map);
        if (2 * boost::matching_size(full_graph, mate_map) !=
            static_cast<std::size_t>(n_vertices)) {
            return matching_failure(HDimerMatchingStatus::no_perfect_matching, true, true);
        }
    }

    if (n_vertices <= SUBSET_DP_MAX_VERTICES) {
        return sample_subset_dp(n_vertices, edges, pair_index, true, true);
    }

    // Subtracting the same half-maximum gauge from every vertex in a connected
    // component multiplies all of its perfect matchings by one common factor.
    // Every exponent is then at most one; underflow is treated as a failure.
    arma::mat kasteleyn_matrix(n_vertices, n_vertices, arma::fill::zeros);
    for (std::size_t edge = 0; edge < edges.size(); edge++) {
        double scaled_log_weight =
            edges[edge].log_weight - gauge[edges[edge].u] - gauge[edges[edge].v];
        double weight = std::exp(scaled_log_weight);
        if (!std::isfinite(weight) || weight == 0.0) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
        }
        kasteleyn_matrix(edges[edge].u, edges[edge].v) = orientation[edge] ? weight : -weight;
        kasteleyn_matrix(edges[edge].v, edges[edge].u) =
            -kasteleyn_matrix(edges[edge].u, edges[edge].v);
    }

    double log_partition = 0.0;
    if (!log_abs_pfaffian(kasteleyn_matrix, log_partition)) {
        return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
    }
    log_partition += std::accumulate(gauge.begin(), gauge.end(), 0.0);
    if (!std::isfinite(log_partition)) {
        return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
    }

    arma::mat inverse;
    if (!arma::inv(inverse, kasteleyn_matrix) || !inverse.is_finite()) {
        return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
    }

    std::vector<int> active_vertices(n_vertices);
    std::iota(active_vertices.begin(), active_vertices.end(), 0);
    std::vector<int> selected_payloads;
    selected_payloads.reserve(n_vertices / 2);
    arma::mat conditional_inverse = inverse;
    double max_probability_error = 0.0;
    double min_probability = 1.0;
    while (!active_vertices.empty()) {
        std::vector<int> candidate_positions;
        std::vector<double> probabilities;
        double probability_sum = 0.0;
        for (int position = 1; position < static_cast<int>(active_vertices.size());
             position++) {
            double probability =
                kasteleyn_matrix(active_vertices[0], active_vertices[position]) *
                conditional_inverse(position, 0);
            if (!std::isfinite(probability) || probability < -NEGATIVE_PROBABILITY_TOLERANCE) {
                return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
            }
            if (probability > 0.0) {
                candidate_positions.push_back(position);
                probabilities.push_back(probability);
                probability_sum += probability;
            }
        }
        double probability_error = std::abs(probability_sum - 1.0);
        max_probability_error = std::max(max_probability_error, probability_error);
        if (candidate_positions.empty() || !(probability_sum > 0.0) ||
            !std::isfinite(probability_sum) || probability_error > PROBABILITY_SUM_TOLERANCE) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
        }
        for (double &probability : probabilities) {
            probability /= probability_sum;
            min_probability = std::min(min_probability, probability);
        }

        double uniform = R::runif(0, 1);
        double cumulative = 0.0;
        int selected_position = candidate_positions.back();
        for (std::size_t candidate = 0; candidate < candidate_positions.size(); candidate++) {
            cumulative += probabilities[candidate];
            if (uniform <= cumulative) {
                selected_position = candidate_positions[candidate];
                break;
            }
        }

        int first_vertex = active_vertices[0];
        int second_vertex = active_vertices[selected_position];
        int selected_edge =
            pair_index[static_cast<std::size_t>(first_vertex) * n_vertices + second_vertex];
        if (selected_edge < 0) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
        }
        int payload = -1;
        if (!sample_parallel_payload(edges[selected_edge], payload)) {
            return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
        }
        selected_payloads.push_back(payload);

        std::vector<int> retained_positions;
        retained_positions.reserve(active_vertices.size() - 2);
        for (int position = 0; position < static_cast<int>(active_vertices.size());
             position++) {
            if (position != 0 && position != selected_position) {
                retained_positions.push_back(position);
            }
        }
        if (!retained_positions.empty()) {
            double pivot = conditional_inverse(0, selected_position);
            if (!std::isfinite(pivot) || pivot == 0.0) {
                return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
            }

            // This skew Schur update is the inverse Kasteleyn matrix after
            // conditioning on the selected edge and deleting its endpoints.
            arma::mat next_inverse(retained_positions.size(), retained_positions.size(),
                                   arma::fill::zeros);
            for (int row = 0; row < static_cast<int>(retained_positions.size()); row++) {
                for (int column = row + 1; column < static_cast<int>(retained_positions.size());
                     column++) {
                    double value =
                        conditional_inverse(retained_positions[row],
                                            retained_positions[column]) +
                        (conditional_inverse(retained_positions[row], 0) *
                             conditional_inverse(selected_position,
                                                 retained_positions[column]) -
                         conditional_inverse(retained_positions[row], selected_position) *
                             conditional_inverse(0, retained_positions[column])) /
                            pivot;
                    if (!std::isfinite(value)) {
                        return matching_failure(HDimerMatchingStatus::numerical_failure, true,
                                                true);
                    }
                    next_inverse(row, column) = value;
                    next_inverse(column, row) = -value;
                }
            }
            conditional_inverse = std::move(next_inverse);
            if (!conditional_inverse.is_finite()) {
                return matching_failure(HDimerMatchingStatus::numerical_failure, true, true);
            }
        }

        std::vector<int> next_active_vertices;
        next_active_vertices.reserve(retained_positions.size());
        for (int position : retained_positions) {
            next_active_vertices.push_back(active_vertices[position]);
        }
        active_vertices = std::move(next_active_vertices);
    }

    HDimerMatchingResult result;
    result.status = HDimerMatchingStatus::success;
    result.planarity_checked = true;
    result.planar = true;
    result.selected_edge_indices = std::move(selected_payloads);
    result.log_partition = log_partition;
    result.max_probability_error = max_probability_error;
    result.min_probability = min_probability;
    return result;
}

// [[Rcpp::export]]
Rcpp::List hdimer_matching_diagnostic(int n_vertices, Rcpp::IntegerVector edge_u,
                                      Rcpp::IntegerVector edge_v, Rcpp::NumericVector weights) {
    if (edge_u.size() != edge_v.size() || edge_u.size() != weights.size()) {
        Rcpp::stop("Edge vectors must have equal lengths.");
    }
    std::vector<HDimerMatchingEdge> edges;
    for (int index = 0; index < weights.size(); index++) {
        if (!std::isfinite(weights[index]) || weights[index] < 0.0) {
            Rcpp::stop("Weights must be finite and non-negative.");
        }
        if (weights[index] > 0.0) {
            edges.push_back(
                {edge_u[index] - 1, edge_v[index] - 1, std::log(weights[index]), index + 1});
        }
    }

    HDimerMatchingResult result = sample_weighted_perfect_matching(n_vertices, edges);
    Rcpp::LogicalVector planar = Rcpp::LogicalVector::create(NA_LOGICAL);
    if (result.planarity_checked) {
        planar[0] = result.planar;
    }
    bool success = result.status == HDimerMatchingStatus::success;
    bool numerical_failure = result.status == HDimerMatchingStatus::numerical_failure;
    return Rcpp::List::create(
        Rcpp::Named("sampled_edge_indices") = result.selected_edge_indices,
        Rcpp::Named("status") = matching_status_label(result.status),
        Rcpp::Named("success") = success,
        Rcpp::Named("planarity_checked") = result.planarity_checked,
        Rcpp::Named("planar") = planar, Rcpp::Named("numerical_failure") = numerical_failure,
        Rcpp::Named("used_subset_dp") = result.used_subset_dp,
        Rcpp::Named("log_partition") = result.log_partition,
        Rcpp::Named("max_probability_error") = result.max_probability_error,
        Rcpp::Named("min_probability") = result.min_probability);
}
