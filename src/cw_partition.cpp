/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * LCT Partition Implementation
 ********************************************************/

#include "cw_partition.h"
#include <queue>
#include <stdexcept>

// ============================================================
// CWEdge Implementation
// ============================================================

CWEdge::CWEdge(int a, int b, double w)
    : u(std::min(a, b)), v(std::max(a, b)), weight(w) {}

bool CWEdge::operator<(const CWEdge& other) const {
    if (u != other.u) return u < other.u;
    if (v != other.v) return v < other.v;
    return weight < other.weight;
}

bool CWEdge::operator==(const CWEdge& other) const {
    return u == other.u && v == other.v;
}

// ============================================================
// LCTPartition Implementation
// ============================================================

const EdgeSet LCTPartition::empty_edge_set;

LCTPartition::LCTPartition(int n_vertices, int n_districts)
    : n_districts(n_districts),
      n_vertices(n_vertices),
      lct(n_vertices),
      district_roots(n_districts, -1),
      node_to_district(n_vertices, -1),
      district_pop(n_districts, 0),
      graph(nullptr),
      pop(nullptr),
      counties(nullptr),
      default_edge_weight_(1.0) {}

int LCTPartition::init_from_plan(const Graph& g,
                                  const arma::uvec& plan,
                                  const arma::uvec& population,
                                  const arma::uvec& county_assignments,
                                  double lower, double upper) {
    graph = &g;
    pop = &population;
    counties = &county_assignments;

    int V = n_vertices;

    // For each district, sample a spanning tree using Wilson's algorithm
    // and load it into the LCT
    Multigraph cg = county_graph(g, county_assignments);

    for (int d = 0; d < n_districts; d++) {
        // Create ignore mask: ignore vertices not in this district
        std::vector<bool> ignore(V);
        int n_in_district = 0;
        for (int i = 0; i < V; i++) {
            if ((int)plan(i) == d + 1) {  // plan is 1-indexed
                ignore[i] = false;
                n_in_district++;
            } else {
                ignore[i] = true;
            }
        }

        if (n_in_district == 0) {
            Rcpp::stop("District %d has no vertices", d + 1);
        }

        // Sample spanning tree for this district
        Tree tree = init_tree(V);
        int root;
        std::vector<bool> visited(V);
        int result = sample_sub_ust(g, tree, V, root, visited, ignore,
                                    population, lower, upper, county_assignments, cg);
        if (result != 0) {
            // Wilson's algorithm failed - return error code
            return 1;
        }

        // Store root
        district_roots[d] = root;

        // Load tree into LCT and assign district labels
        load_tree_into_lct(tree, root, d);

        // Calculate district population
        district_pop[d] = 0;
        for (int i = 0; i < V; i++) {
            if ((int)plan(i) == d + 1) {
                district_pop[d] += population(i);
            }
        }
    }

    // Find all cross-district edges
    find_cross_district_edges();

    return 0;
}

void LCTPartition::load_tree_into_lct(const Tree& tree, int root, int district) {
    // BFS from root to load edges into LCT
    std::queue<int> queue;
    queue.push(root);

    // First mark root's district
    node_to_district[root] = district;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();

        // tree[v] contains children of v (edges pointing away from root)
        for (int child : tree[v]) {
            // Link child to parent v in the LCT
            // First make child its own root (it should be already)
            lct.evert(child);
            lct.link(child, v);

            // Assign district
            node_to_district[child] = district;

            queue.push(child);
        }
    }
}

void LCTPartition::find_cross_district_edges() {
    cross_edges.clear();

    // Scan all edges in the graph
    for (int u = 0; u < n_vertices; u++) {
        int d_u = node_to_district[u];
        for (int v : (*graph)[u]) {
            if (v > u) {  // Only count each edge once
                int d_v = node_to_district[v];
                if (d_u != d_v) {
                    // Cross-district edge
                    DistrictPair key(std::min(d_u, d_v), std::max(d_u, d_v));
                    cross_edges[key].insert(CWEdge(u, v));
                }
            }
        }
    }
}

void LCTPartition::assign_districts_from_root(int root, int district) {
    // Use LCT to traverse all nodes in the tree rooted at root
    lct.evert(root);

    std::queue<LCTNode*> queue;
    LCTNode* root_node = lct.node(root);
    lct.expose(root_node);
    queue.push(root_node);

    while (!queue.empty()) {
        LCTNode* node = queue.front();
        queue.pop();
        node_to_district[node->vertex] = district;

        // Add children and path children
        for (int i = 0; i < 2; i++) {
            if (node->children[i] != nullptr) {
                queue.push(node->children[i]);
            }
        }
        for (LCTNode* pc : node->path_children) {
            queue.push(pc);
        }
    }
}

int LCTPartition::get_district(int v) const {
    return node_to_district[v];
}

int LCTPartition::get_district_pop(int d) const {
    return district_pop[d];
}

const EdgeSet& LCTPartition::get_cross_edges(int d1, int d2) const {
    DistrictPair key(std::min(d1, d2), std::max(d1, d2));
    auto it = cross_edges.find(key);
    if (it != cross_edges.end()) {
        return it->second;
    }
    return empty_edge_set;
}

bool LCTPartition::districts_adjacent(int d1, int d2) const {
    DistrictPair key(std::min(d1, d2), std::max(d1, d2));
    return cross_edges.find(key) != cross_edges.end();
}

std::vector<DistrictPair> LCTPartition::get_adjacent_district_pairs() const {
    std::vector<DistrictPair> pairs;
    pairs.reserve(cross_edges.size());
    for (const auto& kv : cross_edges) {
        pairs.push_back(kv.first);
    }
    return pairs;
}

arma::uvec LCTPartition::get_plan() const {
    arma::uvec plan(n_vertices);
    for (int i = 0; i < n_vertices; i++) {
        plan(i) = node_to_district[i] + 1;  // Convert to 1-indexed
    }
    return plan;
}

void LCTPartition::print_state(int verbosity) const {
    Rcpp::Rcout << "[LCTPartition] " << n_districts << " districts, "
                << n_vertices << " vertices\n";

    if (verbosity >= 1) {
        Rcpp::Rcout << "[LCTPartition] District roots: ";
        for (int i = 0; i < n_districts; i++) {
            Rcpp::Rcout << district_roots[i];
            if (i < n_districts - 1) Rcpp::Rcout << ", ";
        }
        Rcpp::Rcout << "\n";

        Rcpp::Rcout << "[LCTPartition] District populations: ";
        for (int i = 0; i < n_districts; i++) {
            Rcpp::Rcout << district_pop[i];
            if (i < n_districts - 1) Rcpp::Rcout << ", ";
        }
        Rcpp::Rcout << "\n";

        Rcpp::Rcout << "[LCTPartition] Adjacent district pairs: "
                    << cross_edges.size() << "\n";
    }

    if (verbosity >= 2) {
        for (const auto& kv : cross_edges) {
            Rcpp::Rcout << "[LCTPartition] Districts " << kv.first.first + 1
                        << "-" << kv.first.second + 1
                        << ": " << kv.second.size() << " boundary edges\n";
        }
    }
}

// ============================================================
// Edge Weight Methods
// ============================================================

void LCTPartition::set_edge_weights(const Rcpp::List& edge_weights_list) {
    edge_weights_.clear();

    if (edge_weights_list.size() == 0) {
        return;  // No weights specified, use defaults
    }

    // Each element should be a list with $edge and $weight
    for (int i = 0; i < edge_weights_list.size(); i++) {
        Rcpp::List entry = edge_weights_list[i];

        // Extract edge (2-element vector)
        if (!entry.containsElementNamed("edge")) {
            Rcpp::stop("Edge weight entry %d missing 'edge' field", i + 1);
        }
        Rcpp::IntegerVector edge = entry["edge"];
        if (edge.size() != 2) {
            Rcpp::stop("Edge weight entry %d: 'edge' must be length 2", i + 1);
        }

        // Extract weight (scalar)
        if (!entry.containsElementNamed("weight")) {
            Rcpp::stop("Edge weight entry %d missing 'weight' field", i + 1);
        }
        double weight = Rcpp::as<double>(entry["weight"]);

        if (weight <= 0) {
            Rcpp::stop("Edge weight entry %d: weight must be positive", i + 1);
        }

        // Convert from 1-indexed R to 0-indexed C++
        int u = edge[0] - 1;
        int v = edge[1] - 1;

        // Validate vertices
        if (u < 0 || u >= n_vertices || v < 0 || v >= n_vertices) {
            Rcpp::stop("Edge weight entry %d: vertices out of range [1, %d]",
                       i + 1, n_vertices);
        }

        // Check edge exists in graph
        const std::vector<int>& neighbors = (*graph)[u];
        if (std::find(neighbors.begin(), neighbors.end(), v) == neighbors.end()) {
            Rcpp::stop("Edge weight entry %d: edge (%d, %d) not in adjacency graph",
                       i + 1, edge[0], edge[1]);
        }

        // Store weight (canonicalized order)
        auto key = std::make_pair(std::min(u, v), std::max(u, v));
        edge_weights_[key] = weight;
    }
}

double LCTPartition::get_edge_weight(int u, int v) const {
    // Canonicalize edge order (undirected graph)
    auto key = std::make_pair(std::min(u, v), std::max(u, v));

    // Lookup in weight map
    auto it = edge_weights_.find(key);
    if (it != edge_weights_.end()) {
        return it->second;
    }

    // Default weight if not found
    return default_edge_weight_;
}

