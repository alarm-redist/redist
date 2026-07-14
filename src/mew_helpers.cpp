/*
 * MEW Helper Functions Implementation
 *
 * Core utilities for the Marked Edge Walk algorithm
 */

#include "mew_helpers.h"
#include "wilson.h"  // For sample_sub_ust, init_tree, init_multigraph
#include <algorithm>
#include <cmath>

namespace {

void add_tree_edge(Tree &tree, int u, int v) {
    // Add v to u's adjacency list
    if (std::find(tree[u].begin(), tree[u].end(), v) == tree[u].end()) {
        tree[u].push_back(v);
    }
    // Add u to v's adjacency list
    if (std::find(tree[v].begin(), tree[v].end(), u) == tree[v].end()) {
        tree[v].push_back(u);
    }
}

bool incident(const Edge &edge, int vertex) {
    return edge.first == vertex || edge.second == vertex;
}

int old_degree(const MEWTree &tree_new, int vertex,
               const Edge &edge_plus, const Edge &edge_minus) {
    int degree = tree_new.degree(vertex);
    if (edge_plus != edge_minus) {
        degree -= incident(edge_plus, vertex);
        degree += incident(edge_minus, vertex);
    }
    return degree;
}

void dfs_component(const MEWTree &tree, int u, std::vector<bool> &visited,
                   std::vector<int> &component,
                   const MarkedEdgeSet &marked_edges) {
    visited[u] = true;
    component.push_back(u);

    for (int v : tree.neighbors(u)) {
        Edge e = make_edge(u, v);
        if (marked_edges.find(e) == marked_edges.end() && !visited[v]) {
            dfs_component(tree, v, visited, component, marked_edges);
        }
    }
}

} // namespace

MEWTree::MEWTree(Tree tree)
    : lct_(tree.size()), adjacency_(std::move(tree)) {
    for (int u = 0; u < size(); u++) {
        for (int v : adjacency_[u]) {
            if (u < v) {
                edges_.insert(make_edge(u, v));
                lct_.evert(u);
                lct_.link(u, v);
            }
        }
    }
}

int MEWTree::size() const {
    return static_cast<int>(adjacency_.size());
}

const std::vector<int>& MEWTree::neighbors(int vertex) const {
    return adjacency_[vertex];
}

int MEWTree::degree(int vertex) const {
    return static_cast<int>(adjacency_[vertex].size());
}

bool MEWTree::has_edge(int u, int v) const {
    return edges_.find(make_edge(u, v)) != edges_.end();
}

std::vector<Edge> MEWTree::find_cycle(int u, int v) {
    lct_.evert(u);
    std::vector<int> path = lct_.find_path(v);
    if (path.empty() || path.front() != u || path.back() != v) {
        Rcpp::stop("Could not find path between vertices in MEW tree");
    }

    std::vector<Edge> cycle_edges;
    cycle_edges.reserve(path.size());
    for (std::size_t i = path.size() - 1; i > 0; i--) {
        cycle_edges.push_back(make_edge(path[i], path[i - 1]));
    }
    cycle_edges.push_back(make_edge(u, v));
    return cycle_edges;
}

std::size_t MEWTree::remove_neighbor(std::vector<int>& neighbors, int vertex) {
    auto it = std::find(neighbors.begin(), neighbors.end(), vertex);
    if (it == neighbors.end()) {
        Rcpp::stop("MEW tree edge not found");
    }
    std::size_t position = std::distance(neighbors.begin(), it);
    neighbors.erase(it);
    return position;
}

MEWTreeUpdate MEWTree::replace_edge(const Edge& edge_plus,
                                    const Edge& edge_minus) {
    MEWTreeUpdate update{edge_plus, edge_minus, 0, 0,
                         edge_plus != edge_minus};
    if (!update.changed) {
        return update;
    }

    adjacency_[edge_plus.first].push_back(edge_plus.second);
    adjacency_[edge_plus.second].push_back(edge_plus.first);
    update.minus_pos_first = remove_neighbor(
        adjacency_[edge_minus.first], edge_minus.second);
    update.minus_pos_second = remove_neighbor(
        adjacency_[edge_minus.second], edge_minus.first);
    edges_.erase(edge_minus);
    edges_.insert(edge_plus);

    lct_.evert(edge_minus.first);
    lct_.cut(edge_minus.second);
    lct_.evert(edge_plus.first);
    lct_.link(edge_plus.first, edge_plus.second);
    return update;
}

void MEWTree::rollback(const MEWTreeUpdate& update) {
    if (!update.changed) {
        return;
    }

    remove_neighbor(adjacency_[update.edge_plus.first], update.edge_plus.second);
    remove_neighbor(adjacency_[update.edge_plus.second], update.edge_plus.first);
    adjacency_[update.edge_minus.first].insert(
        adjacency_[update.edge_minus.first].begin() + update.minus_pos_first,
        update.edge_minus.second);
    adjacency_[update.edge_minus.second].insert(
        adjacency_[update.edge_minus.second].begin() + update.minus_pos_second,
        update.edge_minus.first);
    edges_.erase(update.edge_plus);
    edges_.insert(update.edge_minus);

    lct_.evert(update.edge_plus.first);
    lct_.cut(update.edge_plus.second);
    lct_.evert(update.edge_minus.first);
    lct_.link(update.edge_minus.first, update.edge_minus.second);
}

std::vector<std::vector<int>> tree_components_list(const MEWTree &tree,
                                                    const MarkedEdgeSet &marked_edges) {
    int V = tree.size();
    std::vector<bool> visited(V, false);
    std::vector<std::vector<int>> components;

    for (int u = 0; u < V; u++) {
        if (!visited[u]) {
            std::vector<int> component;
            dfs_component(tree, u, visited, component, marked_edges);
            components.push_back(component);
        }
    }

    return components;
}

uvec tree_to_partition(const MEWTree &tree, const MarkedEdgeSet &marked_edges,
                       int V, int n_distr) {
    auto components = tree_components_list(tree, marked_edges);

    // Wrong component count can happen when a duplicate marked edge shrinks
    // the set.  The population check in mew_proposal will reject these.

    uvec partition(V);
    for (size_t i = 0; i < components.size(); i++) {
        for (int u : components[i]) {
            partition(u) = i + 1;  // 1-indexed for R
        }
    }

    return partition;
}


/*
 * Transition probability
 */

double transition_probability(const std::vector<Edge> &cycle_edges,
                             const Edge &edge_plus,
                             const Edge &marked_old,
                             const Edge &marked_new,
                             const MarkedEdgeSet &marked_edges_old,
                             const MarkedEdgeSet &marked_edges_new,
                             const Edge &edge_minus,
                             const MEWTree &tree_new) {
    // Check if new marked edge equals edge_plus (invalid proposal)
    if (marked_new == edge_plus) {
        return 0.0;
    }

    // Extract vertices
    int w = marked_old.first;
    int x = marked_old.second;
    int u = marked_new.first;
    int v = marked_new.second;

    // Compute pm: accounts for degree changes in marked edge step
    // Formula from McWhorter & DeFord (2024), Section 3.2
    //
    // NOTE: The paper's formula assumes sampling neighbors from the NEW tree (T').
    // The Julia implementation samples from the OLD tree (T) but uses
    // the NEW tree formula.
    // If correcting this in the future, the formula for OLD tree sampling would
    // invert the degree ratios: d_u/d_u_p instead of d_u_p/d_u.
    double pm = 1.0;

    if ((u == w && v == x) || (u == x && v == w)) {
        // Both endpoints same (identity move on marked edges)
        int d_u = old_degree(tree_new, u, edge_plus, edge_minus);
        int d_u_p = tree_new.degree(u);
        int d_v = old_degree(tree_new, v, edge_plus, edge_minus);
        int d_v_p = tree_new.degree(v);

        pm = ((double)d_u_p / d_u) *
             ((double)(d_u + d_v) / (d_u_p + d_v_p)) *
             ((double)d_v_p / d_v);
    } else {
        // One endpoint shared
        int shared = -1;
        if (u == w || u == x) {
            shared = u;
        } else if (v == w || v == x) {
            shared = v;
        }

        if (shared >= 0) {
            int d_u = old_degree(tree_new, shared, edge_plus, edge_minus);
            int d_u_p = tree_new.degree(shared);
            pm = (double)d_u_p / d_u;
        }
    }

    // Compute pt: accounts for cycle/marked edge intersection
    // Count edges in cycle but not in marked_old
    int l = 0;
    for (const Edge &e : cycle_edges) {
        if (marked_edges_old.find(e) == marked_edges_old.end()) {
            l++;
        }
    }

    // Count edges in cycle but not in marked_new
    int l_p = 0;
    for (const Edge &e : cycle_edges) {
        if (marked_edges_new.find(e) == marked_edges_new.end()) {
            l_p++;
        }
    }

    double pt = (l_p > 0) ? ((double)l / l_p) : 0.0;

    // MH ratio: q(x'->x)/q(x->x') = pt * pm
    // pt = |C\M| / |C\M'| corrects for cycle step asymmetry.
    // pm = deg_T'(v)/deg_T(v) corrects for marked step degree asymmetry.
    return pt * pm;
}

/*
 * Quotient graph spanning tree count
 * Computes log(τ(Q)) where Q is the inter-district multigraph.
 * Uses Kirchhoff's matrix-tree theorem: τ(Q) = det(L*) where L* is
 * any cofactor of the Laplacian of Q.
 */
double log_st_quotient_graph(const Graph &g, const arma::uvec &plan, int n_distr) {
    arma::mat L(n_distr, n_distr, arma::fill::zeros);
    int V = g.size();

    for (int u = 0; u < V; u++) {
        int du = plan(u);
        for (int v : g[u]) {
            if (u < v) {
                int dv = plan(v);
                if (du != dv) {
                    L(du - 1, du - 1) += 1.0;
                    L(dv - 1, dv - 1) += 1.0;
                    L(du - 1, dv - 1) -= 1.0;
                    L(dv - 1, du - 1) -= 1.0;
                }
            }
        }
    }

    // Kirchhoff: τ = det of any (n-1)×(n-1) cofactor of L
    arma::mat L_minor = L.submat(0, 0, n_distr - 2, n_distr - 2);
    double det_val = arma::det(L_minor);
    if (det_val <= 0) return 0.0;
    return std::log(det_val);
}

/*
 * Proposal mechanisms
 */

CycleProposal cycle_basis_step(const Graph &g, MEWTree &tree,
                               const MarkedEdgeSet &marked_edges) {
    // Find edges in g but not in tree
    std::vector<Edge> non_tree_edges;

    for (size_t u = 0; u < g.size(); u++) {
        for (int v : g[u]) {
            if ((int)u < v) {  // Only consider each edge once
                Edge e = make_edge(u, v);
                if (!tree.has_edge(u, v)) {
                    non_tree_edges.push_back(e);
                }
            }
        }
    }

    if (non_tree_edges.empty()) {
        Rcpp::stop("No non-tree edges available - tree spans entire graph");
    }

    // Sample one edge to add
    int idx = r_int(non_tree_edges.size());
    Edge edge_plus = non_tree_edges[idx];

    // Find cycle formed by adding this edge to tree
    std::vector<Edge> cycle_edges = tree.find_cycle(
        edge_plus.first, edge_plus.second);

    // Find edges in cycle that are not marked (including edge_plus).
    // Choosing edge_minus == edge_plus is a valid "identity" move on the tree.
    std::vector<Edge> possible_cuts;
    for (const Edge &e : cycle_edges) {
        if (marked_edges.find(e) == marked_edges.end()) {
            possible_cuts.push_back(e);
        }
    }

    if (possible_cuts.empty()) {
        // All cycle edges are marked - return no-op proposal
        // This is valid - just means we keep the current state
        CycleProposal proposal;
        proposal.cycle_edges = cycle_edges;
        proposal.edge_plus = edge_plus;
        proposal.edge_minus = make_edge(0, 0);  // Invalid edge
        proposal.update = MEWTreeUpdate{
            edge_plus, proposal.edge_minus, 0, 0, false};
        proposal.valid = false;
        return proposal;
    }

    // Sample edge to remove
    int cut_idx = r_int(possible_cuts.size());
    Edge edge_minus = possible_cuts[cut_idx];

    CycleProposal proposal;
    proposal.cycle_edges = cycle_edges;
    proposal.edge_plus = edge_plus;
    proposal.edge_minus = edge_minus;
    proposal.update = tree.replace_edge(edge_plus, edge_minus);
    proposal.valid = true;

    return proposal;
}

MarkedEdgeProposal marked_edge_step(const MEWTree &tree,
                                   const MarkedEdgeSet &marked_edges) {
    if (marked_edges.empty()) {
        Rcpp::stop("No marked edges to update");
    }

    // Sample one marked edge
    int idx = r_int(marked_edges.size());
    auto it = marked_edges.begin();
    std::advance(it, idx);
    Edge old_edge = *it;

    // Sample one endpoint of the marked edge
    int endpoint_choice = r_int(2);
    int chosen_vertex = (endpoint_choice == 0) ? old_edge.first : old_edge.second;

    // Sample one neighbor uniformly at random.
    // If new_edge == old_edge, it's a valid identity move.
    const std::vector<int>& neighbors = tree.neighbors(chosen_vertex);
    if (neighbors.empty()) {
        Rcpp::stop("Chosen vertex has no neighbors in tree");
    }
    int neighbor_idx = r_int(neighbors.size());
    int neighbor = neighbors[neighbor_idx];
    Edge new_edge = make_edge(chosen_vertex, neighbor);

    // Update marked edge set: remove old, insert new
    MarkedEdgeSet marked_new = marked_edges;
    marked_new.erase(old_edge);
    marked_new.insert(new_edge);

    MarkedEdgeProposal proposal;
    proposal.old_edge = old_edge;
    proposal.new_edge = new_edge;
    proposal.marked_new = marked_new;

    return proposal;
}

MEWProposal mew_proposal(const Graph &g, MEWTree &tree,
                        const MarkedEdgeSet &marked_edges,
                        const uvec &pop, int n_distr,
                        double target, double lower, double upper) {
    int V = g.size();
    // Invalid proposals become self-loops rather than being redrawn.
    MEWProposal proposal;
    proposal.n_rejects = 0;
    proposal.valid = false;
    proposal.partition = tree_to_partition(tree, marked_edges, V, n_distr);

    CycleProposal cycle_prop = cycle_basis_step(g, tree, marked_edges);
    proposal.cycle = cycle_prop;
    if (!cycle_prop.valid) {
        proposal.marked.old_edge = make_edge(0, 0);
        proposal.marked.new_edge = make_edge(0, 0);
        proposal.marked.marked_new = marked_edges;
        return proposal;
    }

    MarkedEdgeProposal marked_prop = marked_edge_step(tree, marked_edges);
    proposal.marked = marked_prop;

    if (marked_prop.marked_new.size() != marked_edges.size()) {
        return proposal;
    }

    auto components = tree_components_list(tree, marked_prop.marked_new);
    if ((int)components.size() != n_distr) {
        return proposal;
    }

    uvec partition = tree_to_partition(
        tree, marked_prop.marked_new, V, n_distr);
    proposal.partition = partition;

    std::vector<double> dist_pop(n_distr, 0.0);
    for (int i = 0; i < V; i++) {
        int dist = partition(i) - 1;
        dist_pop[dist] += pop(i);
    }

    for (int d = 0; d < n_distr; d++) {
        if (dist_pop[d] < lower || dist_pop[d] > upper) {
            return proposal;
        }
    }

    proposal.valid = true;
    return proposal;
}

/*
 * Initialization
 */

std::pair<Tree, MarkedEdgeSet> partition_to_tree_marked_edges(
    const Graph &g,
    const uvec &partition,
    int n_distr
) {
    int V = g.size();

    // STEP 0: Convert partition to district lists
    std::vector<std::vector<int>> districts(n_distr);
    for (int i = 0; i < V; i++) {
        int dist = partition(i) - 1;  // Convert to 0-indexed
        if (dist >= 0 && dist < n_distr) {
            districts[dist].push_back(i);
        }
    }

    Tree tree = init_tree(V);

    // STEP 1: Build spanning tree within each district using Wilson's algorithm
    for (int d = 0; d < n_distr; d++) {
        const auto &district = districts[d];
        int dist_size = district.size();

        if (dist_size <= 1) {
            continue;  // Single-vertex district, no edges needed
        }

        // Create vertex mapping: original -> subgraph index
        std::map<int, int> vertex_map;
        std::vector<int> reverse_map(dist_size);
        for (size_t i = 0; i < district.size(); i++) {
            vertex_map[district[i]] = i;
            reverse_map[i] = district[i];
        }

        // Build induced subgraph for this district
        Graph subgraph(dist_size);
        for (int u_orig : district) {
            int u_sub = vertex_map[u_orig];
            for (int v_orig : g[u_orig]) {
                if (vertex_map.count(v_orig) > 0) {  // v also in this district
                    int v_sub = vertex_map[v_orig];
                    if (u_sub < v_sub) {  // Only add each edge once
                        // Add edge in both directions (undirected)
                        if (std::find(subgraph[u_sub].begin(), subgraph[u_sub].end(), v_sub)
                            == subgraph[u_sub].end()) {
                            subgraph[u_sub].push_back(v_sub);
                            subgraph[v_sub].push_back(u_sub);
                        }
                    }
                }
            }
        }

        // Run Wilson's algorithm on subgraph
        Tree subtree = init_tree(dist_size);
        std::vector<bool> visited(dist_size, false);
        std::vector<bool> ignore(dist_size, false);
        uvec dummy_pop = arma::ones<uvec>(dist_size);
        uvec dummy_counties = arma::ones<uvec>(dist_size);
        Multigraph mg = init_multigraph(dist_size);
        int root = 0;  // Root for Wilson's algorithm

        int result = sample_sub_ust(subgraph, subtree, dist_size, root, visited, ignore,
                                   dummy_pop, 0, std::numeric_limits<double>::max(),
                                   dummy_counties, mg);

        if (result != 0) {
            std::ostringstream msg;
            msg << "Failed to generate spanning tree for district " << d
                << " (size " << dist_size << ")";
            Rcpp::stop(msg.str());
        }

        // Add subtree edges to global tree (map back to original vertices)
        for (int u_sub = 0; u_sub < dist_size; u_sub++) {
            for (int v_sub : subtree[u_sub]) {
                int u_orig = reverse_map[u_sub];
                int v_orig = reverse_map[v_sub];
                add_tree_edge(tree, u_orig, v_orig);
            }
        }
    }

    // STEP 2: Find all cut edges (edges crossing district boundaries)
    std::vector<Edge> cut_edges;
    for (int u = 0; u < V; u++) {
        for (int v : g[u]) {
            if (u < v && partition(u) != partition(v)) {
                cut_edges.push_back(make_edge(u, v));
            }
        }
    }

    // STEP 3: Build minimal spanning tree of districts using cut edges
    // This creates exactly k-1 marked edges that connect the k district trees
    std::set<int> connected_districts;
    connected_districts.insert(0);  // Start with district 0

    MarkedEdgeSet marked_edges;

    while ((int)connected_districts.size() < n_distr) {
        bool found = false;

        for (auto it = cut_edges.begin(); it != cut_edges.end(); ++it) {
            Edge e = *it;
            int dist1 = partition(e.first) - 1;   // 0-indexed
            int dist2 = partition(e.second) - 1;

            // Check if exactly one endpoint's district is in connected set (XOR)
            bool in1 = connected_districts.count(dist1) > 0;
            bool in2 = connected_districts.count(dist2) > 0;

            if (in1 != in2) {  // XOR: exactly one is connected
                // This edge connects a new district to the connected component
                marked_edges.insert(e);
                connected_districts.insert(dist1);
                connected_districts.insert(dist2);
                add_tree_edge(tree, e.first, e.second);
                cut_edges.erase(it);  // Remove from candidates
                found = true;
                break;
            }
        }

        if (!found) {
            // This should never happen if partition is valid and graph is connected
            std::ostringstream msg;
            msg << "Failed to connect all districts. Connected: "
                << connected_districts.size() << " of " << n_distr;
            Rcpp::stop(msg.str());
        }
    }

    // Verify we have exactly k-1 marked edges
    if ((int)marked_edges.size() != n_distr - 1) {
        std::ostringstream msg;
        msg << "partition_to_tree_marked_edges: Expected " << (n_distr - 1)
            << " marked edges but got " << marked_edges.size();
        Rcpp::stop(msg.str());
    }

    return std::make_pair(tree, marked_edges);
}
