/*
 * MEW Helper Functions Implementation
 * 
 * Core utilities for the Marked Edge Walk algorithm
 */

#include "mew_helpers.h"
#include "wilson.h"  // For sample_sub_ust, init_tree, init_multigraph
#include <queue>
#include <algorithm>
#include <cmath>

/*
 * Tree edge operations
 */

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

void remove_tree_edge(Tree &tree, int u, int v) {
    // Remove v from u's adjacency list
    auto it_u = std::find(tree[u].begin(), tree[u].end(), v);
    if (it_u != tree[u].end()) {
        tree[u].erase(it_u);
    }
    // Remove u from v's adjacency list
    auto it_v = std::find(tree[v].begin(), tree[v].end(), u);
    if (it_v != tree[v].end()) {
        tree[v].erase(it_v);
    }
}

bool has_tree_edge(const Tree &tree, int u, int v) {
    // Tree is directed with edges pointing away from root
    // Check both directions to see if edge exists
    return std::find(tree[u].begin(), tree[u].end(), v) != tree[u].end() ||
           std::find(tree[v].begin(), tree[v].end(), u) != tree[v].end();
}

std::vector<Edge> tree_to_edges(const Tree &tree) {
    std::vector<Edge> edges;
    for (size_t u = 0; u < tree.size(); u++) {
        for (int v : tree[u]) {
            // Only add each edge once (u < v)
            if ((int)u < v) {
                edges.push_back(make_edge(u, v));
            }
        }
    }
    return edges;
}

/*
 * Cycle detection using BFS
 */

std::vector<Edge> find_cycle(const Tree &tree, int u, int v) {
    // Find path from u to v in tree using BFS, then add edge (u,v) to complete cycle
    // Tree is directed with edges pointing away from root, so we need to traverse bidirectionally
    int V = tree.size();
    std::vector<int> parent(V, -1);
    std::vector<bool> visited(V, false);
    std::queue<int> q;
    
    // Start BFS from u
    q.push(u);
    visited[u] = true;
    parent[u] = u;  // Mark root
    
    // BFS to find v (traverse tree bidirectionally)
    while (!q.empty()) {
        int curr = q.front();
        q.pop();
        
        if (curr == v) {
            // Found target - reconstruct path
            std::vector<int> path;
            int node = v;
            while (node != u) {
                path.push_back(node);
                node = parent[node];
            }
            path.push_back(u);
            
            // Convert path to edges
            std::vector<Edge> cycle_edges;
            for (size_t i = 0; i < path.size() - 1; i++) {
                cycle_edges.push_back(make_edge(path[i], path[i + 1]));
            }
            // Add closing edge (u, v)
            cycle_edges.push_back(make_edge(u, v));
            
            return cycle_edges;
        }
        
        // Explore neighbors (children in directed tree)
        for (int neighbor : tree[curr]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                parent[neighbor] = curr;
                q.push(neighbor);
            }
        }
        
        // Also explore parent direction (who points to curr?)
        for (int other = 0; other < V; other++) {
            if (!visited[other]) {
                for (int child : tree[other]) {
                    if (child == curr) {
                        visited[other] = true;
                        parent[other] = curr;
                        q.push(other);
                        break;
                    }
                }
            }
        }
    }
    
    // Should not reach here if tree is connected and edge is valid
    std::ostringstream msg;
    msg << "Could not find path from vertex " << u << " to vertex " << v << " in tree. ";
    msg << "Tree has " << V << " vertices. ";
    msg << "This suggests the tree is disconnected or the edge is invalid.";
    Rcpp::stop(msg.str());
    return std::vector<Edge>();
}

/*
 * Connected components via DFS
 */

void dfs_component(const Tree &tree, int u, std::vector<bool> &visited,
                  std::vector<int> &component, const MarkedEdgeSet &marked_edges) {
    visited[u] = true;
    component.push_back(u);
    
    // Tree is directed (edges point away from root), but we need to traverse it as undirected
    // Follow edges to children
    for (int v : tree[u]) {
        Edge e = make_edge(u, v);
        // Only traverse if edge is not marked and neighbor not visited
        if (marked_edges.find(e) == marked_edges.end() && !visited[v]) {
            dfs_component(tree, v, visited, component, marked_edges);
        }
    }
    
    // Also follow edges from parents (find who points to u)
    for (size_t other = 0; other < tree.size(); other++) {
        if ((int)other != u && !visited[other]) {
            for (int child : tree[other]) {
                if (child == u) {
                    Edge e = make_edge(other, u);
                    if (marked_edges.find(e) == marked_edges.end()) {
                        dfs_component(tree, other, visited, component, marked_edges);
                    }
                    break;
                }
            }
        }
    }
}

std::vector<std::vector<int>> tree_components_list(const Tree &tree,
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

uvec tree_to_partition(const Tree &tree, const MarkedEdgeSet &marked_edges,
                       int V, int n_distr) {
    auto components = tree_components_list(tree, marked_edges);
    
    // Verify we got the expected number of components
    if ((int)components.size() != n_distr) {
        std::ostringstream msg;
        msg << "tree_to_partition: Expected " << n_distr << " components but got " 
            << components.size() << ". Marked edges: " << marked_edges.size();
        Rcpp::warning(msg.str());
    }
    
    uvec partition(V);
    for (size_t i = 0; i < components.size(); i++) {
        for (int u : components[i]) {
            partition(u) = i + 1;  // 1-indexed for R
        }
    }
    
    return partition;
}

/*
 * Laplacian operations
 */

arma::sp_mat build_laplacian(const std::vector<Edge> &edges, int n_vertices) {
    if (n_vertices <= 0) {
        Rcpp::stop("Number of vertices must be positive");
    }
    
    // Count edges for memory allocation efficiency
    std::vector<int> degree(n_vertices, 0);
    
    // Build adjacency and degree information
    arma::sp_mat A(n_vertices, n_vertices);
    
    for (const Edge &e : edges) {
        int u = e.first;
        int v = e.second;
        
        if (u >= n_vertices || v >= n_vertices || u < 0 || v < 0) {
            std::ostringstream msg;
            msg << "Edge vertex out of bounds: edge (" << u << "," << v << ") ";
            msg << "but n_vertices=" << n_vertices;
            Rcpp::stop(msg.str());
        }
        
        // Add edge to adjacency matrix (with multiplicity for multigraphs)
        A(u, v) += 1.0;
        A(v, u) += 1.0;
        
        degree[u]++;
        degree[v]++;
    }
    
    // Build Laplacian L = D - A
    arma::sp_mat L(n_vertices, n_vertices);
    
    // Set diagonal (degree matrix)
    for (int i = 0; i < n_vertices; i++) {
        L(i, i) = degree[i];
    }
    
    // Subtract adjacency matrix
    L = L - A;
    
    return L;
}

double log_det_laplacian(const arma::sp_mat &L) {
    int n = L.n_rows;
    
    if (n <= 1) {
        return 0.0;  // log(1) = 0 for trivial case
    }
    
    // Extract reduced Laplacian (remove first row and column)
    // L_reduced = L[1:n-1, 1:n-1]
    arma::sp_mat L_red = L.submat(1, 1, n - 1, n - 1);
    
    // Convert to dense for Cholesky (reduced Laplacian should be small enough)
    arma::mat L_dense(L_red);
    
    // Attempt Cholesky decomposition
    // L_red should be positive definite for connected graph
    arma::mat U;
    bool success = arma::chol(U, L_dense);
    
    if (!success) {
        // Laplacian should be positive semi-definite
        // Reduced Laplacian should be positive definite if graph is connected
        // Return -Inf, which will cause proposal rejection in MH step
        return -std::numeric_limits<double>::infinity();
    }
    
    // log(det(L)) = log(det(U^T * U)) = log(det(U)^2) = 2 * log(det(U))
    // log(det(U)) = sum(log(diag(U))) for upper triangular U
    arma::vec diag_U = U.diag();
    
    // Check for numerical issues
    if (arma::any(diag_U <= 0)) {
        // Non-positive diagonal indicates numerical issues
        return -std::numeric_limits<double>::infinity();
    }
    
    double log_det = 2.0 * arma::sum(arma::log(diag_U));
    
    return log_det;
}

double log_spanning_trees(const std::vector<Edge> &edges, int n_vertices) {
    if (n_vertices <= 1) return 0.0;
    if (edges.empty()) return -std::numeric_limits<double>::infinity();
    
    arma::sp_mat L = build_laplacian(edges, n_vertices);
    return log_det_laplacian(L);
}

/*
 * Quotient graph construction
 */

std::vector<Edge> build_quotient_edges(const Graph &g, const uvec &partition,
                                       int n_distr) {
    std::vector<Edge> quotient_edges;
    
    for (size_t u = 0; u < g.size(); u++) {
        for (int v : g[u]) {
            // If edge crosses districts, add to quotient
            if (partition(u) != partition((size_t)v)) {
                Edge e = make_edge(partition(u) - 1, partition((size_t)v) - 1);  // 0-indexed
                quotient_edges.push_back(e);
            }
        }
    }
    
    return quotient_edges;
}

std::vector<std::vector<Edge>> district_induced_edges(const Graph &g,
                                                       const std::vector<std::vector<int>> &components) {
    std::vector<std::vector<Edge>> district_edges(components.size());
    
    // Create map from original vertex to (component_index, remapped_vertex)
    std::vector<int> vertex_to_comp(g.size());
    std::vector<int> vertex_remap(g.size());
    
    for (size_t i = 0; i < components.size(); i++) {
        for (size_t j = 0; j < components[i].size(); j++) {
            int u = components[i][j];
            vertex_to_comp[u] = i;
            vertex_remap[u] = j;  // Remap to 0, 1, 2, ...
        }
    }
    
    // For each edge in g, add to appropriate district if both endpoints in same district
    // Use remapped vertex indices
    for (size_t u = 0; u < g.size(); u++) {
        for (int v : g[u]) {
            if (vertex_to_comp[u] == vertex_to_comp[v] && (int)u < v) {
                Edge e = make_edge(vertex_remap[u], vertex_remap[v]);
                district_edges[vertex_to_comp[u]].push_back(e);
            }
        }
    }
    
    return district_edges;
}

/*
 * Calculate spanning tree ratios (key MEW component)
 */

double calculate_taus(const Graph &g, const Tree &tree_old, const Tree &tree_new,
                     const MarkedEdgeSet &marked_old, const MarkedEdgeSet &marked_new) {
    int V = g.size();
    
    // Get connected components for old and new states
    auto comp_old = tree_components_list(tree_old, marked_old);
    auto comp_new = tree_components_list(tree_new, marked_new);
    
    // Get induced edges within each district for both states
    auto edges_old = district_induced_edges(g, comp_old);
    auto edges_new = district_induced_edges(g, comp_new);
    
    // Compute log spanning trees for each district
    double tau_spanning = 0.0;
    
    for (size_t i = 0; i < edges_old.size(); i++) {
        if (comp_old[i].size() > 1) {  // Only if district has > 1 vertex
            tau_spanning += log_spanning_trees(edges_old[i], comp_old[i].size());
        }
    }
    
    for (size_t i = 0; i < edges_new.size(); i++) {
        if (comp_new[i].size() > 1) {
            tau_spanning -= log_spanning_trees(edges_new[i], comp_new[i].size());
        }
    }
    
    // Build quotient (district-level) graphs
    // Convert component lists to partition vectors
    uvec partition_old(V);
    for (size_t i = 0; i < comp_old.size(); i++) {
        for (int u : comp_old[i]) {
            partition_old(u) = i + 1;  // 1-indexed
        }
    }
    
    uvec partition_new(V);
    for (size_t i = 0; i < comp_new.size(); i++) {
        for (int u : comp_new[i]) {
            partition_new(u) = i + 1;
        }
    }
    
    int n_distr = comp_old.size();
    auto quot_edges_old = build_quotient_edges(g, partition_old, n_distr);
    auto quot_edges_new = build_quotient_edges(g, partition_new, n_distr);
    
    // Compute log spanning trees for quotient graphs
    double tau_quotient = 0.0;
    
    if (n_distr > 1) {
        tau_quotient = log_spanning_trees(quot_edges_old, n_distr) -
                      log_spanning_trees(quot_edges_new, n_distr);
    }
    
    // Combined tau
    double tau = tau_quotient + tau_spanning;
    
    return tau;
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
                             const Tree &tree_old,
                             const Tree &tree_new) {
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
    double pm = 1.0;
    
    // Check if marked edge vertices are same (both endpoints of same edge)
    if ((u == w && v == x) || (u == x && v == w)) {
        // Both endpoints same - more complex degree calculation
        int d_u = tree_old[u].size();
        int d_u_p = tree_new[u].size();
        int d_v = tree_old[v].size();
        int d_v_p = tree_new[v].size();
        
        pm = ((double)(d_u + d_v) / (d_u_p + d_v_p)) * 
             ((double)d_u_p / d_u) * ((double)d_v_p / d_v);
    } else {
        // One endpoint shared
        int shared = -1;
        if (u == w || u == x) {
            shared = u;
        } else if (v == w || v == x) {
            shared = v;
        }
        
        if (shared >= 0) {
            int d_u = tree_old[shared].size();
            int d_u_p = tree_new[shared].size();
            pm = (double)d_u_p / d_u;
        }
    }
    
    // Compute pt: accounts for cycle/marked edge intersection
    std::set<Edge> cycle_set(cycle_edges.begin(), cycle_edges.end());
    
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
    
    return pt * pm;
}

/*
 * Proposal mechanisms
 */

CycleProposal cycle_basis_step(const Graph &g, const Tree &tree,
                               const MarkedEdgeSet &marked_edges) {
    // Find edges in g but not in tree
    std::vector<Edge> non_tree_edges;
    
    for (size_t u = 0; u < g.size(); u++) {
        for (int v : g[u]) {
            if ((int)u < v) {  // Only consider each edge once
                Edge e = make_edge(u, v);
                if (!has_tree_edge(tree, u, v)) {
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
    std::vector<Edge> cycle_edges = find_cycle(tree, edge_plus.first, edge_plus.second);
    
    // Find edges in cycle that are not marked
    std::vector<Edge> possible_cuts;
    for (const Edge &e : cycle_edges) {
        if (e != edge_plus && marked_edges.find(e) == marked_edges.end()) {
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
        proposal.tree_new = tree;  // No change
        proposal.valid = false;
        return proposal;
    }
    
    // Sample edge to remove
    int cut_idx = r_int(possible_cuts.size());
    Edge edge_minus = possible_cuts[cut_idx];
    
    // Create new tree
    Tree tree_new = tree;
    add_tree_edge(tree_new, edge_plus.first, edge_plus.second);
    remove_tree_edge(tree_new, edge_minus.first, edge_minus.second);
    
    CycleProposal proposal;
    proposal.cycle_edges = cycle_edges;
    proposal.edge_plus = edge_plus;
    proposal.edge_minus = edge_minus;
    proposal.tree_new = tree_new;
    proposal.valid = true;
    
    return proposal;
}

MarkedEdgeProposal marked_edge_step(const Tree &tree,
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
    int chosen_vertex = (r_int(2) == 0) ? old_edge.first : old_edge.second;
    
    // Get neighbors of chosen vertex in tree (both parents and children)
    // Tree is directed with edges pointing away from root
    std::vector<int> neighbors;
    for (int child : tree[chosen_vertex]) {
        neighbors.push_back(child);
    }
    // Check other vertices for parent relationship
    for (size_t u = 0; u < tree.size(); u++) {
        if ((int)u != chosen_vertex) {
            for (int child : tree[u]) {
                if (child == chosen_vertex) {
                    neighbors.push_back(u);
                    break;
                }
            }
        }
    }
    
    if (neighbors.empty()) {
        Rcpp::stop("Chosen vertex has no neighbors in tree");
    }
    
    // Sample one neighbor
    int neighbor_idx = r_int(neighbors.size());
    int neighbor = neighbors[neighbor_idx];
    
    // Create new edge
    Edge new_edge = make_edge(chosen_vertex, neighbor);
    
    // Check if new edge is already marked (would reduce set size)
    // If so, pick a different neighbor
    int max_attempts = neighbors.size();
    int attempts = 0;
    while (marked_edges.find(new_edge) != marked_edges.end() && attempts < max_attempts) {
        neighbor_idx = r_int(neighbors.size());
        neighbor = neighbors[neighbor_idx];
        new_edge = make_edge(chosen_vertex, neighbor);
        attempts++;
    }
    
    // If all neighbors are marked edges, just return old edge (no change)
    if (marked_edges.find(new_edge) != marked_edges.end()) {
        MarkedEdgeProposal proposal;
        proposal.old_edge = old_edge;
        proposal.new_edge = old_edge;  // No change
        proposal.marked_new = marked_edges;  // No change
        return proposal;
    }
    
    // Update marked edge set
    MarkedEdgeSet marked_new = marked_edges;
    marked_new.erase(old_edge);
    marked_new.insert(new_edge);
    
    // Verify size is maintained (k-1 marked edges)
    if (marked_new.size() != marked_edges.size()) {
        std::ostringstream msg;
        msg << "marked_edge_step: Set size changed from " << marked_edges.size() 
            << " to " << marked_new.size();
        Rcpp::stop(msg.str());
    }
    
    MarkedEdgeProposal proposal;
    proposal.old_edge = old_edge;
    proposal.new_edge = new_edge;
    proposal.marked_new = marked_new;
    
    return proposal;
}

MEWProposal mew_proposal(const Graph &g, const Tree &tree,
                        const MarkedEdgeSet &marked_edges,
                        const uvec &pop, int n_distr,
                        double target, double lower, double upper) {
    const int MAX_TRIES = 1000;
    int V = g.size();
    
    for (int tries = 0; tries < MAX_TRIES; tries++) {
        // Propose tree update
        CycleProposal cycle_prop = cycle_basis_step(g, tree, marked_edges);
        
        // If cycle proposal is invalid (all edges marked), skip this iteration
        if (!cycle_prop.valid) {
            continue;
        }
        
        // Propose marked edge update on OLD tree (critical for detailed balance)
        // Must use 'tree' not 'cycle_prop.tree_new' to maintain reversibility
        MarkedEdgeProposal marked_prop = marked_edge_step(tree, marked_edges);
        
        // Check if new marked edge is valid
        // Cannot be edge_plus (not yet in tree) or edge_minus (removed from tree)
        if (marked_prop.new_edge == cycle_prop.edge_plus || 
            marked_prop.new_edge == cycle_prop.edge_minus) {
            continue;  // Invalid proposal, try again
        }
        
        // Convert to partition and check population
        uvec partition = tree_to_partition(cycle_prop.tree_new, marked_prop.marked_new, V, n_distr);
        
        // Compute district populations
        std::vector<double> dist_pop(n_distr, 0.0);
        for (int i = 0; i < V; i++) {
            int dist = partition(i) - 1;  // Convert to 0-indexed
            if (dist >= 0 && dist < n_distr) {
                dist_pop[dist] += pop(i);
            }
        }
        
        // Check if all districts meet population constraints
        bool valid = true;
        for (int d = 0; d < n_distr; d++) {
            if (dist_pop[d] < lower || dist_pop[d] > upper) {
                valid = false;
                break;
            }
        }
        
        if (valid) {
            // Valid proposal found
            MEWProposal proposal;
            proposal.cycle = cycle_prop;
            proposal.marked = marked_prop;
            proposal.n_rejects = tries;
            proposal.valid = true;
            return proposal;
        }
    }
    
    // Max tries exceeded - return invalid proposal
    // This will be automatically rejected in MH step
    MEWProposal proposal;
    proposal.cycle.tree_new = tree;
    proposal.cycle.edge_plus = make_edge(0, 0);
    proposal.marked.marked_new = marked_edges;
    proposal.marked.old_edge = make_edge(0, 0);
    proposal.marked.new_edge = make_edge(0, 0);
    proposal.n_rejects = MAX_TRIES;
    proposal.valid = false;
    
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

// DEPRECATED: Old initialization function that doesn't work reliably
// Kept for reference but should not be used
MarkedEdgeSet find_marked_edges_from_plan(const Tree &tree, const uvec &plan,
                                         int n_distr) {
    MarkedEdgeSet marked_edges;
    std::vector<Edge> tree_edges = tree_to_edges(tree);
    
    // Try to find k-1 edges whose removal creates the desired partition
    // Strategy: test each combination of k-1 edges
    
    // For efficiency, first try individual edges if k=2
    if (n_distr == 2) {
        for (const Edge &e : tree_edges) {
            MarkedEdgeSet test_set;
            test_set.insert(e);
            
            uvec test_partition = tree_to_partition(tree, test_set, tree.size(), n_distr);
            
            // Check if partition matches (allowing for relabeling)
            bool matches = true;
            std::map<int, int> label_map;
            
            for (size_t i = 0; i < plan.size(); i++) {
                int plan_label = plan(i);
                int test_label = test_partition(i);
                
                if (label_map.count(plan_label) == 0) {
                    label_map[plan_label] = test_label;
                } else if (label_map[plan_label] != test_label) {
                    matches = false;
                    break;
                }
            }
            
            if (matches) {
                return test_set;
            }
        }
    }
    
    // For k > 2, need to find k-1 edges
    // Strategy: Find all edges that cross districts, then search for valid subset
    
    std::vector<Edge> boundary_edges;
    for (const Edge &e : tree_edges) {
        int u = e.first;
        int v = e.second;
        
        // If this edge crosses districts in the plan, it's a candidate
        if (plan(u) != plan(v)) {
            boundary_edges.push_back(e);
        }
    }
    
    // Try all combinations of k-1 boundary edges
    // For small k this is feasible
    if ((int)boundary_edges.size() >= n_distr - 1) {
        // Try all C(n, k-1) combinations
        std::vector<int> indices(n_distr - 1);
        for (int i = 0; i < n_distr - 1; i++) {
            indices[i] = i;
        }
        
        bool found = false;
        while (!found) {
            // Test this combination
            MarkedEdgeSet test_set;
            for (int idx : indices) {
                test_set.insert(boundary_edges[idx]);
            }
            
            uvec test_partition = tree_to_partition(tree, test_set, tree.size(), n_distr);
            
            // Check if partition matches (allowing for relabeling)
            bool matches = true;
            std::map<int, int> label_map;
            
            for (size_t i = 0; i < plan.size(); i++) {
                int plan_label = plan(i);
                int test_label = test_partition(i);
                
                if (label_map.count(plan_label) == 0) {
                    label_map[plan_label] = test_label;
                } else if (label_map[plan_label] != test_label) {
                    matches = false;
                    break;
                }
            }
            
            if (matches) {
                return test_set;
            }
            
            // Next combination (lexicographic order)
            int i = n_distr - 2;
            while (i >= 0 && indices[i] == (int)boundary_edges.size() - (n_distr - 1 - i)) {
                i--;
            }
            
            if (i < 0) {
                break;  // Tried all combinations
            }
            
            indices[i]++;
            for (int j = i + 1; j < n_distr - 1; j++) {
                indices[j] = indices[j - 1] + 1;
            }
        }
    }
    
    // Fallback: just pick first k-1 edges that cross districts
    marked_edges.clear();
    for (const Edge &e : boundary_edges) {
        marked_edges.insert(e);
        if ((int)marked_edges.size() >= n_distr - 1) {
            break;
        }
    }
    
    // Final fallback: pick any k-1 edges from tree
    if ((int)marked_edges.size() < n_distr - 1) {
        marked_edges.clear();
        int count = 0;
        for (const Edge &e : tree_edges) {
            marked_edges.insert(e);
            count++;
            if (count >= n_distr - 1) break;
        }
    }
    
    return marked_edges;
}

/*
 * Energy functions and constraint evaluation
 */

// Helper: Count cut edges between districts
int count_cut_edges(const Graph &g, const uvec &partition) {
    int cuts = 0;
    for (size_t u = 0; u < g.size(); u++) {
        for (int v : g[u]) {
            if ((int)u < v && partition(u) != partition(v)) {
                cuts++;
            }
        }
    }
    return cuts;
}

// Helper: Count administrative splits
int count_splits(const uvec &partition, const uvec &admin, int n_distr) {
    if (admin.n_elem == 0) return 0;
    
    int n_admin = max(admin);
    int splits = 0;
    
    // For each admin unit, count how many districts it touches
    for (int a = 1; a <= n_admin; a++) {
        std::set<int> dists_in_admin;
        for (size_t i = 0; i < admin.n_elem; i++) {
            if (admin(i) == a) {
                dists_in_admin.insert(partition(i));
            }
        }
        // Number of splits = (number of districts touching - 1)
        if (dists_in_admin.size() > 1) {
            splits += (dists_in_admin.size() - 1);
        }
    }
    
    return splits;
}

// Main constraint evaluation function
// Returns LOG of energy ratio: log(E_new / E_old)
double compute_energy_ratio(const Graph &g,
                           const Tree &tree_old,
                           const MarkedEdgeSet &marked_old,
                           const Tree &tree_new,
                           const MarkedEdgeSet &marked_new,
                           const uvec &pop,
                           double rho,
                           List constraints) {
    // Convert states to partitions
    int V = g.size();
    int n_distr = marked_old.size() + 1;
    
    uvec partition_old = tree_to_partition(tree_old, marked_old, V, n_distr);
    uvec partition_new = tree_to_partition(tree_new, marked_new, V, n_distr);
    
    double log_energy_ratio = 0.0;  // Sum of log energy differences
    
    // 1. Compactness constraint (edge cuts)
    if (rho != 0) {
        int cuts_old = count_cut_edges(g, partition_old);
        int cuts_new = count_cut_edges(g, partition_new);
        
        // Energy: exp(-rho * cuts)
        // Log ratio: log(E_new/E_old) = -rho * (cuts_new - cuts_old)
        log_energy_ratio += -rho * (cuts_new - cuts_old);
    }
    
    // 2. Process other constraints from list
    // Constraints format: list(splits = list(list(strength=..., admin=...)), ...)
    if (constraints.size() > 0) {
        CharacterVector constraint_names = constraints.names();
        
        for (int i = 0; i < constraints.size(); i++) {
            std::string name = Rcpp::as<std::string>(constraint_names[i]);
            List constraint_list = constraints[i];
            
            // Each constraint can have multiple instances
            for (int j = 0; j < constraint_list.size(); j++) {
                List params = constraint_list[j];
                double strength = params["strength"];
                
                if (name == "splits" && params.containsElementNamed("admin")) {
                    uvec admin = params["admin"];
                    int splits_old = count_splits(partition_old, admin, n_distr);
                    int splits_new = count_splits(partition_new, admin, n_distr);
                    
                    // Energy: exp(-strength * splits)
                    // Log ratio: -strength * (splits_new - splits_old)
                    log_energy_ratio += -strength * (splits_new - splits_old);
                    
                } else if (name == "edges_rem") {
                    // Already captured by cuts - same thing
                    int cuts_old = count_cut_edges(g, partition_old);
                    int cuts_new = count_cut_edges(g, partition_new);
                    log_energy_ratio += -strength * (cuts_new - cuts_old);
                }
                // More constraints can be added here as needed
                // For now, focus on compactness and splits as they are most common
            }
        }
    }
    
    return log_energy_ratio;
}
