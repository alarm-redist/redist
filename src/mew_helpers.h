#ifndef MEW_HELPERS_H
#define MEW_HELPERS_H

#include "smc_base.h"
#include <set>
#include <utility>
#include <vector>

/*
 * Data structures for MEW algorithm
 */

// Edge type: normalized pair (u, v) where u < v
typedef std::pair<int, int> Edge;

// Set of marked edges
typedef std::set<Edge> MarkedEdgeSet;

// Helper to create normalized edge
inline Edge make_edge(int u, int v) {
    return (u < v) ? Edge(u, v) : Edge(v, u);
}

/*
 * Proposal structures
 */

struct CycleProposal {
    std::vector<Edge> cycle_edges;  // Edges forming the cycle
    Edge edge_plus;                 // Edge added to tree
    Edge edge_minus;                // Edge removed from tree
    Tree tree_new;                  // New tree after cycle operation
    bool valid;                     // Whether a valid proposal was created (false if all cycle edges marked)
};

struct MarkedEdgeProposal {
    Edge old_edge;                  // Old marked edge
    Edge new_edge;                  // New marked edge
    MarkedEdgeSet marked_new;       // New set of marked edges
};

struct MEWProposal {
    CycleProposal cycle;            // Tree update
    MarkedEdgeProposal marked;      // Marked edge update
    int n_rejects;                  // Number of rejections before valid proposal
    bool valid;                     // Whether proposal meets population constraints
};

/*
 * Core graph/tree operations
 */

// Add edge to tree
void add_tree_edge(Tree &tree, int u, int v);

// Remove edge from tree
void remove_tree_edge(Tree &tree, int u, int v);

// Check if edge exists in tree
bool has_tree_edge(const Tree &tree, int u, int v);

// Convert tree to edge list
std::vector<Edge> tree_to_edges(const Tree &tree);

/*
 * Cycle detection
 */

// Find cycle formed by adding edge (u,v) to tree
// Uses BFS to find path from u to v in tree, then adds edge (u,v)
std::vector<Edge> find_cycle(const Tree &tree, int u, int v);

/*
 * Connected components
 */

// Find connected components after removing marked edges from tree
// Returns partition assignment (0-indexed district labels)
uvec tree_to_partition(const Tree &tree, const MarkedEdgeSet &marked_edges,
                       int V, int n_distr);

// Alternative: return list of vertex sets for each component
std::vector<std::vector<int>> tree_components_list(const Tree &tree,
                                                    const MarkedEdgeSet &marked_edges);

/*
 * Laplacian and spanning tree counting
 */

// Build sparse Laplacian matrix from edge list
arma::sp_mat build_laplacian(const std::vector<Edge> &edges, int n_vertices);

// Compute log(det(L_reduced)) for spanning tree count
// Uses Cholesky decomposition
double log_det_laplacian(const arma::sp_mat &L);

// Compute log number of spanning trees from edge list
double log_spanning_trees(const std::vector<Edge> &edges, int n_vertices);

/*
 * Quotient graph construction
 */

// Build quotient (district-level) edges from graph and partition
std::vector<Edge> build_quotient_edges(const Graph &g, const uvec &partition,
                                       int n_distr);

// Get edges within each district (induced subgraphs)
std::vector<std::vector<Edge>> district_induced_edges(const Graph &g,
                                                       const std::vector<std::vector<int>> &components);

/*
 * MEW-specific computations
 */

// Compute log-ratio of spanning tree counts (key MEW acceptance term)
double calculate_taus(const Graph &g, const Tree &tree_old, const Tree &tree_new,
                     const MarkedEdgeSet &marked_old, const MarkedEdgeSet &marked_new);

// Compute forward/backward transition probability ratio
double transition_probability(const std::vector<Edge> &cycle_edges,
                             const Edge &edge_plus,
                             const Edge &marked_old,
                             const Edge &marked_new,
                             const MarkedEdgeSet &marked_edges_old,
                             const MarkedEdgeSet &marked_edges_new,
                             const Tree &tree_old,
                             const Tree &tree_new);

/*
 * Initialization
 */

// Build tree and marked edges from a partition
// This is the key initialization function that guarantees the tree structure
// matches the partition
std::pair<Tree, MarkedEdgeSet> partition_to_tree_marked_edges(
    const Graph &g,
    const uvec &partition,
    int n_distr
);

/*
 * Proposal mechanisms
 */

// Cycle basis step: update spanning tree
CycleProposal cycle_basis_step(const Graph &g, const Tree &tree,
                               const MarkedEdgeSet &marked_edges);

// Marked edge step: update marked edges
MarkedEdgeProposal marked_edge_step(const Tree &tree,
                                   const MarkedEdgeSet &marked_edges);

// Combined proposal with population constraint checking
MEWProposal mew_proposal(const Graph &g, const Tree &tree,
                        const MarkedEdgeSet &marked_edges,
                        const uvec &pop, int n_distr,
                        double target, double lower, double upper);

/*
 * Initialization
 */

// Find marked edges from a partition
MarkedEdgeSet find_marked_edges_from_plan(const Tree &tree, const uvec &plan,
                                         int n_distr);

/*
 * Energy/constraint functions
 */

// Compute energy ratio for acceptance
double compute_energy_ratio(const Graph &g,
                           const Tree &tree_old,
                           const MarkedEdgeSet &marked_old,
                           const Tree &tree_new,
                           const MarkedEdgeSet &marked_new,
                           const uvec &pop,
                           double rho,
                           List constraints);

#endif
