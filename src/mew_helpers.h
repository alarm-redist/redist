#ifndef MEW_HELPERS_H
#define MEW_HELPERS_H

#include "cw_lct.h"
#include "smc_base.h"
#include <cstddef>
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

struct MEWTreeUpdate {
    Edge edge_plus;
    Edge edge_minus;
    std::size_t minus_pos_first;
    std::size_t minus_pos_second;
    bool changed;
};

class MEWTree {
public:
    explicit MEWTree(Tree tree);
    MEWTree(const MEWTree&) = delete;
    MEWTree& operator=(const MEWTree&) = delete;

    int size() const;
    const std::vector<int>& neighbors(int vertex) const;
    int degree(int vertex) const;
    bool has_edge(int u, int v) const;
    std::vector<Edge> find_cycle(int u, int v);
    MEWTreeUpdate replace_edge(const Edge& edge_plus, const Edge& edge_minus);
    void rollback(const MEWTreeUpdate& update);

private:
    LinkCutTree lct_;
    Tree adjacency_;
    std::set<Edge> edges_;

    static std::size_t remove_neighbor(std::vector<int>& neighbors, int vertex);
};

/*
 * Proposal structures
 */

struct CycleProposal {
    std::vector<Edge> cycle_edges;  // Edges forming the cycle
    Edge edge_plus;                 // Edge added to tree
    Edge edge_minus;                // Edge removed from tree
    MEWTreeUpdate update;           // Applied LCT/adjacency update
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
    int n_rejects;                  // Number of internal redraws
    bool valid;                     // Whether proposal meets population constraints
    uvec partition;                 // Cached partition from population check
};

/*
 * Connected components
 */

// Find connected components after removing marked edges from tree
// Returns partition assignment (1-indexed district labels)
uvec tree_to_partition(const MEWTree &tree, const MarkedEdgeSet &marked_edges,
                       int V, int n_distr);

// Alternative: return list of vertex sets for each component
std::vector<std::vector<int>> tree_components_list(const MEWTree &tree,
                                                    const MarkedEdgeSet &marked_edges);

/*
 * MEW-specific computations
 */

// Compute forward/backward transition probability ratio
double transition_probability(const std::vector<Edge> &cycle_edges,
                             const Edge &edge_plus,
                             const Edge &marked_old,
                             const Edge &marked_new,
                             const MarkedEdgeSet &marked_edges_old,
                             const MarkedEdgeSet &marked_edges_new,
                             const Edge &edge_minus,
                             const MEWTree &tree_new);

// Compute log spanning tree count of the quotient (inter-district) graph
double log_st_quotient_graph(const Graph &g, const arma::uvec &plan, int n_distr);

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
CycleProposal cycle_basis_step(const Graph &g, MEWTree &tree,
                               const MarkedEdgeSet &marked_edges);

// Marked edge step: update marked edges
MarkedEdgeProposal marked_edge_step(const MEWTree &tree,
                                   const MarkedEdgeSet &marked_edges);

// Combined proposal with population constraint checking
MEWProposal mew_proposal(const Graph &g, MEWTree &tree,
                        const MarkedEdgeSet &marked_edges,
                        const uvec &pop, int n_distr,
                        double target, double lower, double upper);


#endif
