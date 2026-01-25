/********************************************************
 * CycleWalk MCMC Redistricting Sampler
 * LCT Partition State Management
 *
 * Maintains the state of the redistricting plan as a collection
 * of spanning trees (one per district) stored in a Link-Cut Tree.
 ********************************************************/

#ifndef CW_PARTITION_H
#define CW_PARTITION_H

#include "cw_lct.h"
#include "smc_base.h"
#include "tree_op.h"
#include "wilson.h"
#include <map>
#include <set>
#include <utility>

/*
 * Edge between two vertices, with weight.
 * Sorted so that first <= second.
 */
struct CWEdge {
    int u, v;
    double weight;

    CWEdge(int a, int b, double w = 1.0);

    bool operator<(const CWEdge& other) const;
    bool operator==(const CWEdge& other) const;
};

/*
 * Key for cross-district edges: ordered pair of districts.
 */
using DistrictPair = std::pair<int, int>;

/*
 * Set of edges between two districts.
 */
using EdgeSet = std::set<CWEdge>;

/*
 * Map from district pairs to edge sets.
 */
using CrossEdgeMap = std::map<DistrictPair, EdgeSet>;

/*
 * LCTPartition: Maintains the redistricting plan state.
 *
 * Each district is represented by a spanning tree stored in the LCT.
 * The partition tracks:
 * - District assignments for each vertex
 * - District roots (one per district)
 * - Cross-district boundary edges (for cycle walk proposals)
 * - Population totals per district
 */
class LCTPartition {
public:
    int n_districts;                      // Number of districts
    int n_vertices;                       // Number of vertices

    LinkCutTree lct;                      // Link-Cut Tree for all vertices
    std::vector<int> district_roots;      // Root vertex of each district's tree
    std::vector<int> node_to_district;    // District assignment for each vertex
    std::vector<int> district_pop;        // Population of each district
    CrossEdgeMap cross_edges;             // Edges between districts

    // References to input data (not owned)
    const Graph* graph;
    const arma::uvec* pop;
    const arma::uvec* counties;

    // Constructor: empty partition (use init_from_plan to populate)
    LCTPartition(int n_vertices, int n_districts);

    /*
     * Initialize the partition from an existing plan.
     * For each district:
     *   1. Use Wilson's algorithm to sample a spanning tree
     *   2. Load tree edges into the LCT
     *   3. Track the district root
     * Then find all cross-district edges.
     */
    int init_from_plan(const Graph& g,
                       const arma::uvec& plan,
                       const arma::uvec& population,
                       const arma::uvec& county_assignments,
                       double lower, double upper);

    /*
     * Get the district of a vertex.
     */
    int get_district(int v) const;

    /*
     * Get the population of a district.
     */
    int get_district_pop(int d) const;

    /*
     * Get all edges between two districts.
     */
    const EdgeSet& get_cross_edges(int d1, int d2) const;

    /*
     * Check if two districts are adjacent (share at least one boundary edge).
     */
    bool districts_adjacent(int d1, int d2) const;

    /*
     * Get all pairs of adjacent districts.
     */
    std::vector<DistrictPair> get_adjacent_district_pairs() const;

    /*
     * Get the current plan as a vector of district assignments.
     */
    arma::uvec get_plan() const;

    /*
     * Debug: Print partition state.
     */
    void print_state(int verbosity = 1) const;

private:
    /*
     * Load a spanning tree into the LCT.
     * tree: adjacency list representation (tree[v] = children of v)
     * root: root vertex of the tree
     * district: district index to assign
     */
    void load_tree_into_lct(const Tree& tree, int root, int district);

    /*
     * Find all cross-district edges by scanning the graph.
     */
    void find_cross_district_edges();

    /*
     * Assign district labels by traversing from root.
     */
    void assign_districts_from_root(int root, int district);

    // Empty edge set for returning when no edges exist
    static const EdgeSet empty_edge_set;
};

#endif // CW_PARTITION_H
