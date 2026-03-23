/********************************************************
 * BUD MCMC Redistricting Sampler
 * BUD Partition State Management
 *
 * Extends LCTPartition with district-level tree, marked edges,
 * cuttable tree info, and temporary cycle LCT.
 ********************************************************/

#ifndef BUD_PARTITION_H
#define BUD_PARTITION_H

#include "cw_partition.h"
#include "cw_lct.h"
#include "cw_forest_walk.h"
#include <map>
#include <set>
#include <vector>
#include <queue>

/*
 * Population interval for cuttable tree closures.
 * Represents a range [lo, hi] of possible population for a subtree configuration.
 */
struct PopInterval {
    double lo, hi;
    PopInterval(double l = 0, double h = 0) : lo(l), hi(h) {}
    bool operator<(const PopInterval& o) const { return lo < o.lo; }
};

/*
 * Cuttable info for one node: maps number_of_cuts -> list of population intervals
 */
using NodeClosure = std::map<int, std::vector<PopInterval>>;

/*
 * CuttableInfo: Stores cuttable tree information for efficient balance checking.
 */
struct CuttableInfo {
    int root;
    int identifier;
    std::vector<NodeClosure> node_closures;  // one per vertex

    CuttableInfo() : root(-1), identifier(-1) {}
    CuttableInfo(int n) : root(-1), identifier(-1), node_closures(n) {}
};

/*
 * Initialize a node closure for a leaf with given population value.
 */
NodeClosure init_node_closure(double val, double lower, double upper);

/*
 * Merge two closures (convolution of population intervals).
 */
void merge_closures(NodeClosure& c1, const NodeClosure& c2, double lower, double upper);

/*
 * Add cut possibilities to a closure.
 */
void add_cuts_to_closures(NodeClosure& closure, double lower, double upper);

/*
 * Clean/merge overlapping intervals in a closure vector.
 */
std::vector<PopInterval> clean_closure(std::vector<PopInterval>& closure, double epsilon);

/*
 * Check if a tree with given root closure can be cut into num_dists balanced pieces.
 */
bool cuttable_huh(const NodeClosure& root_closure, int num_dists);

/*
 * Compute cuttable info for entire tree rooted at root_node.
 */
void set_cuttable_tree_info(LinkCutTree& lct, int root_vertex,
                            const std::vector<int>& pop_values,
                            double lower, double upper,
                            CuttableInfo& info);

/*
 * Update cuttable info along a path (efficient re-rooting).
 */
void set_cuttable_path_info(LinkCutTree& lct,
                            const std::vector<int>& path,
                            const std::vector<int>& pop_values,
                            double lower, double upper,
                            CuttableInfo& info);

/*
 * BUDPartition: Extends LCTPartition with BUD-specific state.
 */
class BUDPartition : public LCTPartition {
public:
    // District-level tree (adjacency list, 0-indexed districts)
    std::vector<std::set<int>> district_tree;

    // Marked edges: for each pair of adjacent districts, one node-level edge
    // Key: ordered pair (min_d, max_d), Value: (node_u, node_v)
    std::map<DistrictPair, std::pair<int, int>> marked_edges;

    // Cuttable tree info for the main spanning forest
    CuttableInfo cuttable_info;

    // Temporary LCT for collapsed cycle operations
    LinkCutTree tmp_cycle_lct;

    // Population values (as int vector for fast access)
    std::vector<int> pop_values;

    // Population bounds
    double pop_lower, pop_upper;

    // Forest edges: all edges in the spanning forest (for LCT rebuild without Wilson's)
    std::vector<std::pair<int,int>> forest_edges;

    BUDPartition(int n_vertices, int n_districts);

    /*
     * Initialize from plan - sets up spanning trees, district tree, marked edges.
     */
    int init_from_plan(const Graph& g,
                       const arma::uvec& plan,
                       const arma::uvec& population,
                       const arma::uvec& county_assignments,
                       double lower, double upper);

    /*
     * Sample a spanning tree on the district graph using Wilson's algorithm.
     * Returns edges of the district-level tree.
     */
    void sample_district_tree();

    /*
     * For each edge in the district tree, sample a marked edge from cross-district edges.
     */
    void sample_marked_edges();

    /*
     * Get the district-level path between two districts using BFS.
     * Returns list of district pairs along the path.
     */
    std::vector<std::pair<int,int>> get_district_path(int d1, int d2) const;

    /*
     * Rebuild district tree adjacency from marked_edges.
     */
    void rebuild_district_tree();

    /*
     * Rebuild LCT from forest_edges (deterministic, preserves configuration).
     */
    void rebuild_lct_from_forest();

    /*
     * Extract forest edges from the current LCT state.
     */
    void extract_forest_edges();

    /*
     * Update cross-district edges for specific affected districts.
     */
    void update_cross_district_edges(const Graph& g,
                                     const std::vector<int>& affected_districts);

    /*
     * Full recomputation of cross-district edges (public wrapper).
     */
    void recompute_cross_district_edges() { find_cross_district_edges(); }

    /*
     * Print BUD-specific state for debugging.
     */
    void print_bud_state(int verbosity = 1) const;
};

#endif // BUD_PARTITION_H
