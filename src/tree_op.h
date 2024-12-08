#include "smc_base.h"

#ifndef TREE_OP_H
#define TREE_OP_H

/*
 * Generate a random vertex (integer) among unvisited vertices
 * `lower` is a lower bound (inclusive) on the index of the first unvisited element
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining, int &lower);

/*
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
int rnbor(const Graph &g, int vtx);

/*
 * Make a county graph from a precinct graph and list of counties
 */
// TESTED
Multigraph county_graph(const Graph &g, const uvec &counties);

/*
 * Make the district adjacency graph for `plan` from the overall precinct graph `g`
 * if `zero`=false then ignore zeros, otherwise map them to `nd`
 */
// TESTED
Graph district_graph(const Graph &g, const uvec &plan, int nd, bool zero=false);


// NOT FULLY TESTED but taken from district_graph function which was tested
//' Creates the region level graph of a plan
//'
//' Given a plan object this returns a graph of the regions in the plan using
//' the region ids as indices
//'
//' @title Get Region-Level Graph
//'
//' @param g The graph of the entire map
//' @param plan A plan object
//'
//' @details No modifications to inputs made
//'
//' @return the log of the probability the specific value of `region_to_split` was chosen
//'
Graph get_region_graph(const Graph &g, const Plan &plan);


/*
 * Initialize empty multigraph structure on graph with `V` vertices
 */
// TESTED
Multigraph init_multigraph(int V);

/*
 * Initialize empty tree structure on graph with `V` vertices
 */
// TESTED
Tree init_tree(int V);

/*
 * Clear a tree
 */
void clear_tree(Tree &tree);

/*
 * Convert R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const List &l);

/*
 * Count population below each node in tree
 */
// TESTED
int tree_pop(Tree &ust, int vtx, const uvec &pop,
             std::vector<int> &pop_below, std::vector<int> &parent);

/*
 * Assign `district` to all descendants of `root` in `ust`
 */
// TESTED
void assign_district(const Tree &ust, subview_col<uword> &districts,
                     int root, int district);

/*
 * Assign `new_region` to all descendants of `root` in `ust`
 */
// TESTED
void assign_region(const Tree &ust, Plan &plan,
                   int root,
                   int new_region_num_id);

// TEMP REMOVE LATER
void assign_region_just_vertex_vec(const Tree &ust, 
                    arma::subview_col<arma::uword> &region_ids,
                    int root,
                    int new_region_num_id);

/*
 * Find the root of a subtree.
 */
// TESTED
int find_subroot(const Tree &ust, const std::vector<bool> &ignore);

#endif
