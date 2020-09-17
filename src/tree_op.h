#include "smc_base.h"

#ifndef TREE_OP_H
#define TREE_OP_H

/*
 * Generate a random vertex (integer) among unvisited vertices
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining);

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
 * Convert R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const List &l);

/*
 * Count population below each node in tree
 */
// TESTED
double tree_pop(Tree &ust, int vtx, const uvec &pop,
                std::vector<int> &pop_below, std::vector<int> &parent);

/*
 * Assign `district` to all descendants of `root` in `ust`
 */
// TESTED
void assign_district(const Tree &ust, subview_col<uword> &districts,
                     int root, int district);

/*
 * Find the root of a subtree.
 */
// TESTED
int find_subroot(const Tree &ust, const std::vector<bool> &ignore);

#endif
