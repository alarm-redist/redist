#ifndef TREE_OP_H
#define TREE_OP_H

#include <vector>
#include <limits>
#include <stack>
#include <queue>
#include <RcppArmadillo.h>
#include "redist_types.h"
#include "smc_base.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;



/*
 * Generate a random vertex (integer) among unvisited vertices
 * `lower` is a lower bound (inclusive) on the index of the first unvisited element
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining, int &lower,
         RNGState &rng_state);

/*
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
int rnbor(const Graph &g, int vtx, RNGState &rng_state);

/*
 * Make the district adjacency graph for `plan` from the overall precinct graph `g`
 * if `zero`=false then ignore zeros, otherwise map them to `nd`
 */
// TESTED
Graph district_graph(const Graph &g, PlanVector const &region_ids, int nd, bool zero=false);




/*
 * Initialize empty tree structure on graph with `V` vertices
 */
// TESTED
Tree init_tree(int V);

/*
 * Clear a tree
 */
void clear_tree(Tree &tree);

// print a tree
void print_tree(Tree const &ust);


/*
 * Count population below each node in tree and get parent
 */
// TESTED
int tree_pop(Tree &ust, int vtx, const arma::uvec &pop,
             std::vector<int> &pop_below, std::vector<int> &parent);


/*
 * Just Count population below each node in tree
 */
// TESTED
int get_tree_pops_below(const Tree &ust, const int vtx, const arma::uvec &pop,
             std::vector<int> &pop_below);




/*
 * Assign `district` to all descendants of `root` in `ust`
 */
// TESTED
void assign_district(const Tree &ust, arma::subview_col<arma::uword> &districts,
                     int root, int district);



/*
 * Assign `new_region_num_id` to all descendants of `root` in `ust`
 */
// NOT TESTED BUT COPIED FROM `assign_district` so should be ok
void assign_region_id_from_tree(const Tree &ust, 
                    PlanVector &region_ids,
                    int root,
                    const int new_region_num_id);


void assign_region_id_and_forest_from_tree(const Tree &ust, 
                    PlanVector &region_ids,
                    VertexGraph &forest_graph,
                    int root,
                    const int new_region_id);

/*
 * Find the root of a subtree.
 */
// TESTED
int find_subroot(const Tree &ust, const std::vector<bool> &ignore);


/*  
 *  Erases an edge from a tree
 * 
 *  Erases the directed edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`)
 *  from the tree `ust`. The directed edge here means we have `child_vertex` being one of 
 *  the values in `ust[parent_vertex]`.
 * 
 * 
 *  @param ust A directed spanning tree passed by reference
 *  @param cut_edge An `EdgeCut` object representing the edge cut
 * 
 *  @details Modifications
 *     - The edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`) 
 *     is removed from `ust`
 * 
 * 
 */ 
void erase_tree_edge(Tree &ust, EdgeCut cut_edge);

#endif
