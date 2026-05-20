#include "tree_op.h"

#ifndef WILSON_H
#define WILSON_H

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
int sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                   std::vector<bool> &visited,
                   const std::vector<bool> &ignore, const uvec &pop,
                   double lower, double upper,
                   const uvec &counties, Multigraph &mg);

/*
 * Sample a hierarchical spanning subtree using multiple nested administrative
 * levels ordered from finest to coarsest.
 */
int sample_sub_ust_hier(const Graph &g, Tree &tree, int V, int &root,
                        std::vector<bool> &visited,
                        const std::vector<bool> &ignore,
                        const std::vector<int> &active_vertices,
                        const uvec &pop,
                        double lower, double upper,
                        const std::vector<uvec> &levels);

#endif
