#include "tree_op.h"

#ifndef WILSON_H
#define WILSON_H

/*
 * Sample a uniform spanning tree using Wilson's algorithm
 */
Tree sample_ust(List g, int &root, const uvec &counties);

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
Tree sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                    const std::vector<bool> &ignore,
                    const uvec &counties, Multigraph &mg);

#endif
