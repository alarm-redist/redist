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

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
std::vector<int> walk_until(const Graph &g, int root,
                            const std::vector<bool> &visited,
                            const std::vector<bool> &ignore,
                            const uvec &counties);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
void loop_erase(std::vector<int> &path, int proposal);

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
std::vector<std::vector<int>> walk_until_cty(Multigraph &mg, int root,
                                             const std::vector<bool> &visited,
                                             const std::vector<bool> &ignore);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::vector<int>> &path, int proposal, int root);

#endif
