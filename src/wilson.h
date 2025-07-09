#ifndef WILSON_H
#define WILSON_H

#include "tree_op.h"

using namespace Rcpp;
using namespace arma;

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
int sample_sub_ust(
    MapParams const &map_params, Tree &tree, int &root,
    double const lower, double const upper,
    std::vector<bool> &visited, const std::vector<bool> &ignore, 
    Tree &cty_tree, arma::uvec &county_pop, std::vector<std::vector<int>> &county_members,
    std::vector<bool> &c_visited, std::vector<int> &cty_pop_below,
    std::vector<std::array<int, 3>> &county_path, std::vector<int> &path,
    RNGState &rng_state
);

#endif