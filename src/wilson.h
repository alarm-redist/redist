#ifndef WILSON_H
#define WILSON_H

#include "tree_op.h"

using namespace Rcpp;
using namespace arma;

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
int sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                   std::vector<bool> &visited,
                   const std::vector<bool> &ignore, const uvec &pop,
                   double lower, double upper,
                   const uvec &counties, const Multigraph &mg,
                   RNGState &rng_state);

#endif