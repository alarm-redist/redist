#include "tree_op.h"

#ifndef WILSON_H
#define WILSON_H

#include <array>
#include <cstdint>

using LevelEdge = std::array<int, 3>;
using LevelPath = std::vector<LevelEdge>;
using ActiveMultigraph = std::vector<std::vector<LevelEdge>>;

struct HierarchicalSamplerWorkspace {
    ActiveMultigraph active_graph;
    std::vector<bool> group_active;
    std::vector<bool> group_has_root;
    std::vector<bool> group_visited;
    std::vector<bool> parent_seen;
    std::vector<int> next_group;
    std::vector<int> next_edge;
    LevelPath level_path;
    Tree prune_group_tree;
    uvec prune_group_pop;
    std::vector<std::vector<int>> prune_members;
    std::vector<int> prune_group_pop_below;
    std::vector<int> prune_group_parent;
    std::vector<int> walk_path;
    std::vector<std::int8_t> status;
    std::vector<int> path_pos;
};

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
                        const std::vector<std::vector<int>> &levels,
                        const std::vector<int> &n_groups,
                        const std::vector<std::vector<int>> &parents,
                        HierarchicalSamplerWorkspace &workspace,
                        const std::vector<int> &finest_adj,
                        const std::vector<int> &finest_off);

#endif
