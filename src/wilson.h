#include "tree_op.h"

#ifndef WILSON_H
#define WILSON_H

#include <array>
#include <cstdint>

using LevelEdge = std::array<int, 3>;
using LevelPath = std::vector<LevelEdge>;
using ActiveMultigraph = std::vector<std::vector<LevelEdge>>;

struct HierarchicalSamplerWorkspace {
    std::vector<ActiveMultigraph> level_graphs;
    std::vector<Graph> internal_level_graphs;
    std::vector<std::vector<std::vector<LevelEdge>>> static_edges_by_vertex;
    std::vector<std::vector<int>> level_group_sizes;
    std::vector<std::vector<int>> level_local_positions;
    std::vector<std::vector<std::vector<int>>> level_group_vertices;
    std::vector<std::vector<std::uint64_t>> internal_neighbor_masks;
    std::vector<std::vector<std::vector<int>>> children_by_parent;
    std::vector<std::vector<std::vector<int>>> static_children_by_parent;
    std::vector<std::vector<int>> active_levels;
    std::vector<int> active_group_counts;
    std::vector<std::vector<int>> active_parents;
    std::vector<int> active_hierarchy_vertices;
    std::vector<int> active_graph_vertices;
    std::vector<int> component_stack;
    std::vector<int> component_parent;
    std::vector<int> component_labels;
    std::vector<std::uint64_t> active_unit_masks;
    std::vector<bool> group_active;
    std::vector<bool> group_has_root;
    std::vector<bool> group_visited;
    std::vector<bool> frontier_queued;
    std::vector<int> frontier_vertices;
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
                        bool relabel_active_components,
                        HierarchicalSamplerWorkspace &workspace,
                        const std::vector<int> &finest_adj,
                        const std::vector<int> &finest_off);

#endif
