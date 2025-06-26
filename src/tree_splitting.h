#pragma once
#ifndef TREE_SPLITTING_H
#define TREE_SPLITTING_H

#include <RcppArmadillo.h>
#include <vector>
#include <utility>
#include <queue>

#include "gredist_constants.h"
#include "redist_types.h"
#include "smc_base.h"
#include "tree_op.h"



// [[Rcpp::depends(RcppArmadillo)]]



/* 
 * Pick a valid tree edges to split with probability ∝ exp(-alpha*larger abs dev)
 *
 * 
 * Returns a valid tree edge to split with probability proporitional to 
 * exp(-alpha*larger dev) where larger abs dev is the bigger absolute deviation
 * from the target of the two regions induced by the cut. If successful returns
 * information on the edge and region sizes associated with the cut. 
 *
 * Note even if a successful cut is found it does not
 * update the plan or the tree.
 *
 *
 * It will only attempt to create regions where the size is between
 * min_potential_d and max_potential_d (inclusive). So the one district
 * split case is `min_potential_d=max_potential_d=1`.
 * 
 * Valid edge here is defined as an edge and region sizes such that the 
 * two induced regions both fall within the population bounds.
 *
 * 
 * @param root The root vertex of the spanning tree
 * @param pop_below The population corresponding to cutting below each vertex. 
 * So `pop_below[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
 * means the vertex is the root or it is not in the tree.
 * @param alpha Used in the exp() term. A larger alpha puts more weight on smaller
 * deviations and smaller makes the weight closer to uniform.
 * @param min_potential_cut_size The smallest potential region size to try for a cut. 
 * @param max_potential_cut_size The largest potential region size it will try for a cut. 
 * Setting this to 1 will result in only 1 district splits. 
 * @param region_ids A vector mapping 0 indexed vertices to their region id number
 * @param region_id_to_split The id of the region in the plan object we're attempting to split
 * @param total_region_pop The total population of the region being split 
 * @param total_region_size The size of the region being split 
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details No modifications made
 *
 * @return <True, information on the edge cut> if two valid regions were 
 * successfully split, false otherwise
 *
 */



// Gets the deviance for each edge in a tree
std::vector<double> get_ordered_tree_cut_devs(Tree &ust, int root,
                             std::vector<int> const &cut_below_pop, double const target,
                             PlanVector const &region_ids,
                             int const region_id, int const region_size, int const region_pop,
                             int const min_potential_cut_size, int const max_potential_cut_size,
                             std::vector<int> const &smaller_cut_sizes_to_try
                             );

arma::vec compute_expo_prob_weights_on_edges(
        std::vector<EdgeCut> valid_edges, double alpha, double target);


arma::vec compute_expo_prob_weights_on_smaller_dev_edges(
        std::vector<EdgeCut> valid_edges, double alpha, double target);


arma::vec compute_almost_best_weights_on_smaller_dev_edges(
        std::vector<EdgeCut> valid_edges, double epsilon, double target);


std::vector<EdgeCut> get_all_valid_edges_in_directed_tree(
    const Tree &a_ust, 
    const int root,
    const arma::uvec &pop,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    const int total_region_pop, const int total_region_size,
   const double lower, const double upper, const double target
);

std::vector<EdgeCut> get_all_valid_edges_in_undirected_tree(
    const VertexGraph &a_ust, 
    const int root,
    const arma::uvec &pop,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    const int total_region_pop, const int total_region_size,
   const double lower, const double upper, const double target
);



std::vector<EdgeCut> get_valid_edges_in_joined_tree(
        MapParams const &map_params,
        VertexGraph const &forest_graph, 
        std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
        const int region1_root, const int region1_pop,
        const int region2_root, const int region2_pop,
        const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try,
        const int total_merged_region_size
);




#endif
