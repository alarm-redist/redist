#pragma once
#ifndef MANUAL_SPLITTING_H
#define MANUAL_SPLITTING_H

// [[Rcpp::depends(redistmetrics)]]

#include "smc_base.h"

#include <string>
#include <cmath>
#include <functional>
#include <cli/progress.h>
#include <RcppThread.h>

#include <kirchhoff_inline.h>
#include "wilson.h"
#include "tree_op.h"
#include "map_calc.h"
#include "splitting.h"
#include "merging.h"





// ' Draws a spanning tree uniformly at random on a region and returns it
// '
// ' Draws a spanning tree uniformly at random on a region of a plan using
// ' Wilson's algorithm. 
// '
// ' @title Draw a uniformly random spanning tree on a region of a plan
// '
// '
// ' @param adj_list A 0-indexed adjacency list representing the undirected graph
// ' which represents the underlying map the plans are to be drawn on
// ' @param counties Vector of county labels of each vertex in `g`
// ' @param pop A vector of the population associated with each vertex in `g`
// ' @param ndists The number of districts the final plans will have
// ' @param num_regions The number of regions in the inputted plan
// ' @param num_districts The number of districts in the inputted plan
// ' @param region_id_to_draw_tree_on The id of the region in the plan to draw
// ' the tree on. 
// ' @param lower Acceptable lower bounds on a valid district's population
// ' @param upper Acceptable upper bounds on a valid district's population
// ' @param region_ids A V by 1 matrix with the region ids of each vertex
// ' @param region_dvals A N by 1 matrix with the sizes of each regions 
// ' @param verbose Whether or not to print out the inputted plan before
// ' attemping to draw a tree. 
//'
//' @returns A list with the following 
//'     - `uncut_tree`: The spanning tree drawn on the region stored as a
//'     0-indexed directed edge adjacency graph.
//'     - `num_attempts`: The number of attempts it took to draw the tree.
//' 
//' @keywords internal
// [[Rcpp::export]]
List draw_a_tree_on_a_region(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int ndists, int num_regions, int num_districts,
    int region_id_to_draw_tree_on,
    double lower, double upper,
    arma::umat region_ids, arma::umat region_dvals,
    bool verbose
);


// [[Rcpp::export]]
List perform_a_valid_region_split(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int N, int num_regions, int num_districts,
    int region_id_to_split,
    double target, double lower, double upper,
    arma::umat region_ids, arma::umat region_dvals,
    int split_dval_min, int split_dval_max, 
    bool verbose = false, int k_param = 1
);

// [[Rcpp::export]]
List perform_merge_split_steps(
        List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int k_param,
        double target, double lower, double upper,
        int N, int num_regions, int num_districts,
        arma::umat region_ids, arma::umat region_dvals,
        std::vector<int> region_pops,
        bool split_district_only, int num_merge_split_steps,
        bool verbose
);

#endif
