#pragma once
#ifndef GRAPH_OP_H
#define GRAPH_OP_H

#include <vector>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "redist_types.h"


// [[Rcpp::export]]
RegionMultigraphCount get_region_multigraph(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

// [[Rcpp::export]]
arma::mat get_region_laplacian(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

double compute_log_region_multigraph_spanning_tree(
    RegionMultigraphCount const &region_multigraph
);

// [[Rcpp::export]]
double get_log_number_linking_edges(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

// [[Rcpp::export]]
double get_merged_log_number_linking_edges(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids,
    int const region1_id, int const region2_id
);

RegionMultigraphCount build_region_multigraph(
    Graph const &g, 
    PlanVector const &region_ids,
    int const num_regions
);



double get_log_merged_region_multigraph_spanning_tree(
    RegionMultigraphCount const &region_multigraph,
    std::vector<int> &merge_index_reshuffle,
    int region1_id, int region2_id
);

#endif
