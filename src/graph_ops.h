#ifndef GRAPH_OP_H
#define GRAPH_OP_H

#include <vector>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "gredist_types.h"


// [[Rcpp::export]]
RegionMultigraph get_region_multigraph(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

// [[Rcpp::export]]
arma::imat get_region_laplacian(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

// [[Rcpp::export]]
double get_log_number_linking_edges(
    Rcpp::List const &adj_list,
    arma::uvec const &region_ids
);

RegionMultigraph get_region_multigraph(
    Graph const &g, 
    arma::subview_col<arma::uword> const &region_ids,
    int const num_regions
);

#endif