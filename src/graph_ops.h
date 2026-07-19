#pragma once
#ifndef GRAPH_OP_H
#define GRAPH_OP_H

#include "advanced_types.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>

// [[Rcpp::export]]
RegionMultigraphCount get_region_multigraph(Rcpp::List const &adj_list,
                                            arma::uvec const &region_ids);

// [[Rcpp::export]]
arma::mat get_region_laplacian(Rcpp::List const &adj_list, arma::uvec const &region_ids);

RegionMultigraphCount build_region_multigraph(Graph const &g, PlanVector const &region_ids,
                                              int const num_regions);

#endif
