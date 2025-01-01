#pragma once
#ifndef SMC_ALG_HELPERS_H
#define SMC_ALG_HELPERS_H

// [[Rcpp::depends(redistmetrics)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <limits>
#include <RcppThread.h>

#include "gredist_types.h"



//' Copies data from an arma Matrix into an Rcpp Matrix
//'
//' Takes an arma matrix subview and copies all the data into an RcppMatrix
//' of the same size using the Rcpp Threadpool to copy in parallel. 
//'
//'
//' @title Copies data from an arma Matrix into an Rcpp Matrix
//'
//' @param pool A threadpool for multithreading
//' @param arma_mat Subview of an arma unsigned integer matrix 
//' @param rcpp_mat A matrix of integers with the same size as the arma_mat
//'
//' @details Modifications
//'    - The `rcpp_mat` is filled in with the data om the arma matrix subview
//'
//' @noRd
//' @keywords internal
void copy_arma_to_rcpp_mat(
    RcppThread::ThreadPool &pool,
    arma::subview<arma::uword> arma_mat,
    Rcpp::IntegerMatrix &rcpp_mat
);



//' Reorders all the plans in the vector by order a region was split
//'
//' Takes a vector of plans and uses the vector of dummy plans to reorder
//' each of the plans by the order a region was split.
//'
//'
//' @title Reorders all the plans in the vector by order a region was split
//'
//' @param pool A threadpool for multithreading
//' @param plans_vec A vector of plans
//' @param dummy_plans_vec A vector of dummy plans 
//'
//' @details Modifications
//'    - Each plan in the `plans_vec` object is reordered by when the region was split
//'    - Each plan is a shallow copy of the plans in `plans_vec`
//'
//' @noRd
//' @keywords internal
void reorder_all_plans(
    RcppThread::ThreadPool &pool,
    std::vector<Plan> &plans_vec, 
    std::vector<Plan> &dummy_plans_vec);


#endif
