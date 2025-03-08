/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Helper functions for all redist algorithm types
********************************************************/

#include "redist_alg_helpers.h"


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
){
    if(rcpp_mat.ncol() != arma_mat.n_cols || rcpp_mat.nrow() != arma_mat.n_rows){
        throw Rcpp::exception("Arma and Rcpp Matrix are not the same size");
    }

    // go by column because both are column major
    int ncols = (int) rcpp_mat.ncol();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, ncols, [&] (int j) {
        // Copy column i into the rcpp matrix
        std::copy(
            arma_mat.colptr(j), // Start of column in subview
            arma_mat.colptr(j) + arma_mat.n_rows, // End of column in subview
            rcpp_mat.column(j).begin() // Start of column in Rcpp::IntegerMatrix
        );
    });

    // Wait for all the threads to finish
    pool.wait();
    
}



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
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec){

    int M = (int) plan_ptrs_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // reorder every plan
        plan_ptrs_vec.at(i)->reorder_plan_by_oldest_split(*dummy_plan_ptrs_vec.at(i));
    });

    // Wait for all the threads to finish
    pool.wait();

    return;
    
}




std::vector<std::unique_ptr<TreeSplitter>> get_tree_splitters(
    MapParams const &map_params,
    SplittingMethodType const splitting_method,
    Rcpp::List const &control,
    int const nsims
){
    // create the pointer 
    std::vector<std::unique_ptr<TreeSplitter>> tree_splitters_ptr_vec; 
    tree_splitters_ptr_vec.reserve(nsims);

    int V = map_params.V;
    double target = map_params.target;

    if(splitting_method == SplittingMethodType::NaiveTopK){
        // set splitting k to -1
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), nsims, [V] {
            return std::make_unique<NaiveTopKSplitter>(V, -1);
        });
    }else if(splitting_method == SplittingMethodType::UnifValid){
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), nsims, [V] {
            return std::make_unique<UniformValidSplitter>(V);
        });
    }else if(splitting_method == SplittingMethodType::ExpBiggerAbsDev){
        double alpha = as<double>(control["splitting_alpha"]);
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), nsims, [V, alpha, target] {
            return std::make_unique<ExpoWeightedSplitter>(V, alpha, target);
        });
    }else if(splitting_method == SplittingMethodType::ExpSmallerAbsDev){
        double alpha = as<double>(control["splitting_alpha"]);
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), nsims, [V, alpha, target] {
            return std::make_unique<ExpoWeightedSmallerDevSplitter>(V, alpha, target);
        });
    }else if(splitting_method == SplittingMethodType::Experimental){
        double epsilon = as<double>(control["splitting_epsilon"]);
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), nsims, [V, epsilon, target] {
            return std::make_unique<ExperimentalSplitter>(V, epsilon, target);
        });
    }else{
        throw Rcpp::exception("Invalid Splitting Method!");
    }

    return tree_splitters_ptr_vec;
}

