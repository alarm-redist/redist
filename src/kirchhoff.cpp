// [[Rcpp::depends(redistmetrics)]]
#include <RcppArmadillo.h>
#include <redistmetrics.h>

// NOTE: I think we no longer need this dependency anymore? There are faster internal versions of this function

// [[Rcpp::export]]
Rcpp::NumericVector log_st_map(const Graph &g, const arma::umat &districts,
                         const arma::uvec &counties, int n_distr) {
    return redistmetrics::log_st_map(g, districts, counties, n_distr);
}

// [[Rcpp::export]]
Rcpp::NumericVector n_removed(const Graph &g, const arma::umat &districts, int n_distr) {
    return redistmetrics::n_removed(g, districts, n_distr);
}
