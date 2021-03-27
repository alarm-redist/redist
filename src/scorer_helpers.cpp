#include <algorithm>
#include <vector>
#include <RcppArmadillo.h>

#include "smc_base.h"

using namespace Rcpp;


// [[Rcpp::export]]
NumericVector k_smallest(NumericMatrix x, int k = 1) {
    int ncol = x.ncol();
    int nrow = x.nrow();
    std::vector<double> col(nrow);
    NumericVector out(ncol);

    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
            col[j] = x(j, i);
        }
        out[i] = col[select_k(col, k)];
    }

    return out;
}

// [[Rcpp::export]]
NumericVector k_biggest(NumericMatrix x, int k = 1) {
    int ncol = x.ncol();
    int nrow = x.nrow();
    std::vector<double> col(nrow);
    NumericVector out(ncol);

    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nrow; j++) {
            col[j] = x(j, i);
        }
        std::nth_element(col.begin(), col.begin() + k - 1,
                         col.end(), std::less<double>());
        out[i] = col[k - 1];
    }

    return out;
}
