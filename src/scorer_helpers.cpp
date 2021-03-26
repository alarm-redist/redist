#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector k_smallest(NumericMatrix x, int k = 1) {
  NumericMatrix sorted(x.nrow(), x.ncol());
  NumericVector temp(x.nrow());
  
  
  for(int i = 0; i < x.ncol(); i++){
    temp = x(_, i);
    sorted(_, i) = temp.sort();
  }
  
  NumericVector out = sorted(k - 1, _);
  
  return out;
}

// [[Rcpp::export]]
NumericVector k_biggest(NumericMatrix x, int k = 1) {
  NumericMatrix sorted(x.nrow(), x.ncol());
  NumericVector temp(x.nrow());
  
  
  for(int i = 0; i < x.ncol(); i++){
    temp = x(_, i);
    sorted(_, i) = temp.sort();
  }
  
  NumericVector out = sorted(x.nrow() - k, _);
  
  return out;
}
