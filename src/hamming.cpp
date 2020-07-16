#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector hamming(IntegerVector v, IntegerMatrix m) {
  int ham = 0, i = 0, c = 0;
  IntegerVector result(m.ncol());
  for(c = 0; c < m.ncol(); c++){
    ham = 0;
    for(i = 0; i < v.size(); i++){
      if(v[i] != m(i, c)){
        ham++;
      }
    }
    result[c] = ham;
  }
  return result;
}