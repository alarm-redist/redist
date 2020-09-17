#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector minkowski(IntegerVector v, IntegerMatrix m, int p) {
  double mink = 0, diff; 
  int i = 0, c = 0;
  NumericVector result(m.ncol());
  for(c = 0; c < m.ncol(); c++){
    mink = 0;
    for(i = 0; i < v.size(); i++){
      diff = pow(abs(v[i] - m(i, c)), p);
      mink+= diff;
    }
    mink = pow(mink, 1.0/p);
    result[c] = mink;
  }
  return result;
}
