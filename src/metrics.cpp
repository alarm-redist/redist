#include "smc_base.h"

// [[Rcpp::export]]
IntegerMatrix agg_p2d(IntegerMatrix dm, IntegerVector vote, int nd) {
  IntegerMatrix mat = IntegerMatrix(nd, dm.ncol());
  for(int j = 0; j < dm.ncol(); j++){
    for(int i = 0; i < dm.nrow(); i++){
      mat(dm(i,j)-1,j) += vote[i];
    }
  }
  return mat;
}
