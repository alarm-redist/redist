#include "smc_base.h"
#include <redistmetrics.h>

// [[Rcpp::export]]
IntegerVector splits(IntegerMatrix dm, IntegerVector community, int nd, int max_split) {
    return redistmetrics::splits(dm, community, nd, max_split);
}
// [[Rcpp::export]]
IntegerMatrix dist_cty_splits(IntegerMatrix dm, IntegerVector community, int nd) {
  IntegerMatrix ret(nd, dm.ncol());
  IntegerVector com_name = sort_unique(community);
  IntegerVector com_found(com_name.size(), 0);

  // by column (aka map)
  for(int c = 0; c < dm.ncol(); c++){
    // by district
   for(int d = 0; d < nd; d++){
     com_found = IntegerVector(com_found.size(), 0);
     // across all rows
     for(int r = 0; r < dm.nrow(); r++){
       if (dm(r,c) == d) {
         com_found(community(r)) = 1;
       }
     }
     ret(d, c) = sum(com_found);
   }
  }
  return ret;
}
