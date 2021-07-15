#include "smc_base.h"

// [[Rcpp::export]]
IntegerVector splits(IntegerMatrix dm, IntegerVector community) {
  IntegerVector ret(dm.ncol());
  IntegerVector com_name = sort_unique(community);
  IntegerVector com_bin(com_name.size());
  List coms(com_name.size());
  IntegerVector idx = seq_len(community.size())-1;
  IntegerVector temp;
  IntegerVector colc(dm.nrow());

// create indices for each community
  for(int i = 0; i < com_name.size(); i++){
    temp = idx[community == com_name(i)];
    coms(i) = temp;
  }
  IntegerVector firstcom(1);
  for(int c =0; c < dm.ncol(); c++){
    colc = dm(_,c);
    for(int i = 0; i < com_name.size(); i++){
      idx = coms(i);
      temp = colc[idx];
      if(is_true(all(temp(0) == temp))){
        com_bin(i) = 0;
      } else{
        com_bin[i] = 1;
      }
    }
    ret[c] = sum(com_bin);
  }
  return ret;
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
