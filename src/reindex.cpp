#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix reindex(IntegerMatrix dm, int nd) {
  IntegerMatrix ordered(dm.nrow(), dm.ncol());
  IntegerVector dict = seq_len(nd);
  int i = 0;
  IntegerVector iv(1);
  
  IntegerVector colc(dm.nrow());
  for(int c = 0; c < ordered.ncol(); c++){
    dict = rep(0, nd);
    colc = dm(_,c);
    dict(0) = colc(0);
    i = 1;
    for(int r = 1; r < ordered.nrow(); r++){
      iv = colc(r);
      if(!(in(iv, dict))(0)){
        dict(i) = colc(r);
        i++;
      }
    }

    colc =  match(colc, dict);
    ordered(_,c) = colc;
  }
  return ordered;
}

