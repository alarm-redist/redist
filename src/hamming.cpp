#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int hamming(IntegerVector x, IntegerVector y) {
   int ham = 0;
  for(int i = 0; i < x.size(); i++){
    if(x[i] != y[i]){
      ham++;
    }
  }
  return ham;
}