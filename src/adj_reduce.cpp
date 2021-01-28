#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List reduce_adj(List adj_list, IntegerVector prec_keep, IntegerVector prec_idx) {
  int n_keep = sum(prec_keep);
  List adj(n_keep);
  IntegerVector full, sub;
  
  
  int fill = 0;
  for(int i = 0; i < adj_list.size(); i++){
    if(prec_keep[i] == 1){
      full = adj_list[i];
      sub = intersect(full, prec_idx);
       sub = match(sub, prec_idx)-1;
       adj[fill] = sub.sort();
    }
    fill += prec_keep[i];
  }

  return adj;
}
