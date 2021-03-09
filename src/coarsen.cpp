#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
List coarsen_adjacency(List adj, IntegerVector groups) {
  List adjclean = clone(adj);
  IntegerVector temp(0);
  IntegerVector temp2;
  
  
  List fill = no_init(max(groups) + 1);
  for(int f = 0; f < fill.size(); f++){
    fill(f) = clone(temp);
  }
  
  for(int i = 0;  i < adj.size(); i++){
    temp = adjclean(i);
    for(int j = 0; j < temp.size(); j++){
      temp(j) = groups(temp(j));
    }
    adjclean(i) = temp;
  }
  
  
  for(int f = 0; f < fill.size(); f++){
    temp = fill(f);
      for(int i = 0;  i < adj.size(); i++){
        if(groups(i) == f){
          temp2 = adjclean(i);
          temp = union_(temp, temp2);
        }
      }
    fill(f) = temp;
  }
  
  for(int f = 0; f < fill.size(); f++){
    temp = fill(f);
    for(int j = 0; j < temp.size(); j++){
      if(temp(j) == f){
        temp.erase(j);
        break;
      }
    }
    fill(f) = temp;
  }
  
  return fill;
}

