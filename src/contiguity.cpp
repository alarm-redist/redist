#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector contiguity(List adj, IntegerVector group) {
  IntegerVector choices = sort_unique(group);
  IntegerVector conncomp(group.size());
  IntegerVector group_cc(choices.size());
  IntegerVector currgroup(1);
  IntegerVector temp;
  int cc, s, r, idx;
  IntegerVector reservoir;
  for(int i = 0; i < group.size(); i++){
    if(conncomp(i) == 0){
      currgroup = group(i);
      idx = match(currgroup, choices)(0) - 1;
      group_cc(idx) ++;  
      cc = group_cc(idx);
      conncomp(i) = cc;
      
      temp = adj(i);
      reservoir = IntegerVector(0);
      s = 0;
      for(int j = 0; j < temp.size(); j++){
        if(group(temp(j)) == group(i) && conncomp(temp(j)) == 0){
          reservoir.push_back(temp(j));
          conncomp(temp(j)) = cc;
          s++;
        } 
      }
      
      if(s > 0){
        r = 0;
        while(r < s){
          temp = adj(reservoir(r));
          for(int j = 0; j < temp.size(); j++){
            if(group(temp(j)) == group(i) && conncomp(temp(j)) == 0){
              reservoir.push_back(temp(j));
              conncomp(temp(j)) = cc;
              s++;
            } 
          }
          r++;
        }
      }
    }
  }
  
  return conncomp;
}

