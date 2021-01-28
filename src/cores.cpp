#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List cores(List adj, IntegerVector dm, int k, List cd_within_k) {
  IntegerVector kvec(dm.size());
  IntegerVector temp;
  IntegerVector conncomp(dm.size());
  int d;
  IntegerVector tempcdk, tempcdj;
  
  
  // Step 1: ID border precincts
  for(int i = 0; i < dm.size(); i++){
    d = 0;
    temp = adj(i);
    for(int j = 0; j < temp.size(); j++){
      if(dm(temp(j)) != dm(i)){
        d += 1;
      }
      if(d > 0){
        kvec(i) = 1;
      }
    } 
  }
  
  // Step 2: Identify k = 2,...,K precincts
  if(k > 1){
    for(int curr = 1; curr < k ; curr++){
      for(int i = 0; i < dm.size(); i++){
        if(kvec(i) == 0){
          temp = adj(i);
          for(int j = 0; j < temp.size(); j++){
            if(dm(temp(j)) == dm(i)){
              if(kvec(temp(j)) == curr){
                kvec(i) = curr + 1;
              }
            }
          }
        }
      }

    }
  }
  
  
  // Updated Step 3: Idnetify connected components for k = 0 precincts
  IntegerVector dist_lookup = sort_unique(dm);
  IntegerVector dist_cc(dist_lookup.size());
  IntegerVector currcd(1);
  int cc, s, r, idx;
  IntegerVector reservoir;
  for(int i = 0; i < dm.size(); i++){
    if(conncomp(i) == 0 && kvec(i) == 0){
      // grab next connected component
      currcd = dm(i);
      idx = match(currcd, dist_lookup)(0) - 1;
      dist_cc(idx) ++;  
      cc = dist_cc(idx);
      conncomp(i) = cc;
      
      // Now loop through friends
      temp = adj(i);
      reservoir = IntegerVector(0);
      s = 0;
      for(int j = 0; j < temp.size(); j++){
        if(dm(temp(j)) == dm(i) && kvec(temp(j)) == 0 && conncomp(temp(j)) == 0){
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
            if(dm(temp(j)) == dm(i) && kvec(temp(j)) == 0 && conncomp(temp(j)) == 0){
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
  
  // (Bonus) Step 4: Identify all cds within k:
  List adj_within_k = clone(cd_within_k);
  IntegerVector curradj;
  IntegerVector nextreservoir;
  // Part A: Get all within k adjacent precincts
  for(int i = 0; i < dm.size(); i++){
    curradj = adj(i);
    reservoir = curradj;
    if(k > 1){
      for(int q = 1; q < k; q++){
        nextreservoir = IntegerVector(0);
        for(int  j = 0; j < reservoir.size(); j++){
          temp = adj(reservoir(j));
          nextreservoir = union_(nextreservoir, temp);
        }
        reservoir = union_(reservoir, nextreservoir);
      }
    }
    adj_within_k(i) = clone(reservoir);
  }
  
  // Transform those into congressional districts
  for(int i = 0; i < dm.size(); i++){
    curradj = adj_within_k(i);
    reservoir = IntegerVector(curradj.size());
    for(int j = 0; j < curradj.size(); j++){
      reservoir(j) = dm(curradj(j));
    }
    cd_within_k(i) = unique(reservoir);
  }
  
  
  List ret;
  ret["dm"] = dm;
  ret["k"] = kvec;
  ret["conncomp"] = conncomp;
  ret["cd_within_k"] = cd_within_k;
  
  return ret;
}


// [[Rcpp::export]]
IntegerVector update_conncomp(IntegerVector dm, IntegerVector kvec, List adj){
  IntegerVector conncomp(dm.size());
  IntegerVector temp;
  IntegerVector dist_lookup = sort_unique(dm);
  IntegerVector dist_cc(dist_lookup.size());
  IntegerVector currcd(1);
  int cc, s, r, idx;
  IntegerVector reservoir;
  for(int i = 0; i < dm.size(); i++){
    if(conncomp(i) == 0 && kvec(i) == 0){
      // grab next connected component
      currcd = dm(i);
      idx = match(currcd, dist_lookup)(0) - 1;
      dist_cc(idx) ++;  
      cc = dist_cc(idx);
      conncomp(i) = cc;
      
      // Now loop through friends
      temp = adj(i);
      reservoir = IntegerVector(0);
      s = 0;
      for(int j = 0; j < temp.size(); j++){
        if(dm(temp(j)) == dm(i) && kvec(temp(j)) == 0 && conncomp(temp(j)) == 0){
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
            if(dm(temp(j)) == dm(i) && kvec(temp(j)) == 0 && conncomp(temp(j)) == 0){
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
