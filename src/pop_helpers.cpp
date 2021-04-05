#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int closest_adj_pop(IntegerVector adj,
                int i_dist,
                NumericVector g_prop){
  
  if(adj.size() == 1){
    return adj[0];
  }
  
  // inits 
  double min_distance = fabs( g_prop(i_dist) - g_prop(adj(0)) ); 
  int min_district = adj[0];
  int curr_distance;
  
  
  // loop
  for(int j = 1; j < adj.size(); j++){
    curr_distance = fabs( g_prop(i_dist) - g_prop(adj(j)) );
    if(min_distance > curr_distance){
      min_distance = curr_distance;
      min_district = adj[j];
    }
  }
  
  return min_district;
}

