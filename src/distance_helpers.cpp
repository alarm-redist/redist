#include <Rcpp.h>
using namespace Rcpp;


double distance(double x1, double x2,
                double y1, double y2){
  double dist = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
  return(dist);
}


NumericMatrix distance_matrix(NumericVector x, NumericVector y) {
  // init matrix
  int dim = x.length();
  NumericMatrix out(dim, dim);
  // loop over rows/columns for symmetric
  for(int r = 0; r < dim; r++){
    for(int c = r+1; c < dim; c++){
      double d = distance(x[r], x[c], y[r], y[c]);
      out(r, c) = d;
      out(c, r) = d;
    }
  }

  return(out);
}


int closest_adj(IntegerVector adj,
                int i_dist,
                NumericVector x,
                NumericVector y){

  if(adj.size() == 1){
    return adj[0];
  }

  // inits
  double min_distance = distance(x[i_dist], x[adj[0]], y[i_dist], y[adj[0]]);
  int min_district = adj[0];
  int curr_distance;


  // loop
  for(int j = 1; j < adj.size(); j++){
    curr_distance = distance(x[i_dist], x[adj[j]], y[i_dist], y[adj[j]]);
    if(min_distance > curr_distance){
      min_distance = curr_distance;
      min_district = adj[j];
    }
  }

  return min_district;
}

// [[Rcpp::export]]
double dist_dist_diff(int p,
                      int i_dist,
                      int j_dist,
                      NumericVector x_center,
                      NumericVector y_center,
                      NumericVector x,
                      NumericVector y){
  double dist_i_dist, dist_j_dist;
  dist_i_dist = distance(x_center[p], x[i_dist], y_center[p], y[i_dist]);
  dist_j_dist = distance(x_center[p], x[j_dist], y_center[p], y[j_dist]);
  return dist_i_dist - dist_j_dist;
}
