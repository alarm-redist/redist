#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

/*
 * Compute the log number of spanning trees for `district`
 */
// TESTED
double log_st_distr(const List &g, const IntegerMatrix &districts, int idx, int district) {
  int V = g.size();
  // number of precincts in this district
  int K = 0;
  std::vector<int> pos(V); // keep track of positions in subgraph
  int start = 0; // where to start loop below, to save time
  for (int i = 0; i < V; i++) {
    pos[i] = K - 1; // minus one because we're dropping 1st row and column
    K += districts(i, idx) == district;
    if (K == 2 && districts(i, idx) == district) start = i; // start 2nd vertex
  }
  if (K <= 1) return 0; 
  mat adj = zeros(K-1, K-1); // adjacency matrix (minus 1st row and column)
  for (int i = start; i < V; i++) { 
    if (districts(i, idx) != district) continue;
    std::vector<int> nbors = as<std::vector<int> >(g[i]);
    int length = nbors.size();
    int degree = 0; // keep track of index within subgraph
    for (int j = 0; j < length; j++) {
      int nbor = nbors[j];
      if (districts(nbor, idx) != district) continue; 
      degree++;
      if (pos[nbor] < 0) continue; 
      adj(pos[i], pos[nbor]) = -1;
    }
    adj(pos[i], pos[i]) = degree;
  }
  double lst, sign;
  log_det(lst, sign, adj);
  return lst;
}



/*
 * Compute the log number of spanning trees which could generate a given set of maps.
 * `districts` should have each column be a map
 */
// TESTED
// [[Rcpp::export]]
NumericVector log_st_map(const List &g, const IntegerMatrix &districts, int n_distr) {
  int N = districts.ncol();
  NumericVector log_st(N);
  for (int i = 0; i < N; i++) {
    double accuml = 0;
    for (int d = 1; d <= n_distr; d++) { // districts are 1-indexed
      accuml += log_st_distr(g, districts, i, d);
    }
    log_st[i] = accuml;
  }
  return log_st;
}

/*
 * Compute the number of edges removed
 */
// TESTED
// [[Rcpp::export]]
NumericVector n_removed(const List &g, const IntegerMatrix &districts, int n_distr) {
  int V = g.size();
  int N = districts.ncol();
  NumericVector n_rem(N);
  for (int n = 0; n < N; n++) {
    double removed = 0.0;
    for (int i = 0; i < V; i++) { 
      int dist = districts(i, n);
      std::vector<int> nbors = as<std::vector<int> >(g[i]);
      int length = nbors.size();
      for (int j = 0; j < length; j++) {
        if (districts(nbors[j], n) != dist) removed += 1.0; 
      }
    }
    n_rem[n] = removed;
  }
  return n_rem;
}