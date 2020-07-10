#include "kirchhoff.h"

/*
 * Compute the log number of spanning trees which could generate a given set of maps.
 * `districts` should have each column be a map
 */
// TESTED
NumericVector log_st_map(const List &g, const IntegerMatrix &districts,
                         const IntegerVector &counties, int n_distr) {
    int N = districts.ncol();
    int n_cty = max(counties);
    NumericVector log_st(N);
    for (int i = 0; i < N; i++) {
        double accuml = 0;
        for (int d = 1; d <= n_distr; d++) { // districts are 1-indexed
            for (int j = 1; j <= n_cty; j++) {
                accuml += log_st_distr(g, districts, counties, i, d, j);
            }
            accuml += log_st_contr(g, districts, counties, n_cty, i, d);
        }
        log_st[i] = accuml;
    }
    return log_st;
}

/*
 * Compute the log number of spanning trees for `district` intersect `county`
 */
// TESTED
double log_st_distr(const List &g, const IntegerMatrix &districts,
                    const IntegerVector &counties, int idx,
                    int district, int county) {
    int V = g.size();
    // number of precincts in this district
    int K = 0;
    std::vector<int> pos(V); // keep track of positions in subgraph
    int start = 0; // where to start loop below, to save time
    for (int i = 0; i < V; i++) {
        pos[i] = K - 1; // minus one because we're dropping 1st row and column
        if (districts(i, idx) == district && counties[i] == county) {
            K++;
            if (K == 2) start = i; // start 2nd vertex
        }
    }
    if (K <= 1) return 0;

    mat adj = zeros(K-1, K-1); // adjacency matrix (minus 1st row and column)
    for (int i = start; i < V; i++) {
        if (districts(i, idx) != district || counties[i] != county) continue;

        int prec = pos[i];
        if (prec < 0) continue;
        std::vector<int> nbors = as<std::vector<int>>(g[i]);
        int length = nbors.size();
        int degree = 0; // keep track of index within subgraph
        for (int j = 0; j < length; j++) {
            int nbor = nbors[j];
            if (districts(nbor, idx) != district || counties[nbor] != county) continue;
            degree++;
            if (pos[nbor] < 0) continue;
            adj(prec, pos[nbor]) = -1;
        }
        adj(prec, prec) = degree;
    }

    double lst, sign;
    log_det(lst, sign, adj);
    return lst;
}

/*
 * Compute the log number of spanning trees for the contracted graph
 */
// TESTED
double log_st_contr(const List &g, const IntegerMatrix &districts,
                    const IntegerVector &counties, int n_cty, int idx,
                    int district) {
    if (n_cty == 1) return 0;
    int V = g.size();
    // number of counties in this district
    int K = 0;
    std::vector<int> pos(V); // keep track of positions in subgraph
    std::vector<int> seen(n_cty, -2); // county lookup
    int start = 0;
    for (int i = 0; i < V; i++) {
        if (districts(i, idx) != district) continue;

        if (seen[counties[i]-1] < 0) {
            pos[i] = K - 1; // minus one because we're dropping 1st row and column
            seen[counties[i]-1] = K;
            K++;
            if (K == 2) start = i; // start 2nd vertex
        } else {
            pos[i] = seen[counties[i]-1] - 1;
        }
    }
    if (K <= 1) return 0;

    mat adj = zeros(K-1, K-1);
    for (int i = start; i < V; i++) {
        if (districts(i, idx) != district) continue;

        int cty = pos[i];
        if (cty < 0) continue; // skip 1st row, col
        std::vector<int> nbors = as<std::vector<int>>(g[i]);
        int length = nbors.size();
        for (int j = 0; j < length; j++) {
            int nbor = nbors[j];
            if (districts(nbor, idx) != district || pos[nbor] == cty) continue;
            adj(cty, cty)++;
            if (pos[nbor] < 0) continue; // ignore 1st row and column
            adj(cty, pos[nbor])--;
        }
    }

    double lst, sign;
    log_det(lst, sign, adj);
    return lst;
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
            std::vector<int> nbors = as<std::vector<int>>(g[i]);
            int length = nbors.size();
            for (int j = 0; j < length; j++) {
                if (districts(nbors[j], n) != dist) removed += 1.0;
            }
        }
        n_rem[n] = removed;
    }

    return n_rem;
}