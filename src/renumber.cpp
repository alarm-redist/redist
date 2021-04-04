#include <Rcpp.h>
using namespace Rcpp;

// row of output is m1 districts
// col of output is m2 districts
// [[Rcpp::export]]
NumericMatrix plan_joint(IntegerVector m1, IntegerVector m2, NumericVector pop) {
    int k = max(m1);
    int V = m1.size();
    NumericMatrix joint(k);
    NumericVector p1(k);
    NumericVector p2(k);

    for (int i = 0; i < V; i++) {
        joint(m1[i]-1, m2[i]-1) += pop[i];
        p1[m1[i]-1] += pop[i];
        p2[m2[i]-1] += pop[i];
    }

    return joint;
}

// [[Rcpp::export]]
IntegerMatrix renumber_matrix(IntegerMatrix plans, IntegerVector renumb) {
    int V = plans.nrow();
    int N = plans.ncol();
    int n_distr = max(plans(_, 1));
    IntegerMatrix out(V, N);

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < V; i++) {
            out(i, j) = renumb(n_distr*j + plans(i, j) - 1);
        }
    }

    return out;
}

// [[Rcpp::export]]
int best_renumber(IntegerMatrix combn, NumericMatrix joint) {
    int n = combn.nrow();
    int n_distr = combn.ncol();
    double best_val = -1;
    int best_idx = -1;

    for (int i = 0; i < n; i++) {
        double diag_pop = 0;
        for (int j = 0; j < n_distr; j++) {
            diag_pop += joint(j, combn(i, j) - 1);
        }

        if (diag_pop > best_val) {
            best_val = diag_pop;
            best_idx = i;
        }
    }

    return best_idx + 1;
}
