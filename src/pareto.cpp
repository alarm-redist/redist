#include "smc_base.h"

// [[Rcpp::export]]
LogicalVector pareto_dominated(arma::mat x) {
    int n = x.n_cols;
    int p = x.n_rows;

    LogicalVector dominated(n); // initialized to FALSE
    // for every elelement
    // In backwards order so that new duplicates don't show as undominated
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < n; j++) { // for other non-dominated el
            if (i == j || dominated[j]) continue;

            bool any_less = false;
            bool all_leq = true;
            bool all_eq = true;
            for (int k = 0; k < p; k++) {
                if (x(k, i) < x(k, j)) { // el of j does not dominate i
                    all_leq = false;
                    all_eq = false;
                    break;
                } else if (x(k, i) > x(k, j)) { // el of j strictly less
                    any_less = true;
                    all_eq = false;
                }
            }

            if (all_eq || (all_leq && any_less)) {
            // if (all_leq && any_less) {
                dominated[i] = true;
                break;
            }
        }
    }

    return dominated;
}
