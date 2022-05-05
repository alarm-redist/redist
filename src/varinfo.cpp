#include <RcppArmadillo.h>
using namespace Rcpp;

double var_info(IntegerVector m1, IntegerVector m2, NumericVector pop, int k) {
    int V = m1.size();
    NumericMatrix joint(k);
    NumericVector p1(k);
    NumericVector p2(k);

    double total_pop = 0;
    for (int i = 0; i < V; i++) {
        joint(m1[i]-1, m2[i]-1) += pop[i];
        p1[m1[i]-1] += pop[i];
        p2[m2[i]-1] += pop[i];
        total_pop += pop[i];
    }

    double varinf = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
            double jo = joint(i, j);
            if (jo < 1) continue;
            varinf -= (jo / total_pop) * (2.0*std::log(jo) - std::log(p1[i]) - std::log(p2[j]));
        }
    }

    if (std::fabs(varinf) <= 1e-9)
        varinf = 0;
    return varinf;
}


/*
 * `m` has rows = precincts, cols = plans
 * `ref` is the plan we want to compute distances to
 * `pop` is population of precincts
 */
// [[Rcpp::export]]
NumericVector var_info_vec(IntegerMatrix m, IntegerVector ref, NumericVector pop) {
    int N = m.ncol();

    NumericVector out(N);
    int k = max(ref);
    for (int j = 0; j < N; j++) {
        out[j] = var_info(ref, m(_, j), pop, k);
    }

    return out;
}
