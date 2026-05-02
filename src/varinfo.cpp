#include <RcppArmadillo.h>
using namespace Rcpp;

double var_info(
    IntegerVector const &m1, IntegerVector const &m2, 
    NumericVector const &pop, 
    int const k1, int const k2
    ) {
    int V = m1.size();
    NumericMatrix joint(k1, k2);
    NumericVector p1(k1);
    NumericVector p2(k2);

    double total_pop = 0;
    for (int i = 0; i < V; i++) {
        joint(m1[i]-1, m2[i]-1) += pop[i];
        p1[m1[i]-1] += pop[i];
        p2[m2[i]-1] += pop[i];
        total_pop += pop[i];
    }

    double varinf = 0;
    for (int i = 0; i < k1; i++) {
        for (int j = 0; j < k2; j++) {
            double jo = joint(i, j);
            if (jo <= 0) continue;
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
 * 
 */
// [[Rcpp::export]]
NumericVector var_info_vec(IntegerMatrix const &m, IntegerVector const &ref, NumericVector const &pop) {
    int N = m.ncol();

    NumericVector out(N);
    int const ref_k = max(ref); // number of districts in reference
    for (int j = 0; j < N; j++) {
        int const plan_k = max(m(_, j));
        out[j] = var_info(ref, m(_, j), pop, ref_k, plan_k);
    }

    return out;
}
