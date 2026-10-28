
#include "map_calc.h"















/*
 * Create the projective distribution of a variable `x`
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix proj_distr_m(Rcpp::IntegerMatrix districts, const Rcpp::NumericVector x,
                           Rcpp::IntegerVector draw_idx,
                           int n_distr) {
    int n = draw_idx.size();
    int V = districts.nrow();

    Rcpp::NumericMatrix out(V, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < V; j++) {
            int idx = draw_idx[i] - 1;
            out(j, i) = x[n_distr * idx + districts(j, idx) - 1];
        }
    }

    return out;
}



/*
 * Column-wise maximum
 */
// [[Rcpp::export]]
Rcpp::NumericVector colmax(const Rcpp::NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    Rcpp::NumericVector out(ncol);
    for (int j = 0; j < ncol; j++) {
        double best = x(0, j);
        for (int i = 1; i < nrow; i++) {
            if (x(i, j) > best) {
                best = x(i, j);
            }
        }
        out[j] = best;
    }

    return out;
}
/*
 * Column-wise minimum
 */
// [[Rcpp::export]]
Rcpp::NumericVector colmin(const Rcpp::NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    Rcpp::NumericVector out(ncol);
    for (int j = 0; j < ncol; j++) {
        double best = x(0, j);
        for (int i = 1; i < nrow; i++) {
            if (x(i, j) < best) {
                best = x(i, j);
            }
        }
        out[j] = best;
    }

    return out;
}









double compute_log_pop_temper(double const target, double const pop_temper, int const ndists,
                              int const region_pop, int const region_size) {
    double region_target = target * region_size;
    // get population deviation
    double const pop_dev =
        std::fabs(static_cast<double>(region_pop) - region_target) / region_target;

    double const pop_pen =
        std::sqrt(static_cast<double>(ndists) - 2) * std::log(1e-12 + pop_dev);

    // now return the values for the old region minus the two new ones
    return pop_pen * pop_temper;
}