// [[Rcpp::depends(redistmetrics)]]
#include "map_calc.h"
#include <redistmetrics.h>



/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const arma::subview_col<arma::uword> &districts, int distr, const arma::uvec &total_pop,
                     arma::mat ssdmat, double denominator = 1.0) {
    arma::uvec idxs = find(districts == distr);
    double ssd = 0.0;

    for (int i = 0; i < idxs.size() - 1; i++) {
        for (int k = i + 1; k < idxs.size(); k++) {
            ssd += (double)ssdmat(idxs(i), idxs(k)) * total_pop(idxs(i)) * total_pop(idxs(k));
        }
    }

    return ssd / denominator;
}

/*
 * Compute the qps penalty for district `distr`
 */
double eval_qps(const arma::subview_col<arma::uword> &districts, int distr, const arma::uvec &total_pop,
                const arma::uvec &cities, int n_city, int nd) {

    arma::vec tally(n_city);
    arma::vec pj(n_city);
    arma::vec j(n_city);
    arma::vec sumpj(n_city);

    arma::uvec idxs_d = find(districts == distr);
    double pop = arma::sum(total_pop(idxs_d));

    for (int i = 0; i < n_city; i++) {
        arma::uvec idxs = find(cities == (i + 1));
        idxs = arma::intersect(idxs_d, idxs);
        tally(i) = sum(total_pop(idxs));
        if (tally(i) > 0) {
            j(i) += 1;
        }
    }

    pj = tally / pop;
    sumpj = pj * (1.0 - pj);
    sumpj = sumpj / (double)nd;

    return sum(sumpj) + log(sum(j));
}

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_log_st(const arma::subview_col<arma::uword> &districts, const Graph g, arma::uvec counties,
                   int ndists) {
    return (double)redistmetrics::log_st_map(g, districts, counties, ndists)[0];
}

/*
 * Compute the edges removed penalty for district `distr`
 */
double eval_er(const arma::subview_col<arma::uword> &districts, const Graph g, int ndists) {
    return (double)redistmetrics::n_removed(g, districts, ndists)[0];
}










/*
 * Create the projective distribution of a variable `x`
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix proj_distr_m(Rcpp::IntegerMatrix districts, const arma::vec x, Rcpp::IntegerVector draw_idx,
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