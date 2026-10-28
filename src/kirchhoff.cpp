#include "map_calc.h"

#include <sstream>
#include <stdexcept>
#include <vector>

// Both of these used to delegate to redistmetrics, which required arma matrices.
// They now use the templated versions in map_calc.h, which work directly on a
// column of an Rcpp::IntegerMatrix.
//
// A county-and-district piece that is not connected has no spanning trees, so
// log_st_map returns -Inf for that plan, matching what redistmetrics did.

// [[Rcpp::export]]
Rcpp::NumericVector log_st_map(const Graph &g, const Rcpp::IntegerMatrix &districts,
                               const Rcpp::IntegerVector &counties, int n_distr) {
    int const V = static_cast<int>(g.size());
    int const N = districts.ncol();

    if (counties.size() != V) {
        std::ostringstream oss;
        oss << "log_st_map: `counties` must have one entry per vertex.\n";
        oss << "counties.size()=" << counties.size() << "\n";
        oss << "g.size()=" << V << "\n";
        throw std::runtime_error(oss.str());
    }

    if (districts.nrow() != V) {
        std::ostringstream oss;
        oss << "log_st_map: `districts` must have one row per vertex.\n";
        oss << "districts.nrow()=" << districts.nrow() << "\n";
        oss << "g.size()=" << V << "\n";
        throw std::runtime_error(oss.str());
    }

    // Counties are 1-indexed and are used unsigned as `counties[i] - 1` to index
    // county-level bookkeeping, so an NA (which is INT_MIN) or a non-positive
    // value would wrap around to a huge index and read out of bounds.
    for (int i = 0; i < V; i++) {
        if (Rcpp::IntegerVector::is_na(counties[i]) || counties[i] < 1) {
            std::ostringstream oss;
            oss << "log_st_map: `counties` must be positive integers with no NA.\n";
            oss << "counties[" << (i + 1) << "]=";
            if (Rcpp::IntegerVector::is_na(counties[i])) {
                oss << "NA\n";
            } else {
                oss << counties[i] << "\n";
            }
            throw std::runtime_error(oss.str());
        }
    }

    std::vector<unsigned int> counties_vec(counties.begin(), counties.end());

    Rcpp::NumericVector out(N);
    for (int n = 0; n < N; n++) {
        out[n] = eval_log_st(districts.column(n), g, counties_vec, n_distr);
    }

    return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector n_removed(const Graph &g, const Rcpp::IntegerMatrix &districts, int n_distr) {
    int const V = static_cast<int>(g.size());
    int const N = districts.ncol();

    if (districts.nrow() != V) {
        std::ostringstream oss;
        oss << "n_removed: `districts` must have one row per vertex.\n";
        oss << "districts.nrow()=" << districts.nrow() << "\n";
        oss << "g.size()=" << V << "\n";
        throw std::runtime_error(oss.str());
    }

    Rcpp::NumericVector out(N);
    for (int n = 0; n < N; n++) {
        out[n] = eval_er(districts.column(n), g);
    }

    return out;
}
