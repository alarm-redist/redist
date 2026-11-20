/********************************************************
 * Purpose: Rank-normalized, folded R-hat convergence diagnostic.
 *
 * This mirrors the R implementation in R/diagnostics.R. Sums are accumulated
 * in long double, so results can differ from R in the last ulp; that is well
 * below anything the diagnostic is sensitive to. The work is parallelized
 * over (statistic, district) pairs, which are independent of each other.
 ********************************************************/

#include <Rcpp.h>
#include <Rmath.h>
#include <RcppThread.h>

#include "threading_helpers.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

/*
 * Mean with long double accumulation and a correction pass, matching the
 * algorithm R's `mean()` uses. The second pass cancels the rounding error
 * left by the first.
 */
double r_mean(std::vector<double> const &x, std::vector<int> const &idx) {
    std::size_t const n = idx.size();

    long double s = 0.0L;
    for (int const i : idx) s += x[i];
    s /= n;

    if (R_FINITE(static_cast<double>(s))) {
        long double t = 0.0L;
        for (int const i : idx) t += (x[i] - s);
        s += t / n;
    }

    return static_cast<double>(s);
}

/*
 * R's `var()` for a single vector: sum of squared deviations over n - 1.
 */
double r_var(std::vector<double> const &x, std::vector<int> const &idx) {
    std::size_t const n = idx.size();
    if (n < 2) return NA_REAL;

    double const m = r_mean(x, idx);

    long double s = 0.0L;
    for (int const i : idx) {
        long double const d = x[i] - m;
        s += d * d;
    }

    return static_cast<double>(s / (n - 1));
}

double r_mean_all(std::vector<double> const &v) {
    std::size_t const n = v.size();

    long double s = 0.0L;
    for (double const d : v) s += d;
    s /= n;

    if (R_FINITE(static_cast<double>(s))) {
        long double t = 0.0L;
        for (double const d : v) t += (d - s);
        s += t / n;
    }

    return static_cast<double>(s);
}

double r_var_all(std::vector<double> const &v) {
    std::size_t const n = v.size();
    if (n < 2) return NA_REAL;

    double const m = r_mean_all(v);

    long double s = 0.0L;
    for (double const value : v) {
        long double const e = value - m;
        s += e * e;
    }

    return static_cast<double>(s / (n - 1));
}

/*
 * R's `median()`: for even n the two middle order statistics are averaged.
 */
double r_median(std::vector<double> scratch) {
    std::size_t const n = scratch.size();
    if (n == 0) return NA_REAL;

    std::size_t const half = n / 2;
    std::nth_element(scratch.begin(), scratch.begin() + half, scratch.end());
    double const hi = scratch[half];

    if (n % 2 == 1) return hi;

    // Even n averages the two middle order statistics.
    double const lo = *std::max_element(scratch.begin(), scratch.begin() + half);
    return static_cast<double>(
        (static_cast<long double>(lo) + static_cast<long double>(hi)) / 2.0L
    );
}

/*
 * Rank normalization: average-ties ranks mapped through qnorm(r / (n + 1)).
 *
 * Ties must take the mean of the ranks they span, matching
 * `rank(x, ties.method = "average")`. Anything else makes the transform
 * depend on row order rather than on the values, which would let storage
 * order change the diagnostic.
 *
 * Because every member of a tied run maps to the same normal score, qnorm is
 * evaluated once per distinct value rather than once per element. On the
 * discrete statistics, where ties dominate, that is most of the speedup.
 */
void ranknorm(std::vector<double> const &x, std::vector<double> &out,
              std::vector<std::pair<double, int>> &buf) {
    int const n = static_cast<int>(x.size());

    buf.resize(n);
    for (int i = 0; i < n; ++i) buf[i] = std::make_pair(x[i], i);

    std::sort(buf.begin(), buf.end(),
              [](std::pair<double, int> const &a,
                 std::pair<double, int> const &b) { return a.first < b.first; });

    out.resize(n);
    double const denom = static_cast<double>(n) + 1.0;

    int i = 0;
    while (i < n) {
        int j = i;
        double const value = buf[i].first;
        while (j + 1 < n && buf[j + 1].first == value) ++j;

        // ranks i+1 through j+1, averaged
        double const avg_rank = (static_cast<double>(i) + j + 2.0) / 2.0;
        double const z = ::Rf_qnorm5(avg_rank / denom, 0.0, 1.0, 1, 0);

        for (int k = i; k <= j; ++k) out[buf[k].second] = z;

        i = j + 1;
    }
}

double calc_rhat(std::vector<double> const &z,
                 std::vector<std::vector<int>> const &groups,
                 double const n_g,
                 std::vector<double> &means,
                 std::vector<double> &vars) {
    std::size_t const k = groups.size();

    means.resize(k);
    vars.resize(k);

    for (std::size_t i = 0; i < k; ++i) {
        means[i] = r_mean(z, groups[i]);
        vars[i] = r_var(z, groups[i]);
    }

    double const var_between = n_g * r_var_all(means);
    double const var_within = r_mean_all(vars);

    return std::sqrt((var_between / var_within + n_g - 1.0) / n_g);
}

} // namespace

// Compute rank-normalized folded R-hats for every (statistic, district) pair
//
// The grouping is built in R and handed over here so that the split-chain
// convention lives in exactly one place.
//
// @param values A list of numeric columns, one per statistic, each of length
//   equal to the number of rows.
// @param row_index 0-indexed rows, grouped by district and laid out
//   contiguously (a CSR-style layout).
// @param district_start Offsets into `row_index`, of length n_districts + 1.
// @param group_index For each entry of `row_index`, the 0-indexed group
//   (chain, or chain half when splitting) that row belongs to.
// @param n_groups Total number of distinct groups.
// @param num_threads Threads to use; 0 means all available.
//
// @returns An n_districts by n_statistics matrix of R-hat values.
//
// @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix compute_rhats_cpp(Rcpp::List const &values,
                                      Rcpp::IntegerVector const &row_index,
                                      Rcpp::IntegerVector const &district_start,
                                      Rcpp::IntegerVector const &group_index,
                                      int const n_groups,
                                      int const num_threads = 0) {
    int const n_cols = values.size();
    int const n_districts = district_start.size() - 1;

    if (n_districts < 1) {
        throw Rcpp::exception("compute_rhats_cpp: need at least one district.");
    }
    if (row_index.size() != group_index.size()) {
        throw Rcpp::exception(
            "compute_rhats_cpp: row_index and group_index must be the same length."
        );
    }
    if (n_groups < 1) {
        throw Rcpp::exception("compute_rhats_cpp: need at least one group.");
    }

    // Hold the columns so their data pointers stay valid, then read through
    // raw pointers inside the parallel region. No R API is touched there.
    std::vector<Rcpp::NumericVector> cols;
    cols.reserve(n_cols);
    for (int c = 0; c < n_cols; ++c) cols.push_back(values[c]);

    std::vector<double const *> col_ptr(n_cols);
    for (int c = 0; c < n_cols; ++c) col_ptr[c] = cols[c].begin();

    /*
     * Group membership depends only on the district, not on the statistic, so
     * build it once here rather than inside every task. Indices are local to
     * the district, matching the gathered buffers below.
     *
     * Empty groups are dropped so that the group count and mean group size
     * match what R's `split()` produces, which only creates levels that are
     * actually present.
     */
    std::vector<std::vector<std::vector<int>>> district_groups(n_districts);
    std::vector<double> district_n_g(n_districts, NA_REAL);

    for (int d = 0; d < n_districts; ++d) {
        int const from = district_start[d];
        int const to = district_start[d + 1];

        std::vector<std::vector<int>> raw(n_groups);
        for (int i = from; i < to; ++i) {
            raw[group_index[i]].push_back(i - from);
        }

        auto &kept = district_groups[d];
        long double total = 0.0L;
        for (int g = 0; g < n_groups; ++g) {
            if (!raw[g].empty()) {
                total += raw[g].size();
                kept.push_back(std::move(raw[g]));
            }
        }

        if (!kept.empty()) {
            district_n_g[d] = static_cast<double>(total / kept.size());
        }
    }

    Rcpp::NumericMatrix out(n_districts, n_cols);
    double *out_ptr = out.begin();

    int const n_tasks = n_cols * n_districts;

    RcppThread::parallelFor(
        0, n_tasks,
        [&](int const task) {
            int const c = task / n_districts;
            int const d = task % n_districts;

            auto const &groups = district_groups[d];
            if (groups.empty()) {
                out_ptr[c * n_districts + d] = NA_REAL;
                return;
            }

            int const from = district_start[d];
            int const to = district_start[d + 1];
            int const n = to - from;

            double const *column = col_ptr[c];

            // Buffers are local to the task, so nothing is shared.
            std::vector<double> buf(n);
            for (int i = 0; i < n; ++i) buf[i] = column[row_index[from + i]];

            std::vector<double> z;
            std::vector<std::pair<double, int>> sort_buf;
            std::vector<double> means, vars;

            ranknorm(buf, z, sort_buf);
            double const plain = calc_rhat(z, groups, district_n_g[d], means, vars);

            double const med = r_median(buf);
            std::vector<double> folded(n);
            for (int i = 0; i < n; ++i) folded[i] = std::fabs(buf[i] - med);

            ranknorm(folded, z, sort_buf);
            double const folded_rhat =
                calc_rhat(z, groups, district_n_g[d], means, vars);

            out_ptr[c * n_districts + d] = std::max(plain, folded_rhat);
        },
        resolve_num_threads(num_threads));

    return out;
}
