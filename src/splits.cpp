#include <cstddef>

#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <vector>

/*
 * Count, for each plan, how many administrative units are split across more
 * than `max_split` districts.
 *
 * `max_split = 1` counts every unit touched by 2 or more districts, and
 * `max_split = 2` counts only those touched by 3 or more (a "multisplit").
 * Both `dm` and `community` are 1-indexed.
 *
 * This is a reimplementation of redistmetrics::splits, which is the only reason
 * this file used to need RcppArmadillo. The one deliberate difference is the
 * bound on the number of units: redistmetrics uses the count of distinct ids,
 * which goes out of bounds if the ids are not contiguous, whereas the maximum
 * id is always a safe bound. The two agree for the contiguous 1..n ids the
 * callers produce with vctrs::vec_group_id().
 */
// [[Rcpp::export]]
IntegerVector splits(IntegerMatrix dm, IntegerVector community, int nd, int max_split) {
    int const n_plans = dm.ncol();
    int const V = dm.nrow();
    int const n_cty = community.size() == 0 ? 0 : Rcpp::max(community);

    IntegerVector out(n_plans);
    // seen[c * nd + d] marks that unit c contains a vertex of district d
    std::vector<char> seen(static_cast<std::size_t>(n_cty) * nd);

    for (int p = 0; p < n_plans; p++) {
        std::fill(seen.begin(), seen.end(), 0);

        for (int i = 0; i < V; i++) {
            std::size_t const cty = static_cast<std::size_t>(community[i] - 1);
            int const distr = dm(i, p) - 1;
            seen[cty * nd + distr] = 1;
        }

        int n_split = 0;
        for (std::size_t c = 0; c < static_cast<std::size_t>(n_cty); c++) {
            int n_distr_in_cty = 0;
            for (int d = 0; d < nd; d++) {
                n_distr_in_cty += seen[c * nd + d];
                if (n_distr_in_cty > max_split) {
                    n_split++;
                    break;
                }
            }
        }

        out[p] = n_split;
    }

    return out;
}

// [[Rcpp::export]]
IntegerMatrix dist_cty_splits(IntegerMatrix dm, IntegerVector community, int nd) {
    IntegerMatrix ret(nd, dm.ncol());
    IntegerVector com_name = sort_unique(community);
    IntegerVector com_found(com_name.size(), 0);

    // by column (aka map)
    for (int c = 0; c < dm.ncol(); c++) {
        // by district
        for (int d = 0; d < nd; d++) {
            com_found = IntegerVector(com_found.size(), 0);
            // across all rows
            for (int r = 0; r < dm.nrow(); r++) {
                if (dm(r, c) == d) {
                    com_found(community(r)) = 1;
                }
            }
            ret(d, c) = sum(com_found);
        }
    }
    return ret;
}
