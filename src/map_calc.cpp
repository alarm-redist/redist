#include "map_calc.h"

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_boundary(const Graph &g, const subview_col<uword> &districts,
                    int distr_root, int distr_other) {
    int V = g.size();

    int count = 0; // number of cuttable edges to create eq-pop districts
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        if (districts(i) != distr_root) continue; // same side of boundary as root
        for (int j = 0; j < length; j++) {
            int nbor = nbors[j];
            if (districts(nbor) != distr_other)
                continue;
            // otherwise, boundary with root -> ... -> i -> nbor
            count++;
        }
    }

    return log((double) count);
}

/*
 * Compute the status quo penalty for district `distr`
 */
double sq_entropy(const subview_col<uword> &districts, const uvec &current,
                  int distr, const uvec &pop, int n_distr, int n_current, int V) {
    double accuml = 0;
    for (int j = 1; j <= n_current; j++) { // 1-indexed districts
        double pop_overlap = 0;
        double pop_total = 0;
        for (int k = 0; k < V; k++) {
            if (current[k] != j) continue;
            pop_total += pop[k];

            if (districts[k] == distr)
                pop_overlap += pop[k];
        }
        double frac = pop_overlap / pop_total;
        if (frac > 0)
            accuml += frac * std::log(frac);
    }

    return -accuml / n_distr / std::log(n_current);
}

/*
 * Compute the VRA penalty for district `distr`
 */
double eval_vra(const subview_col<uword> &districts, int distr, double tgt_min,
                double tgt_other, double pow_vra, const uvec &pop, const uvec &min_pop) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(min_pop(idxs))) / sum(pop(idxs));
    return std::pow(std::abs(frac - tgt_min), pow_vra) *
        std::pow(std::abs(frac - tgt_other), pow_vra);
}

/*
 * Compute the incumbent-preserving penalty for district `distr`
 */
double eval_inc(const subview_col<uword> &districts, int distr, const uvec &incumbents) {
    int n_inc = incumbents.size();
    double inc_in_distr = -1.0; // first incumbent doesn't count
    for (int i = 0; i < n_inc; i++) {
        if (districts[incumbents[i]] == distr)
            inc_in_distr++;
    }

    return inc_in_distr;
}


/*
 * Compute the deviation from the equal population constraint.
 */
// TESTED
NumericMatrix pop_dev(const umat &districts, const uvec &pop, int n_distr) {
    int N = districts.n_cols;
    int V = districts.n_rows;
    double target_pop = sum(pop) / n_distr;
    NumericMatrix dev(n_distr, N);
    for (int i = 0; i < N; i++) {
        std::vector<double> accuml(n_distr);
        for (int j = 0; j < V; j++) {
            int d = districts(j, i) - 1; // districts are 1-indexed
            accuml.at(d) += pop(j) / target_pop;
        }
        for (int d = 0; d < n_distr; d++) {
            dev(d,i) = abs(accuml[d]-1);
        }
    }
    return dev;
}

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// TESTED
NumericVector max_dev(const umat &districts, const uvec &pop, int n_distr) {
    int N = districts.n_cols;
    NumericMatrix dev = pop_dev(districts, pop, n_distr);
    NumericVector max_d(N);
    for (int i = 0; i < N; i++) {
        max_d(i) = max(dev(_, i));
    }
    return max_d;
}

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
std::vector<double> tree_dev(Tree &ust, int root, const uvec &pop,
                             double total_pop, double target) {
    int V = pop.size();
    std::vector<int> pop_below(V);
    std::vector<int> parent(V);
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of candidate edges to cut
    int idx = 0;
    std::vector<double> devs(V-1);
    for (int i = 0; i < V; i++) {
        if (i == root) continue;
        devs.at(idx) = std::min(std::abs(pop_below.at(i) - target),
                std::abs(total_pop - pop_below[i] - target)) / target;
        idx++;
    }

    std::sort(devs.begin(), devs.end());

    return devs;
}
