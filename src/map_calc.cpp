#include "map_calc.h"

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_boundary(const List &g, const IntegerMatrix::Column &districts,
                    int distr_root, int distr_other) {
    int V = g.size();

    int count = 0; // number of cuttable edges to create eq-pop districts
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = as<std::vector<int>>(g[i]);
        int length = nbors.size();
        if (districts[i] != distr_root) continue; // same side of boundary as root
        for (int j = 0; j < length; j++) {
            int nbor = nbors[j];
            if (districts[nbor] != distr_other)
                continue;
            // otherwise, boundary with root -> ... -> i -> nbor
            count++;
        }
    }

    return log(count);
}


/*
 * Compute the deviation from the equal population constraint.
 */
// TESTED
NumericMatrix pop_dev(const IntegerMatrix &districts,
                      const NumericVector &pop, int n_distr) {
    int N = districts.ncol();
    int V = districts.nrow();
    double target_pop = sum(pop) / n_distr;
    NumericMatrix dev(n_distr, N);
    for (int i = 0; i < N; i++) {
        std::vector<double> accuml(n_distr);
        for (int j = 0; j < V; j++) {
            int d = districts(j, i) - 1; // districts are 1-indexed
            accuml[d] += pop[j] / target_pop;
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
NumericVector max_dev(const IntegerMatrix &districts,
                      const NumericVector &pop, int n_distr) {
    int N = districts.ncol();
    NumericMatrix dev = pop_dev(districts, pop, n_distr);
    NumericVector max_d(N);
    for (int i = 0; i < N; i++) {
        max_d[i] = max(dev(_, i));
    }
    return max_d;
}

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
void tree_dev(Tree &ust, int root, vec res,
              const IntegerVector &pop, double total_pop, double target) {
    int V = pop.size();
    std::vector<int> pop_below(V);
    std::vector<int> parent(V);
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of candidate edges to cut
    int idx = 0;
    for (int i = 0; i < V; i++) {
        if (i == root) continue;
        res[idx] = std::min(std::abs(pop_below[i] - target),
                            std::abs(total_pop - pop_below[i] - target)) / target;
        idx++;
    }

    std::sort(res.begin(), res.end());
}


/*
 * Compute the isoperimetric quotient for maps, given the area of each precinct
 * and a matrix of pariwise boundary lengths (with external boundary lengths
 * along the diagonal)
 */
NumericMatrix calc_isoper_quo(const IntegerMatrix &districts,
                              const NumericVector &areas,
                              const NumericMatrix &perims, int n_distr) {
    int N = districts.ncol();
    int V = districts.nrow();
    NumericMatrix comp(n_distr, N);
    for (int n = 0; n < N; n++) {
        std::vector<double> area(n_distr, 0.0);
        std::vector<double> perim(n_distr, 0.0);
        for (int i = 0; i < V; i++) {
            int distr = districts(i, n) - 1;
            area[distr] += areas[i];
            perim[distr] += perims(i, i);
            for (int j = 0; j < V; j++) {
                if (districts(j, n) - 1 != distr) {
                   perim[distr] += perims(i, j);
                }
            }
        }
        for (int k = 0; k < n_distr; k++) {
            comp(k, n) = (4 * M_PI * area[k] / perim[k]) / perim[k];
        }
    }

    return comp;
}
