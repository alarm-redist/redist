// [[Rcpp::depends(redistmetrics)]]
#include "map_calc.h"
#include <redistmetrics.h>

/*
 * Compute the logarithm of the graph theoretic length of the boundary between
 * `distr_root` and `distr_other`, where the root of `ust` is in `distr_root`
 */
double log_boundary(const Graph &g, const subview_col<uword> &districts,
                    int distr_root, int distr_other) {
    int V = g.size();

    double count = 0; // number of cuttable edges to create eq-pop districts
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        if (districts(i) != distr_root) continue; // same side of boundary as root
        for (int nbor : nbors) {
            if (districts(nbor) != distr_other)
                continue;
            // otherwise, boundary with root -> ... -> i -> nbor
            count += 1.0;
        }
    }

    return std::log(count);
}

/*
 * Compute the status quo penalty for district `distr`
 */
double eval_sq_entropy(const subview_col<uword> &districts, const uvec &current,
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
 * Compute the new, hinge group penalty for district `distr`
 */
double eval_grp_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }

    return std::sqrt(std::max(0.0, target - frac));
}

/*
 * Compute the new, hinge group penalty for district `distr`
 */
double eval_grp_inv_hinge(const subview_col<uword> &districts, int distr,
                      const vec &tgts_grp, const uvec &grp_pop, const uvec &total_pop) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }

    return std::sqrt(std::max(0.0, frac - target));
}

/*
 * Compute the power-based group penalty for district `distr`
 */
double eval_grp_pow(const subview_col<uword> &districts, int distr,
                   const uvec &grp_pop, const uvec &total_pop,
                   double tgt_grp, double tgt_other, double pow) {
    uvec idxs = find(districts == distr);
    double frac = ((double) sum(grp_pop(idxs))) / sum(total_pop(idxs));
    return std::pow(std::fabs(frac - tgt_grp) * std::fabs(frac - tgt_other), pow);
}

/*
 * Compute the incumbent-preserving penalty for district `distr`
 */
double eval_inc(const subview_col<uword> &districts, int distr, const uvec &incumbents) {
    int n_inc = incumbents.size();
    double inc_in_distr = -1.0; // first incumbent doesn't count
    for (int i = 0; i < n_inc; i++) {
        if (districts[incumbents[i] - 1] == distr)
            inc_in_distr++;
    }

    if (inc_in_distr < 0.0) {
        inc_in_distr = 0.0;
    }

    return inc_in_distr;
}


// helper function
// calculates districts which appear in each county (but not zeros)
std::vector<std::set<int>> calc_county_dist(const subview_col<uword> &districts,
                                            const uvec &counties, int n_cty,
                                            bool zero_ok) {
    std::vector<std::set<int>> county_dist(n_cty);
    int V = counties.size();
    for (int i = 0; i < n_cty; i++) {
        county_dist[i] = std::set<int>();
    }
    for (int i = 0; i < V; i++) {
        if (zero_ok || districts[i] > 0) {
            county_dist[counties[i]-1].insert(districts[i]);
        }
    }
    return county_dist;
}

/*
 * Compute the county split penalty for district `distr`
 */
double eval_splits(const subview_col<uword> &districts, int distr,
                   const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    int splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 2 : cty_n_distr >= 2;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the county multisplit penalty for district `distr`
 */
double eval_multisplits(const subview_col<uword> &districts, int distr,
                        const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    double splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 3 : cty_n_distr >= 3;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the total splits penalty for district `distr`
 */
double eval_total_splits(const subview_col<uword> &districts, int distr,
                         const uvec &counties, int n_cty, bool smc) {
    std::vector<std::set<int>> county_dist = calc_county_dist(districts, counties, n_cty, distr == 0);

    double splits = 0;
    for (int i = 0; i < n_cty; i++) {
        int cty_n_distr = county_dist[i].size();
        // no over-counting since every split counts
        if (cty_n_distr > 1) {
            if (smc) {
                auto search = county_dist[i].find(distr);
                if (search != county_dist[i].end()) {
                    splits++;
                }
            } else {
                splits++;
            }
        }
    }

    return splits;
}

/*
 * Compute the Polsby Popper penalty for district `distr`
 */
double eval_polsby(const subview_col<uword> &districts, int distr,
                   const ivec &from,
                   const ivec &to,
                   const vec &area,
                   const vec &perimeter) {
    uvec idxs = find(districts == distr);

    double pi4 = 4.0 * 3.14159265;
    double tot_area = sum(area(idxs));

    double tot_perim = 0.0;

    for (int e = 0; e < from.size(); e++) {
        if(from(e) == -1) {
            if(districts(to(e) - 1) == distr) {
                tot_perim += perimeter(e);
            }
        } else {
            if(districts(from(e) - 1) == distr && districts(to(e) - 1) != distr) {
                tot_perim += perimeter(e);
            }
        }
    }

    double dist_peri2 = pow(tot_perim, 2.0);
    return 1.0 - (pi4 * tot_area / dist_peri2);
}

/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const subview_col<uword> &districts, int distr,
                     const uvec &total_pop, mat ssdmat, double denominator = 1.0) {
    uvec idxs = find(districts == distr);
    double ssd = 0.0;

    for (int i = 0; i < idxs.size() - 1; i++) {
        for (int k = i + 1; k < idxs.size(); k++) {
            ssd += (double) ssdmat(idxs(i), idxs(k)) * total_pop(idxs(i)) *
                total_pop(idxs(k));
        }
    }

    return ssd / denominator;
}

/*
 * Compute the population penalty for district `distr`
 */
double eval_pop_dev(const subview_col<uword> &districts, int distr,
                       const uvec &total_pop, double parity) {
    uvec idxs = find(districts == distr);
    double pop = sum(total_pop(idxs));

    return std::pow(pop / parity - 1.0, 2.0);
}


/*
 * Compute the segregation penalty for district `distr`
 */
double eval_segregation(const subview_col<uword> &districts, int distr,
                        const uvec &grp_pop, const uvec &total_pop) {

    int T = sum(total_pop);
    double pAll = (double) sum(grp_pop) / T;
    double denom = (double) 2.0 * T * pAll * (1 - pAll);

    uvec idxs = find(districts == distr);
    double grp = sum(grp_pop(idxs));
    double pop = sum(total_pop(idxs));

    return (double)(pop * std::abs((grp / pop) - pAll) / denom);
}

/*
 * Compute the qps penalty for district `distr`
 */
double eval_qps(const subview_col<uword> &districts, int distr,
                const uvec &total_pop, const uvec &cities, int n_city,
                int nd) {

    vec tally(n_city);
    vec pj(n_city);
    vec j(n_city);
    vec sumpj(n_city);

    uvec idxs_d = find(districts == distr);
    double pop = sum(total_pop(idxs_d));

    for (int i = 0; i < n_city; i++){
        uvec idxs = find(cities == (i + 1));
        idxs = arma::intersect(idxs_d, idxs);
        tally(i) = sum(total_pop(idxs));
        if (tally(i) > 0) {
            j(i) += 1;
        }
    }

    pj = tally / pop;
    sumpj = pj * (1.0 -  pj);
    sumpj = sumpj / (double) nd;

    return sum(sumpj) + log(sum(j));
}

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_log_st(const subview_col<uword> &districts, const Graph g,
                   arma::uvec counties, int ndists) {
    return (double) redistmetrics::log_st_map(g, districts, counties, ndists)[0];
}

/*
 * Compute the edges removed penalty for district `distr`
 */
double eval_er(const subview_col<uword> &districts, const Graph g, int ndists) {
    return (double) redistmetrics::n_removed(g, districts, ndists)[0];
}




/*
 * Compute the cooccurence matrix for a set of precincts indexed by `idxs`,
 * given a collection of plans
 */
mat prec_cooccur(umat m, uvec idxs, int ncores) {
    int v = m.n_rows;
    int n = idxs.n_elem;
    mat out(v, v);

    RcppThread::parallelFor(0, v, [&] (int i) {
        out(i, i) = 1;
        for (int j = 0; j < i; j++) {
            double shared = 0;
            for (int k = 0; k < n; k++) {
                shared += m(i, idxs[k]-1) == m(j, idxs[k]-1);
            }
            shared /= n;
            out(i, j) = shared;
            out(j, i) = shared;
        }
    }, ncores);

    return out;
}

/*
 * Compute the percentage of `group` in each district. Asummes `m` is 1-indexed.
 */
NumericMatrix group_pct(umat m, vec group_pop, vec total_pop, int n_distr) {
    int v = m.n_rows;
    int n = m.n_cols;

    NumericMatrix grp_distr(n_distr, n);
    NumericMatrix tot_distr(n_distr, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < v; j++) {
            int distr = m(j, i) - 1;
            grp_distr(distr, i) += group_pop[j];
            tot_distr(distr, i) += total_pop[j];
        }
    }

    // divide
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n_distr; j++) {
            grp_distr(j, i) /= tot_distr(j, i);
        }
    }

    return grp_distr;
}

/*
 * Compute the percentage of `group` in each district, and return the `k`-th
 * largest such value. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
NumericVector group_pct_top_k(const IntegerMatrix m, const NumericVector group_pop,
                              const NumericVector total_pop, int k, int n_distr) {
    int v = m.nrow();
    int n = m.ncol();
    NumericVector out(n);

    for (int i = 0; i < n; i++) {
         std::vector<double> grp_distr(n_distr, 0.0);
         std::vector<double> tot_distr(n_distr, 0.0);

        for (int j = 0; j < v; j++) {
            int distr = m(j, i) - 1;
            grp_distr[distr] += group_pop[j];
            tot_distr[distr] += total_pop[j];
        }

        for (int j = 0; j < n_distr; j++) {
            grp_distr[j] /= tot_distr[j];
        }

        std::nth_element(grp_distr.begin(), grp_distr.begin() + k - 1,
                         grp_distr.end(), std::greater<double>());

        out[i] = grp_distr[k - 1];
    }

    return out;
}


/*
 * Tally a variable by district.
 */
// TESTED
NumericMatrix pop_tally(IntegerMatrix districts, vec pop, int n_distr) {
    int N = districts.ncol();
    int V = districts.nrow();

    NumericMatrix tally(n_distr, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < V; j++) {
            int d = districts(j, i) - 1; // districts are 1-indexed
            tally(d, i) = tally(d, i) + pop(j);
        }
    }

    return tally;
}

/*
 * Create the projective distribution of a variable `x`
 */
// [[Rcpp::export]]
NumericMatrix proj_distr_m(IntegerMatrix districts, const arma::vec x,
                           IntegerVector draw_idx, int n_distr) {
    int n = draw_idx.size();
    int V = districts.nrow();

    NumericMatrix out(V, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < V; j++) {
            int idx = draw_idx[i] - 1;
            out(j, i) = x[n_distr*idx + districts(j, idx) - 1];
        }
    }

    return out;
}

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// TESTED
NumericVector max_dev(const IntegerMatrix districts, const vec pop, int n_distr) {
    int N = districts.ncol();
    double target_pop = sum(pop) / n_distr;

    NumericMatrix dev = pop_tally(districts, pop, n_distr) / target_pop - 1.0;
    NumericVector res(N);
    for (int j = 0; j < n_distr; j++) {
        for (int i = 0; i < N; i++) {
            if (std::fabs(dev(j, i)) > res(i))
                res(i) = std::fabs(dev(j, i));
        }
    }

    return res;
}

/*
 * Calculate the deviation for cutting at every edge in a spanning tree.
 */
std::vector<double> tree_dev(Tree &ust, int root, const uvec &pop,
                             double total_pop, double target) {
    int V = pop.size();
    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of candidate edges to cut
    int idx = 0;
    std::vector<double> devs(V-1);
    for (int i = 0; i < V; i++) {
        if (i == root) continue;
        devs.at(idx) = std::min(std::fabs(pop_below.at(i) - target),
                std::fabs(total_pop - pop_below[i] - target)) / target;
        idx++;
    }

    std::sort(devs.begin(), devs.end());

    return devs;
}


/*
 * Column-wise maximum
 */
// [[Rcpp::export]]
NumericVector colmax(const NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(ncol);
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
NumericVector colmin(const NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(ncol);
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
