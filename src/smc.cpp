/********************************************************
 * Author: Cory McCartan
 * Institution: Harvard University
 * Date Created: 2019/11
 * Purpose: Sequential Monte Carlo redistricting sampler
 ********************************************************/

#include "smc.h"

/*
 * Main entry point.
 *
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
umat smc_plans(int N, List l, const uvec &counties, const uvec &pop,
               int n_distr, double target, double lower, double upper, double rho,
               double beta_sq, const uvec &current, int n_current,
               double beta_vra, double tgt_min, double tgt_other,
               double pow_vra, const uvec &min_pop,
               double beta_vra_hinge, const vec &tgts_min,
               double beta_inc, const uvec &incumbents,
               vec &lp, double thresh,
               double alpha, double pop_temper, int verbosity) {
    // re-seed MT
    generator.seed((int) Rcpp::sample(INT_MAX, 1)[0]);

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    umat districts(V, N, fill::zeros);
    double total_pop = sum(pop);
    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout << "SEQUENTIAL MONTE CARLO\n";
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Ensuring no more than " << n_distr - 1 << " splits of the "
                  << cg.size() << " administrative units.\n";
    }

    vec pop_left(N);
    pop_left.fill(total_pop);
    vec log_temper(N);
    log_temper.fill(0.0);
    lp.fill(0.0);

    int k;
    vec cum_wgt(N);
    cum_wgt.fill(1.0 / N);
    cum_wgt = cumsum(cum_wgt);
    for (int ctr = 1; ctr < n_distr; ctr++) {
        if (verbosity >= 1)
            Rcout << "Making split " << ctr << " of " << n_distr-1 << "\n";

        // find k and multipliers
        adapt_parameters(g, k, lp, thresh, tol, districts, counties, cg, pop,
                         pop_left, target, verbosity);

        if (verbosity >= 3)
            Rcout << "Using k = " << k << "\n";

        // perform resampling/drawing
        split_maps(g, counties, cg, pop, districts, cum_wgt, lp, pop_left, log_temper,
                   pop_temper, n_distr, ctr, lower, upper, target, rho, k, verbosity);

        // compute weights for next step
        cum_wgt = get_wgts(districts, n_distr, ctr, alpha, lp, pop,
                           beta_sq, current, n_current,
                           beta_vra, tgt_min, tgt_other, pow_vra, min_pop,
                           beta_vra_hinge, tgts_min,
                           beta_inc, incumbents, verbosity);

        Rcpp::checkUserInterrupt();
    }

    lp = lp - log_temper;

    // Set final district label to n_distr rather than 0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < V; j++) {
            districts(j, i) = districts(j, i) == 0 ? n_distr : districts(j, i);
        }
    }

    return districts;
}


/*
 * Add specific constraint weights & return the cumulative weight vector
 */
vec get_wgts(const umat &districts, int n_distr, int distr_ctr,
             double alpha, vec &lp, const uvec &pop,
             double beta_sq, const uvec &current, int n_current,
             double beta_vra, double tgt_min, double tgt_other,
             double pow_vra, const uvec &min_pop,
             double beta_vra_hinge, const vec &tgts_min,
             double beta_inc, const uvec &incumbents, int verbosity) {
    int V = districts.n_rows;
    int N = districts.n_cols;

    for (int i = 0; i < N; i++) {
        if (beta_sq != 0)
            lp[i] += beta_sq * sq_entropy(districts.col(i), current, distr_ctr,
                                          pop, n_distr, n_current, V);
        if (beta_vra != 0)
            lp[i] += beta_vra * eval_vra(districts.col(i), distr_ctr, tgt_min,
                                         tgt_other, pow_vra, pop, min_pop);
        if (beta_vra_hinge != 0)
            lp[i] += beta_vra_hinge * eval_vra_hinge(districts.col(i), distr_ctr,
                                                     tgts_min, pop, min_pop);
        if (beta_inc != 0)
            lp[i] += beta_inc * eval_inc(districts.col(i), distr_ctr, incumbents);
    }

    vec wgt = exp(-alpha * lp);
    if (distr_ctr < n_distr - 1) // not the last iteration
        lp = lp * (1 - alpha);
    vec cuml_wgt = cumsum(wgt);

    if (verbosity >= 1) {
        double neff = cuml_wgt[N-1] * cuml_wgt[N-1]  / sum(square(wgt));
        Rprintf("Resampling effective sample size: %.1f (%.1f%% efficiency).\n", neff, 100*neff/N);
    }

    return cuml_wgt / cuml_wgt[N-1];
}


/*
 * Split off a piece from each map in `districts`,
 * keeping deviation between `lower` and `upper`
 */
void split_maps(const Graph &g, const uvec &counties, Multigraph &cg,
                const uvec &pop, umat &districts, vec &cum_wgt, vec &lp,
                vec &pop_left, vec &log_temper, double pop_temper, int n_distr,
                int dist_ctr, double lower, double upper, double target,
                double rho, int k, int verbosity) {
    int V = districts.n_rows;
    int N = districts.n_cols;
    int new_size = n_distr - dist_ctr;
    int n_cty = max(counties);

    umat districts_new(V, N);
    vec pop_left_new(N);
    vec lp_new(N);
    vec log_temper_new(N);

    int refresh = N / 10; // how often to print update statements
    double iter = 0; // how many actual iterations
    for (int i = 0; i < N; i++, iter++) {
        // resample
        int idx = rint(N, cum_wgt);
        districts_new.col(i) = districts.col(idx);
        double lower_s = std::max(lower, pop_left(idx) - new_size * upper);
        double upper_s = std::min(upper, pop_left(idx) - new_size * lower);

        // split
        double inc_lp = split_map(g, counties, cg, districts_new.col(i), dist_ctr,
                                  pop, pop_left(idx), lower_s, upper_s, target, k);

        // bad sample; try again
        if (!std::isfinite(inc_lp)) {
            i--;
            continue;
        }

        if (rho != 1) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, districts_new, counties, i, dist_ctr, j);
            }
            log_st += log_st_contr(g, districts_new, counties, n_cty, i, dist_ctr);

            if (dist_ctr == n_distr - 1) {
                for (int j = 1; j <= n_cty; j++) {
                    log_st += log_st_distr(g, districts_new, counties, i, 0, j);
                }
                log_st += log_st_contr(g, districts_new, counties, n_cty, i, 0);
            }

            inc_lp += (1 - rho) * log_st;
        }
        // `lower_s` now contains the population of the newly-split district
        pop_left_new(i) = pop_left(idx) - lower_s;
        double pop_pen = sqrt((double) n_distr - 2) * log(std::fabs(lower_s - target)/target);
        log_temper_new(i) = log_temper(idx) - pop_temper*pop_pen;

        lp_new(i) = lp(idx) + inc_lp - pop_temper*pop_pen;

        if (verbosity >= 2 && refresh > 0 && (i+1) % refresh == 0) {
            Rprintf("Iteration %'6d / %'d\n", i+1, N);
            Rcpp::checkUserInterrupt();
        }
    }
    if (verbosity >= 2) {
        Rprintf("%.1f%% acceptance rate.\n", 100.0 * N / iter);
    }

    districts = districts_new;
    pop_left = pop_left_new;
    lp = lp_new;
    log_temper = log_temper_new;
}

/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map(const Graph &g, const uvec &counties, Multigraph &cg,
                 subview_col<uword> districts, int dist_ctr, const uvec &pop,
                 double total_pop, double &lower, double upper, double target, int k) {
    int V = g.size();

    Tree ust = init_tree(V);
    std::vector<bool> ignore(V);
    for (int i = 0; i < V; i++) ignore[i] = districts(i) != 0;

    int root;
    ust = sample_sub_ust(g, ust, V, root, ignore, counties, cg);
    if (ust.size() == 0) return -log(0.0);

    // set `lower` as a way to return population of new district
    lower = cut_districts(ust, k, root, districts, dist_ctr, pop, total_pop,
                          lower, upper, target);
    if (lower == 0) return -log(0.0); // reject sample

    return log_boundary(g, districts, 0, dist_ctr) - log((double)k);
}


/*
 * Cut district into two pieces of roughly equal population
 */
// TESTED
double cut_districts(Tree &ust, int k, int root, subview_col<uword> &districts,
                     int dist_ctr, const uvec &pop, double total_pop,
                     double lower, double upper, double target) {
    int V = ust.size();
    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
    parent[root] = -1;
    tree_pop(ust, root, pop, pop_below, parent);
    // compile a list of:
    std::vector<int> candidates; // candidate edges to cut,
    std::vector<double> deviances; // how far from target pop.
    std::vector<bool> is_ok; // whether they meet constraints
    int distr_root = districts(root);
    for (int i = 0; i < V; i++) {
        if (districts(i) != distr_root || i == root) continue;
        double below = pop_below.at(i);
        double dev1 = std::fabs(below - target);
        double dev2 = std::fabs(total_pop - below - target);
        if (dev1 < dev2) {
            candidates.push_back(i);
            deviances.push_back(dev1);
            is_ok.push_back(lower <= below && below <= upper);
        } else {
            candidates.push_back(-i);
            deviances.push_back(dev2);
            is_ok.push_back(lower <= total_pop - below && total_pop - below <= upper);
        }
    }
    if ((int) candidates.size() < k) return 0.0;

    int idx = rint(k);
    idx = select_k(deviances, idx + 1);
    int cut_at = std::fabs(candidates[idx]);
    // reject sample
    if (!is_ok[idx]) return 0.0;

    // find index of node to cut at
    std::vector<int> *siblings = &ust[parent[cut_at]];
    int length = siblings->size();
    int j;
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_at) break;
    }

    siblings->erase(siblings->begin()+j); // remove edge
    parent[cut_at] = -1;

    if (candidates[idx] > 0) { // if the newly cut district is final
        assign_district(ust, districts, cut_at, dist_ctr);
        return pop_below.at(cut_at);
    } else { // if the root-side district is final
        assign_district(ust, districts, root, dist_ctr);
        return total_pop - pop_below.at(cut_at);
    }
}


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_parameters(const Graph &g, int &k, const vec &lp, double thresh,
                      double tol, const umat &districts, const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      const vec &pop_left, double target, int verbosity) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(20 + ((int) std::sqrt((double)V)), V - 1); // heuristic
    int N_max = districts.n_cols;
    int N_adapt = std::min((int) std::floor(4000.0 / sqrt((double)V)), N_max);

    std::vector<std::vector<double>> devs;
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    for (int i = 0, idx = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (std::isinf(lp(i))) { // skip if not valid
            idx--;
            continue;
        }

        Tree ust = init_tree(V);
        for (int j = 0; j < V; j++) ignore[j] = districts(j, i) != 0;
        ust = sample_sub_ust(g, ust, V, root, ignore, counties, cg);
        if (ust.size() == 0) {
            idx--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, pop_left(i), target));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            n_ok += devs.at(idx).at(j) <= tol;
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max)
            max_ok = n_ok;

        Rcpp::checkUserInterrupt();
    }

    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    uvec idxs(N_adapt);
    for (k = 1; k <= k_max; k++) {
        idxs = as<uvec>(Rcpp::sample(k, N_adapt, true, R_NilValue, false));
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            double dev = devs.at(i).at(idxs[i]);
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k >= k_max) {
        if (verbosity >= 1) {
            Rcout << "Note: maximum hit; falling back to naive k estimator.\n";
        }
        k = max_ok + 1;
    }
}

