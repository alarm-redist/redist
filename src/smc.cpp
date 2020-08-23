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
 * population deviation is within `tol`
 */
IntegerMatrix smc_plans(int N, List l, const uvec &counties,
                        const uvec &pop, int n_distr, double tol,
                        double gamma, NumericVector &log_prob, double thresh,
                        double alpha, int infl, int verbosity) {
    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int N_max = infl * N;
    umat districts(V, N_max, fill::zeros);
    double total_pop = sum(pop);
    double distr_pop = total_pop / n_distr;

    if (verbosity >= 1) {
        Rcout << "Sampling " << N << " " << V << "-unit maps with " << n_distr
              << " districts and population tolerance " << tol*100 << "%.\n";
        if (cg.size() > 1)
            Rcout << "Ensuring no more than " << n_distr - 1 << " splits of the "
                  << cg.size() << " administrative units.\n";
    }

    vec pop_left(N_max);
    pop_left.fill(total_pop);
    vec lp(N_max, fill::zeros);

    int k;
    int N_adapt = std::min((int) std::floor(4000.0 / sqrt(V)), N_max);
    int N_sample = 0;
    int valid = N_adapt;
    double prob;
    for (int ctr = 1; ctr < n_distr; ctr++) {
        if (verbosity >= 1)
            Rcout << "Making split " << ctr << " of " << n_distr-1 << "\n";

        // find k and multipliers
        adapt_parameters(g, k, prob, N_adapt, valid, lp, thresh, tol, districts,
                         counties, cg, pop, pop_left, distr_pop);

        int N_new = std::ceil(N / prob);
        if (ctr == n_distr - 1) // safety margin for last step
            N_new *= n_distr/2;
        N_new = std::min(N_new, N_max);

        if (verbosity >= 3)
            Rprintf("Using k=%d for estimated success probability of %.2f%%\n",
                   k, 100.0*prob);

        // Resample
        resample_maps(N_sample, N_new, alpha, districts, lp, pop_left);
        N_sample = N_new;

        // perform sampling
        split_maps(g, counties, cg, pop, districts, lp, pop_left, N_sample,
                   n_distr, ctr, distr_pop, tol, gamma, k, verbosity);
        valid = 0;
        for (int i = 0; i < N_sample; i++) {
            valid += std::isfinite(lp(i));
        }

        if (valid == 0) {
            Rcout << "No valid samples at stage " << ctr << "; stopping sampling.\n"
            << "Consider increasing `N` and/or `infl`.\n";
            IntegerMatrix distr_ret(V, N);
            return distr_ret;
        }

        if (verbosity >= 2)
            Rprintf("%d valid samples (%.1f%% of total)\n", valid, (100.0*valid)/N_sample);

        Rcpp::checkUserInterrupt();
    }

    resample_maps(N_sample, valid, alpha, districts, lp, pop_left);

    // Output districts (and set final district label to n_distr rather than 0)
    N = std::min(N, valid);
    IntegerMatrix distr_ret(V, N);
    int idx = 0;
    for (int i = 0; i < N; i++) {
        if (std::isinf(lp(i))) continue;
        log_prob(idx) = lp(i);
        for (int j = 0; j < V; j++) {
            distr_ret(j, idx) = districts(j, i) == 0 ? n_distr : districts(j, i);
        }
        idx++;
    }

    return distr_ret;
}

/*
 * Resample partially-drawn maps according to their weights.
 */
void resample_maps(int N_sample, int N_new, double alpha, umat &districts,
                   vec &lp, vec &pop_left) {
    if (N_sample == 0) return;

    vec wgt = exp(-alpha * lp.subvec(0, N_sample - 1));
    lp = lp * (1 - alpha);

    uvec idxs = as<uvec>(Rcpp::sample(N_sample, N_new, true, wrap(wgt), false));

    lp = lp.elem(idxs);
    pop_left = pop_left.elem(idxs);
    districts = districts.cols(idxs);
}


/*
 * Split off a piece from each map in `districts`, keeping deviation within `tol`
 */
void split_maps(const Graph &g, const uvec &counties, Multigraph &cg,
                const uvec &pop, umat &districts, vec &lp, vec &pop_left,
                int N, int n_distr, int dist_ctr, double distr_pop, double tol,
                double gamma, int k, int verbosity) {
    // absolute bounds for district populations
    double upper = distr_pop * (1 + tol);
    double lower = distr_pop * (1 - tol);
    int new_size = n_distr - dist_ctr;
    int n_cty = max(counties);

    int refresh = N / 10; // how often to print update statements

    for (int i = 0; i < N; i++) {
        double lower_s = std::max(lower, pop_left(i) - new_size * upper);
        double upper_s = std::min(upper, pop_left(i) - new_size * lower);
        // split
        double inc_lp = split_map(g, counties, cg, districts.col(i), dist_ctr,
                          pop, pop_left(i), lower_s, upper_s, distr_pop, k);

        if (gamma != 1 && std::isfinite(inc_lp)) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, districts, counties, i, dist_ctr, j);
            }
            log_st += log_st_contr(g, districts, counties, n_cty, i, dist_ctr);

            if (dist_ctr == n_distr - 1) {
                for (int j = 1; j <= n_cty; j++) {
                    log_st += log_st_distr(g, districts, counties, i, 0, j);
                }
                log_st += log_st_contr(g, districts, counties, n_cty, i, 0);
            }

            inc_lp += (1 - gamma) * log_st;
        }
        lp(i) += inc_lp;

        // `lower_s` now contains the population of the newly-split district
        pop_left(i) -= lower_s;

        if (verbosity >= 2 && refresh > 0 && (i+1) % refresh == 0) {
            Rprintf("Iteration %'6d / %'d\n", i+1, N);
            Rcpp::checkUserInterrupt();
        }
    }

    // ensure the leftover slots don't get sampled
    int N_max = lp.size();
    for (int i = N; i < N_max; i++) {
        lp(i) = -log(0.0);
    }
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

    return log_boundary(g, districts, 0, dist_ctr) - log(k);
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
        double dev1 = std::abs(below - target);
        double dev2 = std::abs(total_pop - below - target);
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
    int cut_at = std::abs(candidates[idx]);
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
void adapt_parameters(const Graph &g, int &k, double &prob, int N_adapt, int valid,
                      const vec &lp, double thresh, double tol,
                      const umat &districts, const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      const vec &pop_left, double target) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(20 + ((int) std::sqrt(V)), V - 1); // heuristic
    int N_max = districts.n_cols;
    N_adapt = std::min(N_adapt, valid);

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

    if (k == k_max + 1) {
        Rcout << "Warning: maximum hit; falling back to naive k estimator.\n";
        k = max_ok + 1;
    }

    prob = 0;
    for (int i = 0; i <= k; i++) {
        prob += ((double) i / k) * distr_ok(i);
    }
    if (prob == 0.0) prob = 1.0 / N_max; // happens sometimes
}

