/********************************************************
 * Author: Cory McCartan
 * Institution: Harvard University
 * Date Created: 2021/02
 * Purpose: Merge-split MCMC redistricting sampler
 * (like Carter et al. 2019 but using the SMC proposal)
 ********************************************************/

#include "merge_split.h"

/*
 * Main entry point.
 *
 * USING MCMMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
Rcpp::List ms_plans(int N, List l, const uvec init, const uvec &counties, const uvec &pop,
              int n_distr, double target, double lower, double upper, double rho,
              double beta_sq, const uvec &current, int n_current,
              double beta_vra, double tgt_min, double tgt_other,
              double pow_vra, const uvec &min_pop,
              double beta_vra_hinge, const vec &tgts_min,
              double beta_inc, const uvec &incumbents, double beta_splits,
              double beta_fractures, double thresh, int k, int verbosity) {
    // re-seed MT
    generator.seed((int) Rcpp::sample(INT_MAX, 1)[0]);

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    int n_cty = max(counties);

    umat districts(V, N, fill::zeros);
    districts.col(0) = init;

    Rcpp::IntegerVector mh_decisions(N - 1);

    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout << "MARKOV CHAIN MONTE CARLO\n";
        Rcout << "Sampling " << N-1 << " " << V << "-unit maps with " << n_distr
              << " districts and population between " << lower << " and " << upper << ".\n";
        if (cg.size() > 1)
            Rcout << "Sampling hierarchically with respect to the "
                  << cg.size() << " administrative units.\n";
    }

    // find k and multipliers
    if (k <= 0) {
        adapt_ms_parameters(g, n_distr, k, thresh, tol, init, counties, cg, pop, target);
    }
    if (verbosity >= 2)
        Rcout << "Using k = " << k << "\n";

    int distr_1, distr_2;
    select_pair(n_distr, g, init, distr_1, distr_2);
    int refresh = std::max(N / 20, 1);
    int n_accept = 0;
    int reject_ct;
    for (int i = 1; i < N; i++) {
        districts.col(i) = districts.col(i - 1); // copy over old map

        // make the proposal
        double prop_lp = 0.0;
        reject_ct = 0;
        do {
            select_pair(n_distr, g, districts.col(i), distr_1, distr_2);
            prop_lp = split_map_ms(g, counties, cg, districts.col(i), distr_1,
                                   distr_2, pop, lower, upper, target, k);
            if (reject_ct % 200 == 0) Rcpp::checkUserInterrupt();
            reject_ct++;
        } while (!std::isfinite(prop_lp));

        // tau calculations
        if (rho != 1) {
            double log_st = 0;
            for (int j = 1; j <= n_cty; j++) {
                log_st += log_st_distr(g, districts, counties, i-1, distr_1, j);
                log_st += log_st_distr(g, districts, counties, i-1, distr_2, j);
                log_st -= log_st_distr(g, districts, counties, i, distr_1, j);
                log_st -= log_st_distr(g, districts, counties, i, distr_2, j);
            }
            log_st += log_st_contr(g, districts, counties, n_cty, i-1, distr_1);
            log_st += log_st_contr(g, districts, counties, n_cty, i-1, distr_2);
            log_st -= log_st_contr(g, districts, counties, n_cty, i, distr_1);
            log_st -= log_st_contr(g, districts, counties, n_cty, i, distr_2);

            prop_lp += (1 - rho) * log_st;
        }

        // add gibbs target
        // NOTE: different signs than above b/c of how Metropolis proposal has
        // transition ratio flipped relative to the target density ratio
        prop_lp -= calc_gibbs_tgt(districts.col(i), n_distr, V, distr_1, distr_2,
                                  pop, beta_sq, current, n_current, beta_vra,
                                  tgt_min, tgt_other, pow_vra, min_pop,
                                  beta_vra_hinge, tgts_min,
                                  beta_inc, incumbents,
                                  beta_splits, beta_fractures, counties, n_cty);
        prop_lp += calc_gibbs_tgt(districts.col(i-1), n_distr, V, distr_1, distr_2,
                                  pop, beta_sq, current, n_current, beta_vra,
                                  tgt_min, tgt_other, pow_vra, min_pop,
                                  beta_vra_hinge, tgts_min,
                                  beta_inc, incumbents,
                                  beta_splits, beta_fractures, counties, n_cty);

        double alpha = exp(prop_lp);
        if (alpha >= 1 || unif(generator) <= alpha) { // ACCEPT
            n_accept++;
            // map already stored in districts.col(i);
            mh_decisions(i - 1) = 1;
        } else { // REJECT
            districts.col(i) = districts.col(i - 1); // copy over old map
            mh_decisions(i - 1) = 0;
        }

        if (verbosity >= 2 && refresh > 0 && (i+1) % refresh == 0) {
            Rcout << "Iteration " << i+1 << "/" << N-1 << std::endl;
        }
        Rcpp::checkUserInterrupt();
    }

    if (verbosity >= 1) {
        Rcout << "Acceptance rate: " << (100.0 * n_accept) / (N-1) << std::endl;
    }

    Rcpp::List out;
    out["plans"] = districts;
    out["mhdecisions"] = mh_decisions;

    return out;
}


/*
 * Add specific constraint weights & return the cumulative weight vector
 */
double calc_gibbs_tgt(const subview_col<uword> &plan, int n_distr, int V,
                      int distr_1, int distr_2, const uvec &pop, double beta_sq,
                      const uvec &current, int n_current,
                      double beta_vra, double tgt_min, double tgt_other,
                      double pow_vra, const uvec &min_pop,
                      double beta_vra_hinge, const vec &tgts_min,
                      double beta_inc, const uvec &incumbents,
                      double beta_splits, double beta_fractures,
                      const uvec &counties, int n_cty) {
    double log_tgt = 0;

    if (beta_sq != 0)
        log_tgt += beta_sq * (
            sq_entropy(plan, current, distr_1, pop, n_distr, n_current, V) +
            sq_entropy(plan, current, distr_2, pop, n_distr, n_current, V)
        );
    if (beta_vra != 0)
        log_tgt += beta_vra * (
            eval_vra(plan, distr_1, tgt_min, tgt_other, pow_vra, pop, min_pop) +
            eval_vra(plan, distr_2, tgt_min, tgt_other, pow_vra, pop, min_pop)
        );
    if (beta_vra_hinge != 0)
        log_tgt += beta_vra_hinge * (
            eval_vra_hinge(plan, distr_1, tgts_min, pop, min_pop) +
            eval_vra_hinge(plan, distr_2, tgts_min, pop, min_pop)
        );
    if (beta_inc != 0)
        log_tgt += beta_inc * (
            eval_inc(plan, distr_1, incumbents) +
            eval_inc(plan, distr_2, incumbents)
        );
    if (beta_splits != 0)
        log_tgt += beta_splits * (
            eval_splits(plan, distr_1, counties, n_cty) +
            eval_splits(plan, distr_2, counties, n_cty)
        );
    if (beta_fractures != 0)
        log_tgt += beta_fractures * (
            eval_fractures(plan, distr_1, counties, n_cty) +
                eval_fractures(plan, distr_2, counties, n_cty)
        );

    return log_tgt;
}


/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map_ms(const Graph &g, const uvec &counties, Multigraph &cg,
                    subview_col<uword> districts, int distr_1, int distr_2,
                    const uvec &pop, double lower, double upper, double target,
                    int k) {
    int V = g.size();
    double orig_lb = log_boundary(g, districts, distr_1, distr_2);

    Tree ust = init_tree(V);
    std::vector<bool> ignore(V);
    double total_pop = 0;
    for (int i = 0; i < V; i++) {
        if (districts(i) == distr_1 || districts(i) == distr_2) {
            total_pop += pop(i);
            ignore[i] = false;
        } else {
            ignore[i] = true;
        }
    }

    int root;
    ust = sample_sub_ust(g, ust, V, root, ignore, pop, lower, upper, counties, cg);
    if (ust.size() == 0) return -log(0.0);

    // set `lower` as a way to return population of new district
    bool success = cut_districts_ms(ust, k, root, districts, distr_1, distr_2,
                                    pop, total_pop, lower, upper, target);

    if (!success) return -log(0.0); // reject sample

    return orig_lb - log_boundary(g, districts, distr_1, distr_2);
}


/*
 * Cut district into two pieces of roughly equal population
 */
// TESTED
bool cut_districts_ms(Tree &ust, int k, int root, subview_col<uword> &districts,
                      int distr_1, int distr_2, const uvec &pop, double total_pop,
                      double lower, double upper, double target) {
    int V = ust.size();
    // in case we pick a small-V district
    k = std::max(std::min(k, V-3), 1);
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
        candidates.push_back(i);
        deviances.push_back(std::max(dev1, dev2));
        is_ok.push_back(lower <= below && below <= upper &&
            lower <= total_pop - below && total_pop - below <= upper);
    }
    if ((int) candidates.size() < k) return false;

    int idx = rint(k);
    idx = select_k(deviances, idx + 1);
    int cut_at = candidates[idx];
    // reject sample
    if (!is_ok[idx]) return false;

    // find index of node to cut at
    std::vector<int> *siblings = &ust[parent[cut_at]];
    int length = siblings->size();
    int j;
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_at) break;
    }

    siblings->erase(siblings->begin()+j); // remove edge
    parent[cut_at] = -1;

    if (distr_root == distr_1) {
        assign_district(ust, districts, root, distr_1);
        assign_district(ust, districts, cut_at, distr_2);
    } else {
        assign_district(ust, districts, root, distr_2);
        assign_district(ust, districts, cut_at, distr_1);
    }

    return true;
}


/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void adapt_ms_parameters(const Graph &g, int n_distr, int &k, double thresh,
                         double tol, const uvec &plan, const uvec &counties,
                         Multigraph &cg, const uvec &pop, double target) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(20 + ((int) std::sqrt(V)), V - 1); // heuristic
    int N_adapt = (int) std::floor(4000.0 / sqrt((double) V));

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);

    std::vector<std::vector<double>> devs;
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    int distr_1, distr_2;
    for (int i = 0; i < N_adapt; i++) {
        Tree ust = init_tree(V);

        double joint_pop = 0;
        select_pair(n_distr, g, plan, distr_1, distr_2);
        for (int j = 0; j < V; j++) {
            if (plan(j) == distr_1 || plan(j) == distr_2) {
                joint_pop += pop(j);
                ignore[j] = false;
            } else {
                ignore[j] = true;
            }
        }

        ust = sample_sub_ust(g, ust, V, root, ignore, pop, lower, upper, counties, cg);
        if (ust.size() == 0) {
            i--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, joint_pop, target));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            if (ignore[j]) devs.at(i).at(j) = 2; // force not to work
            n_ok += devs.at(i).at(j) <= tol;
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
        Rcerr << "Warning: maximum hit; falling back to naive k estimator.\n";
        k = max_ok + 1;
    }
}

/*
 * Select a pair of neighboring districts i, j
 */
void select_pair(int n, const Graph &g, const uvec &plan, int &i, int &j) {
    int V = g.size();
    i = 1 + rint(n);

    std::set<int> neighboring;
    for (int k = 0; k < V; k++) {
        if (plan(k) != i) continue;
        std::vector<int> nbors = g[k];
        int length = nbors.size();
        for (int l = 0; l < length; l++) {
            int nbor = nbors[l];
            if (plan(nbor) == i) continue;
            neighboring.insert(plan[nbor]);
        }
    }

    int n_nbor = neighboring.size();
    j = *std::next(neighboring.begin(), rint(n_nbor));

    return;
}
