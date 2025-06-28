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
List smc_plans(int N, List l, const uvec &counties, const uvec &pop,
               int n_distr, double target, double lower, double upper, double rho,
               umat districts, int n_drawn, int n_steps,
               List constraints, List control, int verbosity) {
    // re-seed MT so that `set.seed()` works in R
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // unpack control params
    double thresh = (double) control["adapt_k_thresh"];
    double alpha = (double) control["seq_alpha"];
    double pop_temper = (double) control["pop_temper"];
    double final_infl = (double) control["final_infl"];
    std::vector<int> lags = as<std::vector<int>>(control["lags"]);

    int cores = (int) control["cores"];
    if (cores <= 0) cores = std::thread::hardware_concurrency();
    if (cores == 1) cores = 0;

    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    if (districts.n_rows != V || districts.n_cols != N)
        throw std::range_error("Initialization districts have wrong dimensions.");
    double total_pop = sum(pop);
    bool check_both = total_pop/n_distr > lower && total_pop/n_distr < upper;
    double tol = std::max(target - lower, upper - target) / target;

    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "SEQUENTIAL MONTE CARLO\n";
        Rcout << "Sampling " << N << " " << V << "-unit ";
        if (n_drawn + n_steps + 1 == n_distr) {
            Rcout << "maps with " << n_distr << " districts and population between "
                << lower << " and " << upper << ".\n";
        } else {
            Rcout << "partial maps of " << n_drawn + n_steps + 1
                << " districts and population between "
                << lower << " and " << upper << ".\n";
        }
        if (cg.size() > 1)
            Rcout << "Ensuring no more than " << n_distr - 1 << " splits of the "
                  << cg.size() << " administrative units.\n";
        if (verbosity >= 3) {
            if (cores == 0) {
                Rcout << "Using 1 core.\n";
            } else {
                Rcout << "Using " << cores << " cores.\n";
            }

            if (check_both) {
                Rcout << "Using two-sided population checks.\n";
            } else {
                Rcout << "Using one-sided population checks.\n";
            }
        }
    }

    vec pop_left(N);
    if (n_drawn == 0) {
        pop_left.fill(total_pop);
    } else {
        // compute population not assigned (i.e., in district '0')
        pop_left.fill(0.0);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < V; j++) {
                if (districts(j, i) == 0) {
                    pop_left[i] += pop[j];
                }
            }
        }
    }

    vec log_temper(N, fill::zeros);
    vec lp(N, fill::zeros);
    umat ancestors(N, lags.size(), fill::zeros);

    std::vector<int> cut_k(n_steps);
    std::vector<int> n_unique(n_steps);
    std::vector<double> n_eff(n_steps);
    std::vector<double> accept_rate(n_steps);
    std::vector<double> sd_lp(n_steps);
    std::vector<double> sd_temper(n_steps);
    vec cum_wgt(N, fill::value(1.0 / N));
    cum_wgt = cumsum(cum_wgt);

    RcppThread::ThreadPool pool(cores);

    std::string bar_fmt = "Split [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} | ETA{cli::pb_eta}";
    RObject bar = cli_progress_bar(n_steps, cli_config(false, bar_fmt.c_str()));
    try {
    for (int ctr = n_drawn + 1; ctr <= n_drawn + n_steps; ctr++) {
        int i_split = ctr - n_drawn - 1;
        if (verbosity >= 3) {
            Rcout << "Making split " << ctr - n_drawn << " of " << n_steps;
        }

        // find k and multipliers
        int last_k = i_split == 0 ? std::max(1, V - 5) : cut_k[i_split - 1];
        adapt_parameters(g, cut_k[i_split], last_k, lp, thresh, tol, districts,
                         counties, cg, pop, pop_left, target, verbosity);

        if (verbosity >= 3) {
            Rcout << " (using k = " << cut_k[i_split] << ")\n";
        }

        // perform resampling/drawing
        bool final = ctr == n_drawn + n_steps;
        if (final) {
            lower = target - (target - lower) * final_infl;
            upper = target + (upper - target) * final_infl;
        }
        split_maps(g, counties, cg, pop, districts, cum_wgt, lp, pop_left,
                   log_temper, pop_temper, accept_rate[i_split],
                   n_distr, ctr, ancestors, lags, n_unique[i_split],
                   lower, upper, target,
                   rho, cut_k[i_split], check_both, pool, verbosity);

        sd_lp[i_split] = stddev(lp);
        sd_temper[i_split] = stddev(log_temper);

        // compute weights for next step
        cum_wgt = get_wgts(districts, n_distr, ctr, final, alpha, lp,
                           n_eff[i_split], pop, target, g, constraints,
                           verbosity);

        if (verbosity == 1 && CLI_SHOULD_TICK)
            cli_progress_set(bar, i_split);
        Rcpp::checkUserInterrupt();
    } // end for
    } catch (Rcpp::internal::InterruptedException e) {
        cli_progress_done(bar);
        return R_NilValue;
    }
    cli_progress_done(bar);

    lp = lp - log_temper;

    if (n_drawn + n_steps + 1 == n_distr) {
        // Set final district label to n_distr rather than 0
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < V; j++) {
                districts(j, i) = districts(j, i) == 0 ? n_distr : districts(j, i);
            }
        }
    }

    List out = List::create(
        _["plans"] = districts,
        _["lp"] = lp,
        _["ancestors"] = ancestors,
        _["sd_lp"] = sd_lp,
        _["sd_temper"] = sd_temper,
        _["est_k"] = cut_k,
        _["step_n_eff"] = n_eff,
        _["unique_survive"] = n_unique,
        _["accept_rate"] = accept_rate
    );

    return out;
}


/*
 * Helper function to iterate over constraints and apply them
 */
double add_constraint(const std::string& name, List constraints,
                      std::function<double(List)> fn_constr) {
    if (!constraints.containsElementNamed(name.c_str())) return 0;

    List constr = constraints[name];
    double val = 0;
    for (int i = 0; i < constr.size(); i++) {
        List constr_inst = constr[i];
        double strength = constr_inst["strength"];
        if (strength != 0) {
            val += strength * fn_constr(constr_inst);
        }
    }
    return val;
}

/*
 * Add specific constraint weights & return the cumulative weight vector
 */
vec get_wgts(const umat &districts, int n_distr, int distr_ctr, bool final,
             double alpha, vec &lp, double &neff,
             const uvec &pop, double parity, const Graph g,
             List constraints, int verbosity) {
    int V = districts.n_rows;
    int N = districts.n_cols;

    std::vector<int> distr_calc;
    if (final) {
        distr_calc = {distr_ctr, 0};
    } else {
        distr_calc = {distr_ctr};
    }

    if (constraints.size() > 0) {
    for (int i = 0; i < N; i++) {
        for (int j : distr_calc) {
            lp[i] += add_constraint("pop_dev", constraints,
                                      [&] (List l) -> double {
                                          return eval_pop_dev(districts.col(i), j,
                                                                 pop, parity);
                                      });

            lp[i] += add_constraint("status_quo", constraints,
                [&] (List l) -> double {
                    return eval_sq_entropy(districts.col(i), as<uvec>(l["current"]),
                                           j, pop, n_distr,
                                           as<int>(l["n_current"]), V);
                });

            lp[i] += add_constraint("segregation", constraints,
                                      [&] (List l) -> double {
                                          return eval_segregation(districts.col(i), j,
                                                                  as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                                      });

            lp[i] += add_constraint("grp_pow", constraints,
                [&] (List l) -> double {
                    return eval_grp_pow(districts.col(i), j,
                                        as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]),
                                        as<double>(l["tgt_group"]), as<double>(l["tgt_other"]),
                                        as<double>(l["pow"]));
                });

            lp[i] += add_constraint("compet", constraints,
                [&] (List l) -> double {
                    uvec dvote = l["dvote"];
                    uvec total = dvote + as<uvec>(l["rvote"]);
                    return eval_grp_pow(districts.col(i), j,
                                        dvote, total, 0.5, 0.5, as<double>(l["pow"]));
                });

            lp[i] += add_constraint("grp_hinge", constraints,
                [&] (List l) -> double {
                    return eval_grp_hinge(districts.col(i), j, as<vec>(l["tgts_group"]),
                                          as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                });

            lp[i] += add_constraint("grp_inv_hinge", constraints,
                                    [&] (List l) -> double {
                                        return eval_grp_hinge(districts.col(i), j, as<vec>(l["tgts_group"]),
                                                              as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                                    });

            lp[i] += add_constraint("incumbency", constraints,
                [&] (List l) -> double {
                    return eval_inc(districts.col(i), j, as<uvec>(l["incumbents"]));
                });

            lp[i] += add_constraint("splits", constraints,
                [&] (List l) -> double {
                    return eval_splits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
                });

            lp[i] += add_constraint("multisplits", constraints,
                [&] (List l) -> double {
                    return eval_multisplits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
                });

            lp[i] += add_constraint("total_splits", constraints,
                [&] (List l) -> double {
                    return eval_total_splits(districts.col(i), j, as<uvec>(l["admin"]), l["n"], true);
                });

            lp[i] += add_constraint("polsby", constraints,
                                      [&] (List l) -> double {
                                          return eval_polsby(districts.col(i), j,
                                                             as<ivec>(l["from"]),
                                                             as<ivec>(l["to"]), as<vec>(l["area"]),
                                                             as<vec>(l["perimeter"]));
                                      });

            lp[i] += add_constraint("fry_hold", constraints,
                                      [&] (List l) -> double {
                                          return eval_fry_hold(districts.col(i), j,
                                                               as<uvec>(l["total_pop"]),
                                                               as<mat>(l["ssdmat"]),
                                                               as<double>(l["denominator"]));
                                      });

            lp[i] += add_constraint("qps", constraints,
                                      [&] (List l) -> double {
                                          return eval_qps(districts.col(i), j,
                                                          as<uvec>(l["total_pop"]),
                                                          as<uvec>(l["cities"]), as<int>(l["n_city"]),
                                                          n_distr);
                                      });

            lp[i] += add_constraint("custom", constraints,
                [&] (List l) -> double {
                    Function fn = l["fn"];
                    return as<NumericVector>(fn(districts.col(i), j))[0];
                });
        }
    } // for
    } // if

    vec wgt = exp(-alpha * lp);
    if (!final) // not the last iteration
        lp = lp * (1 - alpha);
    vec cuml_wgt = cumsum(wgt);

    neff = cuml_wgt[N-1] * cuml_wgt[N-1]  / sum(square(wgt));
    if (verbosity >= 3) {
        Rcout << std::setprecision(1) << 100*neff/N <<  "% efficiency.\n";
    }

    return cuml_wgt / cuml_wgt[N-1];
}


/*
 * Split off a piece from each map in `districts`,
 * keeping deviation between `lower` and `upper`
 */
void split_maps(const Graph &g, const uvec &counties, Multigraph &cg,
                const uvec &pop, umat &districts, vec &cum_wgt, vec &lp,
                vec &pop_left, vec &log_temper, double pop_temper,
                double &accept_rate, int n_distr, int dist_ctr,
                umat &ancestors, const std::vector<int> &lags, int &n_unique,
                double lower, double upper, double target,
                double rho, int k, bool check_both,
                RcppThread::ThreadPool &pool, int verbosity) {
    const int V = districts.n_rows;
    const int N = districts.n_cols;
    const int new_size = n_distr - dist_ctr;
    const int n_cty = max(counties);
    const int n_lags = lags.size();

    umat districts_new(V, N);
    vec pop_left_new(N);
    vec lp_new(N);
    vec log_temper_new(N);
    umat ancestors_new(N, n_lags);
    uvec uniques(N);

    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50; // check for interrupts every _ iterations
    uvec iters(N, fill::zeros); // how many actual iterations

    RcppThread::ProgressBar bar(N, 1);
    pool.parallelFor(0, N, [&] (int i) {
        int reject_ct = 0;
        bool ok = false;
        int idx = i;
        double inc_lp;
        double lower_s = lower;
        double upper_s = upper;

        Tree ust = init_tree(V);
        std::vector<bool> visited(V);
        std::vector<bool> ignore(V);
        while (!ok) {
            // resample
            idx = r_int_wgt(N, cum_wgt);
            districts_new.col(i) = districts.col(idx);
            iters[i]++;

            if (check_both) {
                lower_s = std::max(lower, pop_left(idx) - new_size * upper);
                upper_s = std::min(upper, pop_left(idx) - new_size * lower);
            }

            if (lower_s >= upper_s) {
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);
                continue;
            }
            inc_lp = split_map(g, ust, counties, cg, districts_new.col(i),
                               dist_ctr, visited, ignore,
                               pop, pop_left(idx), lower_s, upper_s, target, k);

            // bad sample; try again
            if (!std::isfinite(inc_lp)) {
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);
                continue;
            }

            ok = true;
        }
        uniques[i] = idx;
        clear_tree(ust);

        // save ancestors/lags
        for (int j = 0; j < n_lags; j++) {
            if (dist_ctr <= lags[j]) {
                ancestors_new(i, j) = i;
            } else {
                ancestors_new(i, j) = ancestors(idx, j);
            }
        }

        // backwards kernel adjustment (for label counts)
        std::vector<bool> gr_bool(n_distr, false);
        for (int j = 0; j < V; j++) {
            if (districts_new.col(i)[j] != 0) continue;
            std::vector<int> nbors = g[j];
            for (int nbor : nbors) {
                int dist_k = districts_new.col(i)[nbor];
                if (dist_k != 0) {
                    gr_bool[dist_k] = true;
                }
            }
        }
        int nbors_zero = 0;
        for (int j = 0; j < n_distr; j++) {
            nbors_zero += gr_bool[j];
        }
        inc_lp += std::log(nbors_zero);


        // handle log_st compactness calc when needed
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
        double dev = std::fabs(lower_s - target)/target;
        double pop_pen = std::sqrt((double) n_distr - 2) * std::log(1e-12 + dev);
        log_temper_new(i) = log_temper(idx) + pop_temper*pop_pen;

        lp_new(i) = lp(idx) + inc_lp + pop_temper*pop_pen;

        if (verbosity >= 3) {
            bar++;
        }
        RcppThread::checkUserInterrupt(i % check_int == 0);
    });
    pool.wait();

    accept_rate = N / (1.0 * sum(iters));
    if (verbosity >= 3) {
        Rcout << "  " << std::setprecision(2) << 100.0 * accept_rate << "% acceptance rate, ";
    }

    districts = districts_new;
    pop_left = pop_left_new;
    lp = lp_new;
    log_temper = log_temper_new;
    ancestors = ancestors_new;
    n_unique = ((uvec) find_unique(uniques)).n_elem;
}

/*
 * Split a map into two pieces with population lying between `lower` and `upper`
 */
double split_map(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                 subview_col<uword> districts, int dist_ctr,
                 std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                 double total_pop, double &lower, double upper, double target, int k) {
    int V = g.size();

    for (int i = 0; i < V; i++) ignore[i] = districts(i) != 0;

    int root;
    clear_tree(ust);
    int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
    if (result != 0) return -std::log(0.0);

    double new_pop = cut_districts(ust, k, root, districts, dist_ctr, pop, total_pop,
                          lower, upper, target);

    if (new_pop == 0) {
        return -std::log(0.0); // reject sample
    } else {
        lower = new_pop;  // set `lower` as a way to return population of new district
        return log_boundary(g, districts, 0, dist_ctr);// - log((double) k); (k is constant)
    }
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
    candidates.reserve(k);
    deviances.reserve(k);
    is_ok.reserve(k);
    int distr_root = districts(root);
    for (int i = 1; i <= V; i++) { // 1-indexing here
        if (districts(i - 1) != distr_root || i - 1 == root) continue;
        double below = pop_below.at(i - 1);
        double dev1 = std::fabs(below - target);
        double dev2 = std::fabs(total_pop - below - target);
        if (dev1 < dev2) {
            candidates.push_back(i);
            deviances.push_back(dev1);
            is_ok.push_back(lower < below && below < upper);
        } else {
            candidates.push_back(-i);
            deviances.push_back(dev2);
            is_ok.push_back(lower < total_pop - below && total_pop - below < upper);
        }
    }
    if ((int) candidates.size() < k) return 0.0;

    int idx = r_int(k);
    idx = select_k(deviances, idx + 1);
    int cut_at = std::fabs(candidates[idx]) - 1;
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
void adapt_parameters(const Graph &g, int &k, int last_k, const vec &lp, double thresh,
                      double tol, const umat &districts, const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      const vec &pop_left, double target, int verbosity) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(10 + (int) (2.0 * V * tol), last_k + 4); // heuristic
    int N_max = districts.n_cols;
    int N_adapt = std::min(100, N_max);

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);

    std::vector<std::vector<double>> devs;
    devs.reserve(N_adapt);
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    int idx = 0;
    int max_V = 0;
    Tree ust = init_tree(V);
    for (int i = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (std::isinf(lp(i))) { // skip if not valid
            idx--;
            continue;
        }

        int n_vtx = V;
        for (int j = 0; j < V; j++) {
            if (districts(j, i) != 0) {
                ignore[j] = true;
                n_vtx--;
            }
        }
        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                    pop, lower, upper, counties, cg);
        if (result != 0) {
            idx--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, pop_left(i), target));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            if (devs.at(idx).at(j) <= tol) { // sorted
                n_ok++;
            } else {
                break;
            }
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max)
            max_ok = n_ok;

        Rcpp::checkUserInterrupt();
    }

    if (idx < N_adapt) N_adapt = idx; // if rejected too many in last step
    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    for (k = 1; k <= k_max; k++) {
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            double dev = devs.at(i).at(r_int(k));
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k >= k_max) {
        if (verbosity >= 3) {
            Rcout << " [maximum hit; falling back to naive k estimator]";
        }
        k = max_ok;
    }

    if (last_k < k_max && k < last_k * 0.6) k = (int) (0.5*k + 0.5*last_k);

    k = std::min(std::max(max_ok + 1, k) + 1 - (distr_ok(k) > 0.99) + (thresh == 1),
                 max_V - 1);
}
