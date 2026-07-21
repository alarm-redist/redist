#pragma once
#ifndef MAP_CALC_H
#define MAP_CALC_H


#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
#include "redist_types.h"




/*
 * Compute the Fryer-Holden penalty for district `distr`
 */
double eval_fry_hold(const arma::subview_col<arma::uword> &districts, int distr, const arma::uvec &total_pop,
                     arma::mat ssdmat, double denominator);

/*
 * Compute the qps penalty for district `distr`
 */
double eval_qps(const arma::subview_col<arma::uword> &districts, int distr, const arma::uvec &total_pop,
                const arma::uvec &cities, int n_city, int nd);

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_log_st(const arma::subview_col<arma::uword> &districts, const Graph g, arma::uvec counties,
                   int ndists);

/*
 * Compute the log spanning tree penalty for district `distr`
 */
double eval_er(const arma::subview_col<arma::uword> &districts, const Graph g, int ndists);






/*
 * Non-templated constraint functions
 */

double compute_log_pop_temper(double const target, double const pop_temper, int const ndists,
                              int const region_pop, int const region_size);

/************************
 * Templated Constraint Functions
 *************************/

/*
 * Compute the population penalty for district `distr`
 */
template <typename PlanID>
double eval_pop_dev(const PlanID &region_ids, int const region1, int const region2,
                    std::vector<unsigned int> const &total_pop, double const parity) {
    double pop = 0.0;

    for (size_t i = 0; i < total_pop.size(); ++i) {
        if (region_ids[i] == region1 || region_ids[i] == region2) {
            pop += total_pop[i];
        }
    }

    double frac = pop / parity;
    return std::pow(frac - 1.0, 2.0);
}

/*
 * Compute the power-based group penalty for district `distr`
 */
template <typename PlanID>
double eval_grp_pow(const PlanID &region_ids, int const V, int const region1_id,
                    int const region2_id, arma::uvec const &grp_pop,
                    arma::uvec const &total_pop, double const tgt_grp, double const tgt_other,
                    double const pow) {
    double sum_grp = 0.0;
    double sum_total = 0.0;

    for (size_t i = 0; i < V; ++i) {
        if (region_ids[i] == region1_id || region_ids[i] == region2_id) {
            sum_grp += grp_pop[i];
            sum_total += total_pop[i];
        }
    }

    if (sum_total == 0.0)
        return 0.0; // avoid div-by-zero

    double frac = sum_grp / sum_total;
    return std::pow(std::fabs(frac - tgt_grp) * std::fabs(frac - tgt_other), pow);
}

/*
 * Compute the new, hinge group penalty for district `distr`
 *
 */
template <typename PlanID>
double eval_grp_hinge(PlanID const &region_ids, int const V, int const region1_id,
                      int const region2_id, arma::vec const &tgts_grp,
                      const arma::uvec &grp_pop, const arma::uvec &total_pop) {
    double subsetted_grp_pop_sum = 0.0;
    double subsetted_total_pop_sum = 0.0;
    // get the sum of the two columns in region 1 or 2
    for (size_t i = 0; i < V; i++) {
        auto const region_id = region_ids[i];
        if (region_id == region1_id || region_id == region2_id) {
            subsetted_grp_pop_sum += grp_pop[i];
            subsetted_total_pop_sum += total_pop[i];
        }
    }
    // do subsetted_grp_pop_sum/subsetted_total_pop_sum
    double frac = std::exp(std::log(subsetted_grp_pop_sum) - std::log(subsetted_total_pop_sum));
    // REprintf("Frac %f\n", frac);
    // figure out which to compare it to
    double target;
    double diff = 1;
    int n_tgt = tgts_grp.size();
    for (int i = 0; i < n_tgt; i++) {
        double new_diff = std::fabs(tgts_grp[i] - frac);
        // REprintf("Target %d - Diff %f, Target %f \n", i, new_diff, tgts_grp[i]);
        if (new_diff <= diff) {
            diff = new_diff;
            target = tgts_grp[i];
        }
    }
    // REprintf("%f\n", target - frac);

    return std::sqrt(std::max(0.0, target - frac));
}

/*
 * Compute the incumbent-preserving penality
 */
template <typename PlanID>
double eval_inc(PlanID const &region_ids, int const region1_id, int const region2_id,
                const arma::uvec &incumbents) {
    int n_inc = incumbents.size();
    double inc_in_distr = -1.0; // first incumbent doesn't count
    for (int i = 0; i < n_inc; i++) {
        if (region_ids[incumbents[i] - 1] == region1_id ||
            region_ids[incumbents[i] - 1] == region2_id)
            inc_in_distr++;
    }

    if (inc_in_distr < 0.0) {
        inc_in_distr = 0.0;
    }

    return inc_in_distr;
};

/*
 * Compute the status quo penalty for district `distr`, maybe needs to be revised ...
 */
template <typename PlanID>
double eval_sq_entropy(PlanID const &region_ids, std::vector<unsigned int> const &current,
                       int const region1_id, int const region2_id, std::vector<unsigned int> const &pop,
                       int const ndists, int const n_current, int const V) {
    double accuml = 0;
    for (int j = 0; j < n_current; j++) { // 0-indexed districts, last one is remainder
        double pop_overlap = 0;
        double pop_total = 0;
        for (int k = 0; k < V; k++) {
            if (current[k] != j)
                continue;
            pop_total += pop[k];

            if (region_ids[k] == region1_id || region_ids[k] == region2_id)
                pop_overlap += pop[k];
        }
        double frac = pop_overlap / pop_total;
        if (frac > 0)
            accuml += frac * std::log(frac);
    }

    return -accuml / ndists / std::log(n_current);
}

// helper function
// calculates districts which appear in each county (but not zeros)
template <typename PlanID>
std::vector<std::set<int>> calc_county_dist(PlanID const &region_ids,
                                            arma::uvec const &counties, int const n_cty,
                                            bool const zero_ok) {
    std::vector<std::set<int>> county_dist(n_cty);
    int V = counties.size();
    for (int i = 0; i < n_cty; i++) {
        county_dist[i] = std::set<int>();
    }
    for (int i = 0; i < V; i++) {
        if (zero_ok || region_ids[i] > 0) {
            county_dist[counties[i] - 1].insert(region_ids[i]);
        }
    }
    return county_dist;
}

/*
 * Compute the county split penalty for region region1
 * No merged version because I'm not sure what exactly its doing so instead
 * the region constraint wrapper just duplicates the plan id and makes one
 * where the regions are actually merged
 */
template <typename PlanID>
double eval_splits(PlanID const &region_ids, int const region_id, arma::uvec const &admin_units,
                   int const n_admin_units, bool const smc) {
    std::vector<std::set<int>> county_dist =
        calc_county_dist(region_ids, admin_units, n_admin_units, region_id == 0);

    int splits = 0;
    for (int i = 0; i < n_admin_units; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 2 : cty_n_distr >= 2;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(region_id);
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
 * Compute the county multisplit penalty for region `region_id`
 * No merged version because I'm not sure what exactly its doing so instead
 * the region constraint wrapper just duplicates the plan id and makes one
 * where the regions are actually merged
 */
template <typename PlanID>
double eval_multisplits(PlanID const &region_ids, int const region_id,
                        const arma::uvec &admin_units, int const n_admin_units,
                        bool const smc) {
    std::vector<std::set<int>> county_dist =
        calc_county_dist(region_ids, admin_units, n_admin_units, region_id == 0);

    double splits = 0;
    for (int i = 0; i < n_admin_units; i++) {
        int cty_n_distr = county_dist[i].size();
        // for SMC, just count the split when it crosses the threshold
        // for MCMC there is no sequential nature
        bool cond = smc ? cty_n_distr == 3 : cty_n_distr >= 3;
        if (cond) {
            if (smc) {
                auto search = county_dist[i].find(region_id);
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
 * Compute the total splits penalty for region `region_id`
 * No merged version because I'm not sure what exactly its doing so instead
 * the region constraint wrapper just duplicates the plan id and makes one
 * where the regions are actually merged
 */
template <typename PlanID>
double eval_total_splits(PlanID const &region_ids, int const region_id,
                         arma::uvec const &admin_units, int const n_admin_units,
                         bool const smc) {
    std::vector<std::set<int>> county_dist =
        calc_county_dist(region_ids, admin_units, n_admin_units, region_id == 0);

    double splits = 0;
    for (int i = 0; i < n_admin_units; i++) {
        int cty_n_distr = county_dist[i].size();
        // no over-counting since every split counts
        if (cty_n_distr > 1) {
            if (smc) {
                auto search = county_dist[i].find(region_id);
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
 * Compute the Polsby Popper penalty for region `region1`
 * or the region formed by merging region1 and region2 if region1 != region2
 */
template <typename PlanID>
double eval_polsby(PlanID const &region_ids, int const region1_id, int const region2_id,
                   int const V, arma::ivec const &from, arma::ivec const &to,
                   arma::vec const &area, arma::vec const &perimeter) {
    double tot_area = 0.0;
    double tot_perim = 0.0;
    constexpr double pi4 = 4.0 * 3.14159265;

    // Sum area directly for vertices in the district
    for (int v = 0; v < V; ++v) {
        if (region_ids[v] == region1_id || region_ids[v] == region2_id) {
            tot_area += area[v];
        }
    }

    // Sum perimeter contributions from boundary edges
    int E = to.n_elem;
    for (int e = 0; e < E; ++e) {
        auto out_vertex = to[e] - 1;
        if (region_ids[out_vertex] == region1_id || region_ids[out_vertex] == region2_id) {
            int src_vertex = from[e] - 1;
            // only count if its on the permiter
            // need to subtract one because R is 1 indexed
            if (src_vertex <= -1 || (region_ids[src_vertex] != region1_id &&
                                     region_ids[src_vertex] != region2_id)) {
                tot_perim += perimeter[e];
            }
        }
    }

    double dist_peri2 = std::pow(tot_perim, 2.0);

    return 1.0 - (pi4 * tot_area / dist_peri2);
}

/*
 * Compute the segregation penalty for district `distr`
 */
template <typename PlanID>
double eval_segregation(const PlanID &region_ids, int const region1_id, int const region2_id,
                        int const V, const arma::uvec &grp_pop, const arma::uvec &total_pop) {
    // Step 1: compute overall group share (pAll) and total population
    double total_grp = arma::sum(grp_pop);
    double total_pop_sum = arma::sum(total_pop);

    if (total_pop_sum == 0.0)
        return 0.0;

    double pAll = total_grp / total_pop_sum;
    double denom = 2.0 * total_pop_sum * pAll * (1.0 - pAll);
    if (denom == 0.0)
        return 0.0;

    // Step 2: sum population and group population in the target district
    double local_grp = 0.0;
    double local_pop = 0.0;

    for (int i = 0; i < V; ++i) {
        // count if in either region
        if (region_ids[i] == region1_id || region_ids[i] == region2_id) {
            local_grp += grp_pop[i];
            local_pop += total_pop[i];
        }
    }

    if (local_pop == 0.0)
        return 0.0;

    // Step 3: compute and return segregation score
    return local_pop * std::abs((local_grp / local_pop) - pAll) / denom;
}

#endif
