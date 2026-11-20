#pragma once
#ifndef MAP_CALC_H
#define MAP_CALC_H


#include <cstddef>
#include <cstdlib>

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <vector>
#include "redist_types.h"
#include "sparse_logdet.h"




// NOTE: eval_fry_hold is now templated and lives with the other templated
// constraint functions at the bottom of this header.

// NOTE: eval_log_st is now templated and lives with the other templated
// constraint functions at the bottom of this header.

// NOTE: eval_er is now templated and lives with the other templated constraint
// functions at the bottom of this header.






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
                    int const region2_id, std::vector<unsigned int> const &grp_pop,
                    std::vector<unsigned int> const &total_pop, double const tgt_grp, double const tgt_other,
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
                      int const region2_id, std::vector<double> const &tgts_grp,
                      const std::vector<unsigned int> &grp_pop, const std::vector<unsigned int> &total_pop) {
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
                const std::vector<unsigned int> &incumbents) {
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
                                            std::vector<unsigned int> const &counties, int const n_cty,
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
double eval_splits(PlanID const &region_ids, int const region_id, std::vector<unsigned int> const &admin_units,
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
                        const std::vector<unsigned int> &admin_units, int const n_admin_units,
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
                         std::vector<unsigned int> const &admin_units, int const n_admin_units,
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
                   int const V, std::vector<int> const &from, std::vector<int> const &to,
                   std::vector<double> const &area, std::vector<double> const &perimeter) {
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
    int E = static_cast<int>(to.size());
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
 * Compute the log number of spanning trees of the subgraph induced by the
 * vertices that are both in region `region_id` and in county `county`.
 *
 * This is Kirchhoff's theorem: the number of spanning trees equals any cofactor
 * of the Laplacian, so we build the Laplacian with its first row and column
 * dropped and take its log determinant. Mirrors
 * `compute_log_region_and_county_spanning_tree` in base_plan_type.cpp.
 */
// `check_connectivity` guards against region-and-county pieces that are not
// connected. Such a piece has no spanning trees at all, so the count is 0 and
// this returns -infinity. It has to be caught rather than left to the
// factorization: the Laplacian minor is positive definite only when the
// subgraph is connected, so compute_log_det_from_entries would throw.
//
// Callers that already guarantee connected pieces should pass false and skip
// the check. Tree-based samplers are in that position: they cut a
// county-constrained spanning tree, in which every region-and-county piece is
// already a connected subtree, and removing one edge leaves each side
// connected. Flip-style samplers move one precinct at a time and check only
// whole-district contiguity, so they can fragment a district's presence inside
// a county and must pass true.
template <typename PlanID>
double eval_log_region_county_st(PlanID const &region_ids, Graph const &g,
                                 std::vector<unsigned int> const &counties,
                                 int const region_id, int const county,
                                 bool const check_connectivity = true) {
    int const V = static_cast<int>(g.size());

    // number of vertices in this region-and-county intersection
    int K = 0;
    std::vector<int> pos(V); // position of each vertex within the subgraph
    int start = 0;           // where to start the second loop, to save time
    for (int i = 0; i < V; i++) {
        pos[i] = K - 1; // minus one because we drop the 1st row and column
        if (region_ids[i] == region_id && static_cast<int>(counties[i]) == county) {
            K++;
            if (K == 2) start = i; // 2nd vertex
        }
    }
    // a single vertex has exactly one spanning tree, and log(1) = 0
    if (K <= 1) return 0.0;

    if (check_connectivity) {
        std::vector<bool> visited(V, false);
        std::vector<int> to_visit;
        for (int i = 0; i < V; i++) {
            if (region_ids[i] == region_id && static_cast<int>(counties[i]) == county) {
                to_visit.push_back(i);
                visited[i] = true;
                break;
            }
        }

        int reached = 1;
        while (!to_visit.empty()) {
            int const u = to_visit.back();
            to_visit.pop_back();
            for (auto const nbor : g[u]) {
                if (visited[nbor] || region_ids[nbor] != region_id ||
                    static_cast<int>(counties[nbor]) != county) {
                    continue;
                }
                visited[nbor] = true;
                reached++;
                to_visit.push_back(nbor);
            }
        }

        // disconnected: no spanning trees, so log(0) = -infinity
        if (reached < K) return -std::numeric_limits<double>::infinity();
    }

    int const minor_V = K - 1;
    std::vector<SparseEntry> entries;
    // degrees are accumulated first and appended at the end, since editing
    // entries after the fact is expensive
    std::vector<int> vertex_degrees(minor_V, 0);

    for (int i = start; i < V; i++) {
        if (region_ids[i] != region_id || static_cast<int>(counties[i]) != county) continue;

        int const prec = pos[i];
        if (prec < 0) continue; // the dropped first row

        for (auto const nbor : g[i]) {
            if (region_ids[nbor] != region_id || static_cast<int>(counties[nbor]) != county) {
                continue;
            }
            vertex_degrees[prec]++;
            if (pos[nbor] < 0) continue; // the dropped first column
            // only the upper triangle, ie row <= column
            if (prec <= pos[nbor]) {
                entries.emplace_back(prec, pos[nbor], -1.0);
            }
        }
    }

    for (int v = 0; v < minor_V; v++) {
        entries.emplace_back(v, v, static_cast<double>(vertex_degrees[v]));
    }

    return compute_log_det_from_entries(entries, minor_V);
}

/*
 * Compute the log number of spanning trees of the county-level multigraph for
 * region `region_id`, ie the graph whose vertices are the counties the region
 * touches and whose edges are the region's edges that cross county lines.
 *
 * Mirrors `compute_log_county_level_spanning_tree` in base_plan_type.cpp.
 */
template <typename PlanID>
double eval_log_county_level_st(PlanID const &region_ids, Graph const &g,
                                std::vector<unsigned int> const &counties, int const n_cty,
                                int const region_id) {
    // with a single county there is nothing to contract, and log(1) = 0
    if (n_cty == 1) return 0.0;

    int const V = static_cast<int>(g.size());

    int K = 0;
    std::vector<int> pos(V, -2);      // position of each vertex's county in the subgraph
    std::vector<int> seen(n_cty, -2); // county lookup
    int start = 0;
    for (int i = 0; i < V; i++) {
        if (region_ids[i] != region_id) continue;

        if (seen[counties[i] - 1] < 0) {
            pos[i] = K - 1; // minus one because we drop the 1st row and column
            seen[counties[i] - 1] = K;
            K++;
            if (K == 2) start = i; // 2nd vertex
        } else {
            pos[i] = seen[counties[i] - 1] - 1;
        }
    }
    if (K <= 1) return 0.0;

    int const minor_V = K - 1;
    std::vector<SparseEntry> entries;
    std::vector<int> vertex_degrees(minor_V, 0);

    for (int i = start; i < V; i++) {
        if (region_ids[i] != region_id) continue;

        int const cty = pos[i];
        if (cty < 0) continue; // the dropped first row

        for (auto const nbor : g[i]) {
            // skip if outside the region, or inside the same county
            if (region_ids[nbor] != region_id || pos[nbor] == cty) continue;

            vertex_degrees[cty]++;
            if (pos[nbor] < 0) continue; // the dropped first column
            if (cty <= pos[nbor]) {
                entries.emplace_back(cty, pos[nbor], -1.0);
            }
        }
    }

    for (int v = 0; v < minor_V; v++) {
        entries.emplace_back(v, v, static_cast<double>(vertex_degrees[v]));
    }

    return compute_log_det_from_entries(entries, minor_V);
}

/*
 * Compute the log number of spanning trees of a whole plan.
 *
 * For each region this is the sum over the counties it touches of the spanning
 * trees of the region-and-county intersection, plus the spanning trees of the
 * region's county-level multigraph. That matches the convention used by
 * `Plan::compute_log_region_spanning_tree` in base_plan_type.cpp.
 *
 * Regions are 1-indexed, as are counties.
 *
 * `check_connectivity` is passed through to eval_log_region_county_st; see the
 * note there for when it can be skipped.
 */
template <typename PlanID>
double eval_log_st(PlanID const &region_ids, Graph const &g,
                   std::vector<unsigned int> const &counties, int const ndists,
                   bool const check_connectivity = true) {
    int const n_cty = counties.empty()
                          ? 0
                          : static_cast<int>(*std::max_element(counties.begin(), counties.end()));

    double log_st = 0.0;
    for (int region_id = 1; region_id <= ndists; region_id++) {
        for (int county = 1; county <= n_cty; county++) {
            log_st += eval_log_region_county_st(region_ids, g, counties, region_id, county,
                                                check_connectivity);
        }
        log_st += eval_log_county_level_st(region_ids, g, counties, n_cty, region_id);
    }

    return log_st;
}

/*
 * Compute the Fryer-Holden penalty for region `region_id`
 *
 * `ssdmat` is the V x V matrix of squared distances between vertices. It is
 * taken as an Rcpp matrix rather than an arma one so that it is referenced in
 * place: converting to an arma::mat copies the whole V x V matrix on every
 * call, and this is called once per region per iteration.
 */
template <typename PlanID>
double eval_fry_hold(PlanID const &region_ids, int const region_id, int const V,
                     std::vector<unsigned int> const &total_pop,
                     Rcpp::NumericMatrix const &ssdmat, double const denominator = 1.0) {
    // collect the vertices belonging to this region
    std::vector<int> idxs;
    for (int i = 0; i < V; i++) {
        if (region_ids[i] == region_id) idxs.push_back(i);
    }

    double ssd = 0.0;
    // `i + 1 < size` rather than `i < size - 1` so an empty region does not
    // underflow the unsigned bound and walk off the end
    for (std::size_t i = 0; i + 1 < idxs.size(); i++) {
        for (std::size_t k = i + 1; k < idxs.size(); k++) {
            ssd += ssdmat(idxs[i], idxs[k]) *
                   static_cast<double>(total_pop[idxs[i]]) *
                   static_cast<double>(total_pop[idxs[k]]);
        }
    }

    return ssd / denominator;
}

/*
 * Maps every region to itself.
 *
 * Lets the plain and merged edges removed scores share one implementation:
 * the plain caller passes this and the compiler folds the lookup away, while
 * the merged caller passes a real reindex vector.
 */
struct IdentityReindex {
    template <typename T>
    constexpr T operator[](T const region_id) const noexcept { return region_id; }
};


/*
 * Compute the edges removed penalty, ie the number of edges in the graph whose
 * two endpoints fall in different regions.
 *
 * `reindex` maps each region id to the region it should be counted as. Passing
 * `IdentityReindex` scores the plan as it stands; passing a vector that sends
 * one region's id to another scores the plan as though those two regions had
 * been merged, without having to build the merged plan.
 *
 * Mirrors `redistmetrics::n_removed` for a single plan: the loop sees every
 * crossing edge once from each of its endpoints, so the total is halved. The
 * number of districts is not needed, which is why `redistmetrics::n_removed`
 * never uses the `n_distr` it takes.
 */
template <typename PlanID, typename Reindex>
double eval_er(PlanID const &region_ids, Reindex const &reindex, Graph const &g) {
    int const V = static_cast<int>(g.size());

    double removed = 0.0;
    for (int i = 0; i < V; i++) {
        // hoisted out of the neighbor loop, so the reindex costs one lookup
        // per edge rather than two
        auto const region_id = reindex[region_ids[i]];
        for (int const nbor : g[i]) {
            if (reindex[region_ids[nbor]] != region_id) removed += 1.0;
        }
    }

    return removed / 2.0;
}


/*
 * Edges removed for the plan exactly as it is.
 */
template <typename PlanID>
double eval_er(PlanID const &region_ids, Graph const &g) {
    return eval_er(region_ids, IdentityReindex{}, g);
}


/*
 * Number of edges running between two regions.
 *
 * Merging two regions can only ever un-cut edges that ran between them, so
 * this is exactly the drop in the edges removed score from merging them. Used
 * to check the merged score against the plain one.
 */
template <typename PlanID>
double count_edges_between_regions(PlanID const &region_ids, int const region1_id,
                                   int const region2_id, Graph const &g) {
    int const V = static_cast<int>(g.size());

    double between = 0.0;
    for (int i = 0; i < V; i++) {
        auto const region_id = region_ids[i];
        if (region_id != region1_id && region_id != region2_id) continue;

        for (int const nbor : g[i]) {
            auto const nbor_region = region_ids[nbor];
            if (nbor_region != region_id &&
                (nbor_region == region1_id || nbor_region == region2_id)) {
                between += 1.0;
            }
        }
    }

    return between / 2.0;
}

/*
 * Compute the segregation penalty for district `distr`
 */
template <typename PlanID>
double eval_segregation(const PlanID &region_ids, int const region1_id, int const region2_id,
                        int const V, const std::vector<unsigned int> &grp_pop, const std::vector<unsigned int> &total_pop) {
    // Step 1: compute overall group share (pAll) and total population
    double total_grp = std::accumulate(grp_pop.begin(), grp_pop.end(), 0.0);
    double total_pop_sum = std::accumulate(total_pop.begin(), total_pop.end(), 0.0);

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
