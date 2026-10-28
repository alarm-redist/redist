/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2025/3
 * Purpose: Various functions called on plans 
 ********************************************************/



// #include <Rcpp.h>
#include <atomic>
#include <functional>
#include <sstream>
#include <thread>

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <unordered_set>
#include "advanced_types.h"
#include "random.h"
#include "scoring.h"
#include "base_plan_type.h"
#include "redist_alg_helpers.h"
#include "threading_helpers.h"


namespace {

/*
 * The canonical labelling of a plan written out as "1,1,2,3,...", with labels
 * numbered from 1 in order of first appearance.
 *
 * Only ever called once per distinct plan. Building one of these per particle
 * is exactly what `canonical_plan_hash` exists to avoid.
 */
template <typename PlanID>
std::string canonical_plan_string(PlanID const &region_ids, int const V,
                                  int const num_regions,
                                  std::vector<int> &label_scratch) {
    label_scratch.assign(num_regions, UNSEEN_REGION);

    std::string out;
    out.reserve(static_cast<std::size_t>(V) * 3);

    int next_label = 0;

    for (int v = 0; v < V; ++v) {
        int const region = static_cast<int>(region_ids[v]);

        int label = label_scratch[region];
        if (label == UNSEEN_REGION) {
            label = next_label++;
            label_scratch[region] = label;
        }

        if (v > 0) out.push_back(',');
        out += std::to_string(label + 1);
    }

    return out;
}

} // namespace



/*
 * Compute the cooccurence matrix for a set of precincts indexed by `idxs`,
 * given a collection of plans
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix prec_cooccur(Rcpp::IntegerMatrix m, Rcpp::IntegerVector idxs, int ncores = 1) {
    int v = m.nrow();
    int n = idxs.size();
    Rcpp::NumericMatrix out(v, v);

    RcppThread::parallelFor(
        0, v,
        [&](int i) {
            out(i, i) = 1;
            for (int j = 0; j < i; j++) {
                double shared = 0;
                for (int k = 0; k < n; k++) {
                    shared += m(i, idxs[k] - 1) == m(j, idxs[k] - 1);
                }
                shared /= n;
                out(i, j) = shared;
                out(j, i) = shared;
            }
        },
        resolve_num_threads(ncores));

    return out;
}



/*
 * Compute the percentage of `group` in each district. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix group_pct(Rcpp::IntegerMatrix const &plans_mat, Rcpp::NumericVector const &group_pop,
                        Rcpp::NumericVector const &total_pop, int const n_distr, int const ncores = 1) {
    int V = plans_mat.nrow();
    int num_plans = plans_mat.ncol();

    Rcpp::NumericMatrix grp_distr(n_distr, num_plans);
    Rcpp::NumericMatrix tot_distr(n_distr, num_plans);

    RcppThread::ThreadPool pool = get_thread_pool(ncores);

    pool.parallelFor(0, num_plans, [&](unsigned int i) {
        for (int j = 0; j < V; j++) {
            int distr = plans_mat(j, i) - 1;
            grp_distr(distr, i) += group_pop[j];
            tot_distr(distr, i) += total_pop[j];
        }
    });

    pool.wait();

    // divide
    pool.parallelFor(0, num_plans, [&](unsigned int i) {
        for (int j = 0; j < n_distr; j++) {
            grp_distr(j, i) /= tot_distr(j, i);
        }
    });

    pool.wait();

    return grp_distr;
}


/*
 * Compute the percentage of `group` in each district, and return the `k`-th
 * largest such value. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
Rcpp::NumericVector group_pct_top_k(const Rcpp::IntegerMatrix m, const Rcpp::NumericVector group_pop,
                              const Rcpp::NumericVector total_pop, int k, int n_distr) {
    int v = m.nrow();
    int n = m.ncol();
    Rcpp::NumericVector out(n);

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

        std::nth_element(grp_distr.begin(), grp_distr.begin() + k - 1, grp_distr.end(),
                         std::greater<double>());

        out[i] = grp_distr[k - 1];
    }

    return out;
}


// We assume that the population deviations are less than 1
// If they are 1 or greater then you cannot uniquely infer sizes
// [[Rcpp::export]]
Rcpp::IntegerMatrix infer_region_seats(Rcpp::IntegerMatrix const &region_pops,
                                       double const lower, double const upper,
                                       int const total_seats, int const num_threads = 1) {
    //
    int const num_plans = region_pops.ncol();
    int const num_regions = region_pops.nrow();

    bool bounds_issues = false;
    // warn if population bounds aren't tight
    for (int a_size = 2; a_size <= total_seats; a_size++) {
        if (upper * (a_size - 1) >= lower * a_size) {
            // REprintf("WARNING: Population bounds are not tight for size %d and %d\n",
            // a_size-1, a_size);
            // Rcpp::warning("Population bounds are not tight, inferring a unique number of
            // seats"
            //     " may not be possible.\n");
            bounds_issues = true;
        }
    }
    if (bounds_issues) {
        Rcpp::warning("Population bounds are not tight, inferring a unique number of seats"
                      " may not be possible.\n");
    }

    Rcpp::IntegerMatrix region_sizes(num_regions, num_plans);

    // parallel for loop over each plan
    RcppThread::parallelFor(
        0, num_plans,
        [&](unsigned int i) {
            // loop over each region
            for (int j = 0; j < num_regions; j++) {
                // get region pop
                auto region_pop = region_pops(j, i);
                int region_size;
                bool size_selected = false;
                // find the first instance in which the region is in bounds
                for (int potential_size = 1; potential_size <= total_seats; potential_size++) {
                    // see if this size works
                    if (lower * potential_size <= region_pop &&
                        region_pop <= upper * potential_size) {
                        region_size = potential_size;
                        size_selected = true;
                    }
                }
                if (!size_selected) {
                    std::ostringstream oss;
                    oss << "Could not infer the number of seats for region " << (j + 1)
                        << " of plan " << (i + 1) << ": its population " << region_pop
                        << " is not between " << lower << " * seats and " << upper
                        << " * seats for any number of seats from 1 to " << total_seats
                        << ". Pass the seats in explicitly or check the population bounds.\n";
                    throw std::runtime_error(oss.str());
                }

                region_sizes(j, i) = region_size;
            }
        },
        resolve_num_threads(num_threads));

    return region_sizes;
}


/*
 * Tally a variable by district.
 */
// TESTED
// NOTE: Maybe can make parallel version of this? Not sure
// [[Rcpp::export]]
Rcpp::NumericMatrix pop_tally(Rcpp::IntegerMatrix const &districts, Rcpp::NumericVector const &pop, int const n_distr,
                        int const ncores = 1) {
    int const num_plans = districts.ncol();
    int const V = districts.nrow();

    Rcpp::NumericMatrix tally(n_distr, num_plans);

    // parallel for loop over each plan
    RcppThread::parallelFor(
        0, num_plans,
        [&](unsigned int i) {
            for (int j = 0; j < V; j++) {
                int d = districts(j, i) - 1; // districts are 1-indexed
                tally(d, i) = tally(d, i) + pop(j);
            }
        },
        resolve_num_threads(ncores));

    return tally;
}

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// [[Rcpp::export]]
Rcpp::NumericVector max_dev(const Rcpp::IntegerMatrix &districts, const Rcpp::NumericVector &pop,
                            int const n_distr, bool const multimember_districts = false,
                            int const nseats = -1,
                            Rcpp::IntegerMatrix const &seats_matrix = Rcpp::IntegerMatrix(1, 1),
                            int const num_threads = 1) {
    int const num_plans = districts.ncol();

    Rcpp::NumericVector res(num_plans);
    Rcpp::NumericMatrix district_pops = pop_tally(districts, pop, n_distr, num_threads);

    if (multimember_districts) {
        double const target_pop = Rcpp::sum(pop) / nseats;
        RcppThread::parallelFor(
            0, num_plans,
            [&](unsigned int i) {
                for (int j = 0; j < n_distr; j++) {
                    double target_seat_pop = target_pop * seats_matrix(j, i);
                    double dev = std::fabs(district_pops(j, i) / target_seat_pop - 1.0);
                    // If deviation at district j bigger then record that
                    if (dev > res(i)) {
                        res(i) = dev;
                    }
                }
            },
            resolve_num_threads(num_threads));
    } else {
        double const target_pop = Rcpp::sum(pop) / n_distr;
        RcppThread::parallelFor(
            0, num_plans,
            [&](unsigned int i) {
                for (int j = 0; j < n_distr; j++) {
                    double dev = std::fabs(district_pops(j, i) / target_pop - 1.0);
                    // If deviation at district j bigger then record that
                    if (dev > res(i)) {
                        res(i) = dev;
                    }
                }
            },
            resolve_num_threads(num_threads));
    }

    return res;
}


// Given a numeric vector of statistics computed on each district this
// sorts the statistics within each plan.
// the length district_stats must be a multiple of ndists
// [[Rcpp::export]]
Rcpp::NumericVector order_district_stats(Rcpp::NumericVector const &district_stats,
                                         int const ndists, int const num_threads = 0) {

    if (district_stats.size() % ndists != 0) {
        throw Rcpp::exception("The length of the vector of district statistics must be a "
                              "multiple of the number of districts\n");
    } else if (ndists <= 1) {
        throw Rcpp::exception("Number of districts must be at least 2!\n");
    }

    int num_plans = district_stats.size() / ndists;

    Rcpp::NumericVector ordered_district_stats = Rcpp::clone(district_stats);

    RcppThread::parallelFor(
        0, num_plans,
        [&](unsigned int i) {
            // sort each chunk
            int start_index = i * ndists;
            // we don't subtract 1 since the end index is exclusive!!
            int end_index = start_index + ndists;

            std::sort(ordered_district_stats.begin() + start_index,
                      ordered_district_stats.begin() + end_index);
        },
        resolve_num_threads(num_threads));

    return ordered_district_stats;
}

// [[Rcpp::export]]
Rcpp::DataFrame order_columns_by_district(Rcpp::DataFrame const &df,
                                          Rcpp::CharacterVector const &columns,
                                          int const ndists, int const num_threads = 0) {

    Rcpp::List out(df.size()); // same number of columns
    Rcpp::CharacterVector names = df.names();

    /*
     * This loop is deliberately serial. Every statement in it touches the R
     * API: pulling a column out of the data frame, allocating the sorted
     * result, and assigning into the output list. R's allocator and garbage
     * collector are not thread safe, so running this across threads risks
     * heap corruption rather than a speedup.
     *
     * The parallelism belongs one level down instead. `order_district_stats`
     * sorts disjoint chunks of a plain numeric buffer and makes no R calls,
     * so it is safe to thread, and there are normally far more plans than
     * columns to spread across cores.
     */
    for (int i = 0; i < df.size(); ++i) {
        Rcpp::String colname = names[i];
        bool should_process = false;

        // check if this column is in the selected set
        for (int j = 0; j < columns.size(); ++j) {
            if (columns[j] == colname) {
                should_process = true;
                break;
            }
        }

        if (should_process) {
            Rcpp::NumericVector col = df[i];
            out[i] = order_district_stats(col, ndists, num_threads);
        } else {
            out[i] = df[i]; // leave unchanged
        }
    }

    out.attr("names") = names;
    out.attr("class") = df.attr("class");
    out.attr("row.names") = df.attr("row.names");

    return out;
}


RegionMultigraphCount build_region_multigraph(Graph const &g, PlanVector const &region_ids,
                                              int const num_regions) {
    RegionMultigraphCount region_multigraph(num_regions);
    int const V = g.size();

    for (int v = 0; v < V; v++) {
        // Find out which region this vertex corresponds to
        int v_region_num = region_ids[v];

        // now iterate over its neighbors
        for (int v_nbor : g[v]) {
            // find which region neighbor corresponds to
            int v_nbor_region_num = region_ids[v_nbor];

            // to avoid double counting only count when v less u
            if (v_region_num < v_nbor_region_num) {
                // we increase the count of edges
                region_multigraph[v_region_num][v_nbor_region_num]++;
                region_multigraph[v_nbor_region_num][v_region_num]++;
            }
        }
    }

    return region_multigraph;
}

Rcpp::NumericMatrix build_region_laplacian(RegionMultigraphCount const &region_multigraph) {
    int num_regions = region_multigraph.size();
    Rcpp::NumericMatrix laplacian_mat(num_regions, num_regions);
    // iterate over the multigraph
    for (size_t region_id = 0; region_id < num_regions; region_id++) {
        int vertex_degree = 0;
        // iterate over neighbors
        for (auto const &it : region_multigraph[region_id]) {
            // add number of edges to degree
            vertex_degree += it.second;
            laplacian_mat(region_id, it.first) = -it.second;
            laplacian_mat(it.first, region_id) = -it.second;
        }
        laplacian_mat(region_id, region_id) = vertex_degree;
    }

    return (laplacian_mat);
}

// Can call from R
// [[Rcpp::export]]
RegionMultigraphCount get_region_multigraph(Rcpp::List const &adj_list,
                                            Rcpp::IntegerVector const &region_ids) {
    std::unordered_set<int> uniqueElements;
    for (int element : region_ids) {
        uniqueElements.insert(element);
    }

    AllPlansVector underlying_id_vec(region_ids.begin(), region_ids.end());

    PlanVector region_id_vec(underlying_id_vec, 0, underlying_id_vec.size());

    int num_regions = uniqueElements.size();
    return (build_region_multigraph(list_to_graph(adj_list), region_id_vec, num_regions));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_region_laplacian(Rcpp::List const &adj_list,
                                         Rcpp::IntegerVector const &region_ids) {
    return (build_region_laplacian(get_region_multigraph(adj_list, region_ids)));
}


/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
// [[Rcpp::export]]
Rcpp::IntegerVector resample_lowvar(Rcpp::NumericVector wgts) {
    int N = wgts.size();

    double r = GLOBAL_RNG.r_unif() / N;
    double cuml = wgts[0];
    Rcpp::IntegerVector out(N);

    int i = 0;
    for (int n = 0; n < N; n++) {
        double u = r + n / (double)N;
        while (u > cuml) {
            cuml += wgts[++i]; // increment then access
        }
        out[n] = i + 1;
    }

    return out;
}


// [[Rcpp::export]]
double get_log_number_linking_edges(Rcpp::List const &adj_list, Rcpp::IntegerVector const &counties,
                                    Rcpp::List const &constraints, int const ndists,
                                    int const nseats, int const num_regions,
                                    Rcpp::IntegerVector const &region_ids) {
    int V = adj_list.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector)adj_list[i]));
    }

    MapParams const map_params(g, Rcpp::as<std::vector<unsigned int>>(counties), 
    {}, ndists, nseats, std::vector<int>{1}, 0,
                               0, 0, SamplingSpace::LinkingEdgeSpace);

    PlanMultigraph plan_multigraph(map_params, true);

    ScoringFunction scoring_function(map_params, constraints, 0, true);

    // need to make arma vector into plan multigraph
    std::vector<RegionID> flattened_all_plans(region_ids.begin(), region_ids.end());
    PlanVector plan_region_ids(flattened_all_plans, 0, plan_multigraph.map_params.V);

    plan_multigraph.build_plan_multigraph(plan_region_ids, num_regions);

    return plan_multigraph.compute_log_multigraph_tau(num_regions, scoring_function);
}

// [[Rcpp::export]]
double get_merged_log_number_linking_edges(Rcpp::List const &adj_list,
                                           Rcpp::IntegerVector const &counties,
                                           Rcpp::List const &constraints, int const ndists,
                                           int const nseats, int const num_regions,
                                           Rcpp::IntegerVector const &region_ids, int const region1_id,
                                           int const region2_id) {
    int V = adj_list.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector)adj_list[i]));
    }

    MapParams const map_params(g, 
        Rcpp::as<std::vector<unsigned int>>(counties), {}, ndists, nseats, std::vector<int>{1}, 0,
                               0, 0, SamplingSpace::LinkingEdgeSpace);

    PlanMultigraph plan_multigraph(map_params, true);

    ScoringFunction scoring_function(map_params, constraints, 0, true);

    // need to make arma vector into plan multigraph
    std::vector<RegionID> flattened_all_plans(region_ids.begin(), region_ids.end());
    PlanVector plan_region_ids(flattened_all_plans, 0, plan_multigraph.map_params.V);

    plan_multigraph.build_plan_multigraph(plan_region_ids, num_regions);

    return plan_multigraph.compute_merged_log_multigraph_tau(num_regions, region1_id,
                                                             region2_id, scoring_function);
}


// Check which plans are hierarchically valid
//
// Given a matrix of 1-indexed plans this checks, for each plan, whether it is
// hierarchically valid with respect to `counties`, meaning
//   - every region intersect county is a single connected piece;
//   - the number of county region components is small enough; and
//   - the administratively adjacent quotient graph has no cycles.
//
// This matters because the hierarchical Wilson sampler can never draw a tree
// on a region that splits a county into disconnected pieces, so an invalid
// plan makes the sampler retry until it gives up rather than fail fast.
//
// @param adj_list Zero indexed adjacency list of the map
// @param counties 1-indexed vector of county labels, one per vertex
// @param plans_mat A matrix of 1-indexed plans, one plan per column
// @param ndists The number of districts in each plan
//
// @returns A logical vector with one entry per column of `plans_mat`
//
// @keywords internal
// [[Rcpp::export]]
Rcpp::LogicalVector plans_are_hierarchically_valid(
    Rcpp::List const &adj_list,
    Rcpp::IntegerVector const &counties,
    Rcpp::IntegerMatrix const &plans_mat,
    int const ndists
) {
    int const V = adj_list.size();

    if (plans_mat.nrow() != V) {
        throw Rcpp::exception(
            "plans_are_hierarchically_valid: plans matrix must have one row "
            "per vertex in the adjacency list."
        );
    }
    if (counties.size() != V) {
        throw Rcpp::exception(
            "plans_are_hierarchically_valid: counties must have one entry "
            "per vertex in the adjacency list."
        );
    }
    if (ndists <= 0) {
        throw Rcpp::exception(
            "plans_are_hierarchically_valid: ndists must be positive."
        );
    }

    Graph g;
    g.reserve(V);
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector)adj_list[i]));
    }

    // Population and the population bounds play no part in this check, so
    // they are left empty. ndists is real because `PlanMultigraph` sizes its
    // per-region vectors by it.
    MapParams const map_params(
        g, Rcpp::as<std::vector<unsigned int>>(counties), {},
        ndists, ndists, std::vector<int>{1}, 0, 0, 0,
        SamplingSpace::GraphSpace
    );

    PlanMultigraph plan_multigraph(map_params);

    int const num_plans = plans_mat.ncol();
    Rcpp::LogicalVector out(num_plans);

    // Reused across plans; `is_hierarchically_valid` clears it itself.
    std::vector<bool> const component_lookup(
        static_cast<std::size_t>(ndists) * map_params.num_counties, false
    );

    std::vector<RegionID> region_id_buffer(V);

    for (int j = 0; j < num_plans; j++) {
        for (int v = 0; v < V; v++) {
            int const region = plans_mat(v, j);

            if (region < 1 || region > ndists) {
                std::ostringstream oss;
                oss << "plans_are_hierarchically_valid: plan " << (j + 1)
                    << " has region id " << region << " at vertex " << (v + 1)
                    << ", but region ids must be between 1 and " << ndists
                    << ".";

                throw Rcpp::exception(oss.str().c_str());
            }

            // internal region ids are 0-indexed
            region_id_buffer[v] = static_cast<RegionID>(region - 1);
        }

        PlanVector region_ids(region_id_buffer, 0, V);

        out[j] = plan_multigraph.is_hierarchically_valid(
            region_ids, ndists, component_lookup
        );
    }

    return out;
}


// Get canonically relabeled plans matrix
//
// Given a matrix of 1-indexed plans (or partial plans) this function
// returns a new plans matrix with all the plans labeled canonically.
// The canonical labelling of a plan is the one where the region of the
// first vertex gets mapped to 1, the region of the next smallest vertex
// in a different region than the first gets mapped to 2, and so on. This
// is guaranteed to result in the same labelling for any plan where the
// region ids have been permuted.
//
//
// @param plans_mat A matrix of 1-indexed plans
// @param num_regions The number of regions in the plan
// @param num_threads The number of threads to use. Defaults to number of machine threads.
//
// @details Modifications
//    - None
//
// @returns A matrix of canonically labelled plans
//
// @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerMatrix get_canonical_plan_labelling(Rcpp::IntegerMatrix const &plans_mat,
                                                 int const num_regions, int const ncores = 0) {
    int const V = plans_mat.nrow();
    int const nsims = plans_mat.ncol();
    // check the plan isn't zero indexed
    for (size_t i = 0; i < V; i++) {
        if (plans_mat(i, 0) == 0) {
            throw Rcpp::exception(
                "Plans matrix in `get_canonical_plan_labelling` must be 1-indexed!\n");
        }
    }

    Rcpp::IntegerMatrix relabelled_plan_mat(V, nsims);

    int const num_threads = ncores <= 0 ? std::thread::hardware_concurrency() : ncores;
    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);

    // trick to give each thread a unique id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};

    // make vectors which maps old region ids to the new canonical one
    std::vector<std::vector<int>> reindex_vecs(num_threads, std::vector<int>(num_regions));

    // now relabel
    pool.parallelFor(0, nsims, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }
        // reset the vector indices
        std::fill(reindex_vecs[thread_id].begin(), reindex_vecs[thread_id].end(), -1);

        int current_region_relabel_counter = 1;

        for (size_t v = 0; v < V; v++) {
            // check if this region has been relabelled yet
            if (reindex_vecs[thread_id][plans_mat(v, i) - 1] <= 0) {
                // if not then we haven't set a relabel for this region
                reindex_vecs[thread_id][plans_mat(v, i) - 1] = current_region_relabel_counter;
                ++current_region_relabel_counter;
            }

            // now relabel
            relabelled_plan_mat(v, i) = reindex_vecs[thread_id][plans_mat(v, i) - 1];
        }
    });

    pool.wait();

    return relabelled_plan_mat;
}



// Count how many times each plan appears in a plans matrix
//
// Given a matrix of 1-indexed plans (or partial plans) this function
// returns a list mapping plan vectors as a giant concatened string to
// the count of how many times the plan appears.
//
// If `use_canonical_ordering` is set to true then the plans will be
// reordered using the canonical reordering function
// `get_canonical_plan_labelling`. This guarantees that the same plan
// will not be incorrectly counted if there are different permutations
// of its labels. If `use_canonical_ordering` is not set to true then
// its possible the count will be incorrect because of different
// permutations of the same underlying plan.
//
//
// @param plans_mat A matrix of 1-indexed plans
// @param num_regions The number of regions in the plan
// @param use_canonical_ordering Whether or not to reorder the plans using the
// canonical ordering on plans.
// @param num_threads The number of threads to use. Defaults to number of machine threads.
//
// @details Modifications
//    - None
//
// @returns A list mapping plans (stored as a string concatened vector) to
// how many times they appear in the matrix
//
// @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame get_plan_counts(Rcpp::IntegerMatrix const &input_plans_mat,
                                int const num_regions,
                                bool const use_canonical_ordering = true,
                                bool const return_plan_strings = false) {
    int const V = input_plans_mat.nrow();
    int const nsims = input_plans_mat.ncol();

    // plans arrive 1-indexed, so region ids run 1..num_regions
    int const label_slots = num_regions + 1;

    auto const plan_at = [&](int const i) { return input_plans_mat.column(i); };

    std::vector<PlanTallyEntry> tally;

    if (use_canonical_ordering) {
        // canonicalising happens inside the hash, so no relabelled copy of
        // the plans matrix is ever built
        tally = tally_unique_plans(plan_at, nsims, V, label_slots);
    } else {
        // count region id permutations of the same partition separately
        std::unordered_map<std::uint64_t, int> heads;
        heads.reserve(nsims * 2);
        tally.reserve(nsims);

        for (int i = 0; i < nsims; ++i) {
            auto const plan = plan_at(i);

            std::uint64_t hash = 0xcbf29ce484222325ULL;
            for (int v = 0; v < V; ++v) {
                hash ^= static_cast<std::uint64_t>(plan[v]) + 1ULL;
                hash *= 0x100000001b3ULL;
            }
            hash = mix64_hash(hash);

            auto const inserted = heads.emplace(hash, static_cast<int>(tally.size()));
            if (inserted.second) {
                tally.push_back(PlanTallyEntry{i, 1, -1});
            } else {
                ++tally[inserted.first->second].count;
            }
        }
    }

    int const n_unique = static_cast<int>(tally.size());

    Rcpp::IntegerVector counts(n_unique);
    // 1-indexed column of the first plan carrying each distinct plan
    Rcpp::IntegerVector representative(n_unique);

    for (int i = 0; i < n_unique; ++i) {
        counts[i] = tally[i].count;
        representative[i] = tally[i].representative + 1;
    }

    if (!return_plan_strings) {
        return Rcpp::DataFrame::create(Rcpp::Named("count") = counts,
                                       Rcpp::Named("representative") = representative);
    }

    // Built once per distinct plan rather than once per particle.
    Rcpp::CharacterVector plan_string(n_unique);
    std::vector<int> label_scratch(label_slots);

    for (int i = 0; i < n_unique; ++i) {
        auto const plan = plan_at(tally[i].representative);
        plan_string[i] = canonical_plan_string(plan, V, label_slots, label_scratch);
    }

    return Rcpp::DataFrame::create(Rcpp::Named("plan_string") = plan_string,
                                   Rcpp::Named("count") = counts,
                                   Rcpp::Named("representative") = representative,
                                   Rcpp::Named("stringsAsFactors") = false);
}


// Checks a matrix of seat counts is valid
//
// Checks that a matrix of seat counts associated with a plan is valid
// meaning that every region has a positive seat value and for each plan
// the sum of seats is equal to the total number of seats (`nseats`).
// If anything is not correct an error will be thrown.
//
// @param init_seats A matrix of 1-indexed plans
// @param num_regions The number of regions in the plan.
// @param nseats The total number of seats in the map
// @param seats_range Vector of number of seats a district is allowed to have
// @param split_districts_only Whether or not to check that all but the last region are
// districts or not. (Allows for the possibility the last region is a district too).
// @param num_threads The number of threads to use. Defaults to number of machine threads.
//
// @details Modifications
//    - None
//
// @keywords internal
// @noRd
// [[Rcpp::export]]
void validate_init_seats_cpp(Rcpp::IntegerMatrix const &init_seats, int const num_regions,
                             int const nseats, Rcpp::IntegerVector const &seats_range,
                             bool const split_districts_only, int const num_threads) {
    // create thread pool
    RcppThread::ThreadPool pool = get_thread_pool(num_threads);

    // check matrix dimensions
    if (init_seats.nrow() != num_regions) {
        REprintf("Expected init_seats to have %d rows but actually had %d!\n", num_regions,
                 init_seats.nrow());
        throw Rcpp::exception("`init_seats` matrix did not have `num_regions` rows!\n");
    }

    int num_cols = init_seats.ncol();

    // get minimum district size
    int min_district_size = *std::min_element(seats_range.begin(), seats_range.end());
    std::vector<bool> is_district(nseats + 1, false);
    for (auto const a_size : seats_range) {
        is_district[a_size] = true;
    }

    // now check each column
    // `init_seats` is regions by plans, so a plan is a column: region j of plan
    // i is init_seats(j, i). Everything below reads through `seat` rather than
    // indexing again, since indexing it the other way round both checks the
    // wrong entry and runs off the end of the matrix whenever there are more
    // plans than regions.
    pool.parallelFor(0, num_cols, [&](int i) {
        // check each value is positive and sums to nseats
        int seat_sum = 0;
        for (int j = 0; j < num_regions; j++) {
            int const seat = init_seats(j, i);
            std::ostringstream oss;
            if (seat <= 0) {
                oss << "Region " << (j + 1) << " of plan " << (i + 1)
                    << " does not have a positive seat count (" << seat << ")! "
                    << "Non-positive seat values in `init_seats`!\n";
                throw Rcpp::exception(oss.str().c_str());
            } else if (seat < min_district_size) {
                oss << "Region " << (j + 1) << " of plan " << (i + 1) << " has " << seat
                    << " seats, smaller than the smallest district seat size ("
                    << min_district_size << ")!\n";
                throw Rcpp::exception(oss.str().c_str());
            } else if (seat > nseats) {
                oss << "Region " << (j + 1) << " of plan " << (i + 1) << " has " << seat
                    << " seats, more than `nseats` (" << nseats << ")!\n";
                throw Rcpp::exception(oss.str().c_str());
            }

            if (
                split_districts_only &&
                j + 1 != num_regions &&
                !is_district[seat]
            ) {
                oss << "Region " << (j + 1) << " of plan " << (i + 1) << " has " << seat
                    << " seats, which is not a district size, but only the final region may "
                       "be a non-district remainder!\n";
                throw Rcpp::exception(oss.str().c_str());
            }
            seat_sum += seat;
        }

        if (seat_sum != nseats) {
            std::ostringstream oss;
            oss << "The sum of seats in plan " << (i + 1) << " is " << seat_sum
                << ", which is not equal to `nseats` (" << nseats << ")!\n";
            throw Rcpp::exception(oss.str().c_str());
        }
    });

    pool.wait();

    return;
}

