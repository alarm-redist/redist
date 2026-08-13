/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2025/3
 * Purpose: Various functions called on plans 
 ********************************************************/



// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>
#include <unordered_set>
#include "advanced_types.h"
#include "random.h"
#include "scoring.h"
#include "base_plan_type.h"
#include "threading_helpers.h"



/*
 * Compute the cooccurence matrix for a set of precincts indexed by `idxs`,
 * given a collection of plans
 */
// [[Rcpp::export]]
arma::mat prec_cooccur(arma::umat m, arma::uvec idxs, int ncores) {
    int v = m.n_rows;
    int n = idxs.n_elem;
    arma::mat out(v, v);

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
        ncores);

    return out;
}



/*
 * Compute the percentage of `group` in each district. Asummes `m` is 1-indexed.
 */
// [[Rcpp::export]]
Rcpp::NumericMatrix group_pct(Rcpp::IntegerMatrix const &plans_mat, arma::vec const &group_pop,
                        arma::vec const &total_pop, int const n_distr, int const ncores = 0) {
    int V = plans_mat.nrow();
    int num_plans = plans_mat.ncol();

    Rcpp::NumericMatrix grp_distr(n_distr, num_plans);
    Rcpp::NumericMatrix tot_distr(n_distr, num_plans);

    // 0 or 1 ncores means no threading
    RcppThread::ThreadPool pool(ncores > 1 ? ncores : 0);

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
                                       int const total_seats, int const num_threads = 0) {
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
                    REprintf("No valid size could be found for Plan %i\n", i + 1);
                    throw Rcpp::exception("No valid size could be inferred!\n");
                }

                region_sizes(j, i) = region_size;
            }
        },
        num_threads > 0 ? num_threads : 0);

    return region_sizes;
}


/*
 * Tally a variable by district.
 */
// TESTED
// NOTE: Maybe can make parallel version of this? Not sure
// [[Rcpp::export]]
Rcpp::NumericMatrix pop_tally(Rcpp::IntegerMatrix const &districts, arma::vec const &pop, int const n_distr,
                        int const ncores = 0) {
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
        ncores > 1 ? ncores : 0);

    return tally;
}

/*
 * Compute the maximum deviation from the equal population constraint.
 */
// [[Rcpp::export]]
Rcpp::NumericVector max_dev(const Rcpp::IntegerMatrix &districts, const arma::vec &pop,
                            int const n_distr, bool const multimember_districts = false,
                            int const nseats = -1,
                            Rcpp::IntegerMatrix const &seats_matrix = Rcpp::IntegerMatrix(1, 1),
                            int const num_threads = 1) {
    int const num_plans = districts.ncol();

    Rcpp::NumericVector res(num_plans);
    Rcpp::NumericMatrix district_pops = pop_tally(districts, pop, n_distr, num_threads);

    if (multimember_districts) {
        double const target_pop = arma::sum(pop) / nseats;
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
            num_threads > 0 ? num_threads : 0);
    } else {
        double const target_pop = arma::sum(pop) / n_distr;
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
            num_threads > 0 ? num_threads : 0);
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
        num_threads > 0 ? num_threads : 0);

    return ordered_district_stats;
}

// [[Rcpp::export]]
Rcpp::DataFrame order_columns_by_district(Rcpp::DataFrame const &df,
                                          Rcpp::CharacterVector const &columns,
                                          int const ndists, int const num_threads = 0) {

    Rcpp::List out(df.size()); // same number of columns
    Rcpp::CharacterVector names = df.names();

    RcppThread::parallelFor(
        0, df.size(),
        [&](unsigned int i) {
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
                out[i] = order_district_stats(col, ndists, 1);
            } else {
                out[i] = df[i]; // leave unchanged
            }
        },
        num_threads > 0 ? num_threads : 0);

    out.attr("names") = names;
    out.attr("class") = df.attr("class");
    out.attr("row.names") = df.attr("row.names");

    return out;
}

namespace{

/*
 * Convert zero-indxed R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const Rcpp::List &l) {
    int V = l.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector)l[i]));
    }
    return g;
}

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

arma::mat build_region_laplacian(RegionMultigraphCount const &region_multigraph) {
    int num_regions = region_multigraph.size();
    arma::mat laplacian_mat(num_regions, num_regions, arma::fill::zeros);
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
                                            arma::uvec const &region_ids) {
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
arma::mat get_region_laplacian(Rcpp::List const &adj_list, arma::uvec const &region_ids) {
    return (build_region_laplacian(get_region_multigraph(adj_list, region_ids)));
}



// [[Rcpp::export]]
Rcpp::IntegerVector resample_plans_lowvar(Rcpp::NumericVector const &normalized_weights,
                                          Rcpp::IntegerMatrix &plans_mat,
                                          Rcpp::IntegerMatrix &region_pops_mat,
                                          Rcpp::IntegerMatrix &region_sizes_mat,
                                          bool const reorder_sizes_mat) {
    // generate resampling index
    int const nsims = normalized_weights.size();

    int rng_seed = (int)Rcpp::sample(INT_MAX, 1)[0];
    RNGState rng_state(rng_seed, 42);
    double r = rng_state.r_unif() / nsims;
    double cuml = normalized_weights[0];
    Rcpp::IntegerVector resample_index(nsims);
    std::vector<bool> index_unchanged(nsims);

    int i = 0;
    for (int n = 0; n < nsims; n++) {
        double u = r + n / (double)nsims;
        while (u > cuml) {
            cuml += normalized_weights[++i]; // increment then access
        }
        // resample_index maps entry i to its new value
        // `resample_index[i] = k` means you should replace plan i with plan k
        resample_index[n] = i;
    }

    std::vector<int> buffer(nsims);
    // Now we're going to reorder things one row at a time
    // makes algout$plans[i] now equal to algout$plans[rs_idx[i]]
    // so we're mapping i -> rs_idx[i]

    int const V = plans_mat.nrow();
    for (int row = 0; row < V; ++row) {
        // copy current row
        for (int col = 0; col < nsims; ++col) {
            buffer[col] = plans_mat(row, col);
        }
        // write back in new order
        for (int col = 0; col < nsims; ++col) {
            plans_mat(row, col) = buffer[resample_index[col]];
        }
    }

    // reorder the region populations
    int const num_regions = region_pops_mat.nrow();
    for (int row = 0; row < num_regions; ++row) {
        // copy current row
        for (int col = 0; col < nsims; ++col) {
            buffer[col] = region_pops_mat(row, col);
        }
        // write back in new order
        for (int col = 0; col < nsims; ++col) {
            region_pops_mat(row, col) = buffer[resample_index[col]];
        }
    }

    // Reorder region sizes if needed
    if (reorder_sizes_mat) {
        for (int row = 0; row < num_regions; ++row) {
            // copy current row
            for (int col = 0; col < nsims; ++col) {
                buffer[col] = region_sizes_mat(row, col);
            }
            // write back in new order
            for (int col = 0; col < nsims; ++col) {
                region_sizes_mat(row, col) = buffer[resample_index[col]];
            }
        }
    }

    // make resampling thing 1 indexed
    std::transform(resample_index.begin(), resample_index.end(), resample_index.begin(),
                   [](int x) { return x + 1; });

    return resample_index;
}


/*
 * Generate an integer vector of resampling indices with a low-variance resampler.
 */
// [[Rcpp::export]]
arma::ivec resample_lowvar(arma::vec wgts) {
    int N = wgts.n_elem;

    double r = GLOBAL_RNG.r_unif() / N;
    double cuml = wgts[0];
    arma::ivec out(N);

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
                                    arma::uvec const &region_ids) {
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
                                           arma::uvec const &region_ids, int const region1_id,
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
                                                 int const num_regions, int const ncores) {
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
    RcppThread::ThreadPool pool(num_threads > 1 ? num_threads : 0);

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
                                int const num_regions, bool const use_canonical_ordering = true,
                                int const num_threads = 1) {

    Rcpp::IntegerMatrix plans_mat =
        use_canonical_ordering
            ? get_canonical_plan_labelling(input_plans_mat, num_regions, num_threads)
            : input_plans_mat;

    RcppThread::ThreadPool pool = get_thread_pool(num_threads);
    int const V = plans_mat.nrow();
    int const nsims = plans_mat.ncol();

    std::vector<std::unordered_map<std::string, int>> plan_count_maps_vec(
        pool.getNumThreads() == 0 ? 1 : pool.getNumThreads());

    // trick to give each thread a unique id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};

    pool.parallelFor(0, nsims, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }

        std::ostringstream oss;

        for (int row = 0; row < V; ++row) {
            oss << plans_mat(row, i);
            if (row < V - 1) {
                oss << ",";
            }
        }
        auto key = oss.str();
        plan_count_maps_vec[thread_id][key]++;
    });

    pool.wait();

    // now combine into one map
    std::unordered_map<std::string, int> pattern_counts;
    for (size_t i = 0; i < plan_count_maps_vec.size(); i++) {
        for (const auto &pair : plan_count_maps_vec[i]) {
            pattern_counts[pair.first] += pair.second;
        }
    }

    // Convert to R named vector
    // Fill vectors for DataFrame
    std::vector<std::string> patterns;
    std::vector<int> counts;

    for (const auto &pair : pattern_counts) {
        patterns.push_back(pair.first);
        counts.push_back(pair.second);
    }

    return Rcpp::DataFrame::create(Rcpp::Named("plan_string") = patterns,
                                   Rcpp::Named("count") = counts);
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
    pool.parallelFor(0, num_cols, [&](int i) {
        // check each value is positive and sums to nseats
        int seat_sum = 0;
        for (size_t j = 0; j < num_regions; j++) {
            int const seat = init_seats(static_cast<int>(j), i);
            if (init_seats(i, j) <= 0) {
                REprintf("Region %zu of plan %i does not have a positive seat count (%d)!\n",
                         j + 1, i + 1, init_seats(i, j));
                throw Rcpp::exception("Non-positive seat values in `init_seats`!\n");
            } else if (init_seats(i, j) < min_district_size) {
                REprintf("Region %zu of plan %i has a seat size smaller than the smallest "
                         "district seat size!\n",
                         j + 1, i + 1);
                throw Rcpp::exception(
                    "Seat values in `init_seats` smaller than smallest district seat size!\n");
            } else if (init_seats(i, j) > nseats) {
                REprintf(
                    "Region %zu of plan %i has %d seats, more than `nseats` (%d) number of "
                    "seats!\n",
                    j + 1, i + 1, init_seats(i, j), nseats);
                throw Rcpp::exception("Seat values greater than `nseats` in `init_seats`!\n");
            }

            if (
                split_districts_only &&
                j + 1 != num_regions &&
                !is_district[seat]
            ) {
                throw Rcpp::exception(
                    "Non-remainder region is not a district!\n"
                );
            }
            seat_sum += init_seats(i, j);
        }

        if (seat_sum != nseats) {
            REprintf(
                "The sum of seats in plan %i is %d, which is not equal to `nseats` values %d\n",
                i + 1, seat_sum, nseats);
            throw Rcpp::exception(
                "Sum of seat values in a plan is not equal to `nseats` in `init_seats`!\n");
        }
    });

    pool.wait();

    return;
}

