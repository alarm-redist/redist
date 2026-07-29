/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2024/10
 * Purpose: Helper functions for all redist algorithm types
 ********************************************************/

#include "redist_alg_helpers.h"

#include "forest_plan_type.h"
#include "graph_plan_type.h"
#include "lct_graph_plan_type.h"
#include "linking_edge_plan_type.h"
#include "splitting_schedule_types.h"
#include "tree_splitting.h"
#include "threading_helpers.h"
#include "wilson.h"
#include "random.h"
#include "utils.h"
#include <RcppThread.h>

Graph list_to_graph(const Rcpp::List &l) {
    int V = l.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(Rcpp::as<std::vector<int>>((Rcpp::IntegerVector)l[i]));
    }
    return g;
}

// [[Rcpp::export]]
Rcpp::List maximum_input_sizes() {
    // Return results
    Rcpp::List out = Rcpp::List::create(Rcpp::_["max_V"] = MAX_SUPPORTED_NUM_VERTICES,
                            Rcpp::_["max_districts"] = MAX_SUPPORTED_NUM_DISTRICTS,
                            Rcpp::_["max_counties"] = MAX_SUPPORTED_NUM_COUNTIES);

    return out;
}



void set_merged_region_reindex_vec(int const num_regions, std::vector<int> &region_reindex_vec,
                                   int const region1_id, int const region2_id) {
    for (int region_id = 0; region_id < num_regions; region_id++) {
        if (region_id == region2_id) {
            region_reindex_vec[region2_id] = region1_id;
        } else {
            region_reindex_vec[region_id] = region_id;
        }
    }

    return;
}



// creates plan ensemble of blank plans
PlanEnsemble::PlanEnsemble(MapParams const &map_params, int const total_pop, int const nsims,
                           SamplingSpace const sampling_space, RcppThread::ThreadPool &pool,
                           int const verbosity)
    : nsims(nsims), V(map_params.V), ndists(map_params.ndists),
      total_seats(map_params.total_seats), sampling_space(sampling_space),
      flattened_all_plans(V * nsims, 0),
      flattened_all_region_sizes(ndists * nsims, 0),
      flattened_all_region_pops(ndists * nsims, 0),
      flattened_all_region_order_added(ndists * nsims, -1), 
      num_forest_edge_bit_words_per_plan(
      sampling_space == SamplingSpace::ForestSpace ||
      sampling_space == SamplingSpace::LinkingEdgeSpace
          ? map_params.num_edge_bit_words
          : 0
    ),
    flattened_all_forest_edge_bits(
        nsims * num_forest_edge_bit_words_per_plan,
        EdgeBitWord{0}
    ),
      plan_ptr_vec(nsims) {
    if (ndists < 2)
        throw Rcpp::exception("Tried to create a plan with fewer than 2 districts!");

    bool const use_graph_space = sampling_space == SamplingSpace::GraphSpace;
    bool const use_forest_space = sampling_space == SamplingSpace::ForestSpace;
    bool const use_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    bool const use_lct_graph_space = sampling_space == SamplingSpace::LCTGraphSpace;
    // create the plans
    if (verbosity >= 3) {
        Rcpp::Rcout << "Creating Blank Plans!" << std::endl;
    }

    RcppThread::ProgressBar bar(nsims, 1);
    pool.parallelFor(0, nsims, [&](int i) {
        // create the plan attributes for this specific plan
        PlanVector plan_region_ids(flattened_all_plans, V * i, V * (i + 1));
        RegionSizes plan_sizes(flattened_all_region_sizes, ndists * i, ndists * (i + 1));
        IntPlanAttribute plan_pops(flattened_all_region_pops, ndists * i, ndists * (i + 1));
        IntPlanAttribute plan_region_order_added(flattened_all_region_order_added, ndists * i,
                                                 ndists * (i + 1));
        PlanEdgeBits plan_forest_edge_bits(
            flattened_all_forest_edge_bits,
            i * num_forest_edge_bit_words_per_plan,
            (i + 1) * num_forest_edge_bit_words_per_plan
        );

        // create the plans
        if (use_graph_space) {
            plan_ptr_vec[i] =
                std::make_unique<GraphPlan>(total_seats, total_pop, plan_region_ids, plan_sizes,
                                            plan_pops, plan_region_order_added, plan_forest_edge_bits);
        } else if (use_lct_graph_space) {
            plan_ptr_vec[i] =
                std::make_unique<LCTGraphPlan>(total_seats, total_pop, plan_region_ids,
                                               plan_sizes, plan_pops, plan_region_order_added, plan_forest_edge_bits);
        } else if (use_forest_space) {
            plan_ptr_vec[i] =
                std::make_unique<ForestPlan>(total_seats, total_pop, plan_region_ids,
                                             plan_sizes, plan_pops, plan_region_order_added, plan_forest_edge_bits);
        } else if (use_linking_edge_space) {
            plan_ptr_vec[i] = std::make_unique<LinkingEdgePlan>(
                total_seats, total_pop, plan_region_ids, plan_sizes, plan_pops,
                plan_region_order_added, plan_forest_edge_bits);
        } else {
            throw Rcpp::exception("Input is invalid\n");
        }
        if (verbosity >= 3) {
            ++bar;
        }
    });

    pool.wait();

}

// creates plan ensemble of partial plans
PlanEnsemble::PlanEnsemble(MapParams const &map_params,
                           SplittingSchedule const &splitting_schedule, int const num_regions,
                           int const nsims, SamplingSpace const sampling_space,
                           Rcpp::IntegerMatrix const &plans_mat,
                           Rcpp::IntegerMatrix const &region_sizes_mat,
                           std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool,
                           int const verbosity)
    : nsims(nsims), V(map_params.V), ndists(map_params.ndists),
      total_seats(map_params.total_seats), sampling_space(sampling_space),
      flattened_all_plans(plans_mat.begin(), plans_mat.end()),
      flattened_all_region_sizes(ndists * nsims, 0),
      flattened_all_region_pops(ndists * nsims, 0),
      flattened_all_region_order_added(ndists * nsims, -1), 
      num_forest_edge_bit_words_per_plan(
      sampling_space == SamplingSpace::ForestSpace ||
      sampling_space == SamplingSpace::LinkingEdgeSpace
          ? map_params.num_edge_bit_words
          : 0
    ),
    flattened_all_forest_edge_bits(
        nsims * num_forest_edge_bit_words_per_plan,
        EdgeBitWord{0}
    ),
      plan_ptr_vec(nsims) {
    // make sure 0 indexed plans were not passed in
    if (*std::min_element(flattened_all_plans.begin(), flattened_all_plans.end()) <= 0) {
        throw Rcpp::exception(
            "The initial plans passed in are zero-indexed. They should only be 1 indexed!\n");
    }

    // Now subtract 1 from plans
    std::transform(flattened_all_plans.begin(), flattened_all_plans.end(),
                   flattened_all_plans.begin(), [](int x) { return x - 1; });

    // check matrix dimensions
    if (plans_mat.ncol() != nsims) {
        REprintf("The number of columns (%d) in the initial plan matrix was not equal to nsims "
                 "(%d)!\n",
                 plans_mat.ncol(), nsims);
        throw Rcpp::exception(
            "The number of columns in the initial plan matrix was not equal to nsims!\n");
    }
    if (region_sizes_mat.ncol() != nsims) {
        REprintf("The number of columns (%d) in the initial sizes matrix was not equal to "
                 "nsims (%d)!\n",
                 region_sizes_mat.ncol(), nsims);
        throw Rcpp::exception(
            "The number of columns in the initial sizes matrix was not equal to nsims!\n");
    }
    if (plans_mat.nrow() != V) {
        REprintf(
            "The number of rows (%d) in the initial plan matrix , was not equal to V (%d)!\n",
            plans_mat.nrow(), V);
        throw Rcpp::exception(
            "The number of rows in the initial plan matrix , was not equal to V!\n");
    }
    if (region_sizes_mat.nrow() != num_regions) {
        REprintf("The number of rows (%d) in the initial sizes matrix, was not equal to "
                 "initial number of regions (%d)!\n",
                 region_sizes_mat.nrow(), num_regions);
        throw Rcpp::exception(
            "The number of rows in the initial sizes matrix was not equal to ndists!\n");
    }

    // check num_regions and num_districts inputs make sense
    if (ndists < 2)
        throw Rcpp::exception("Tried to create a plan with ndists < 2 regions!");
    if (num_regions > ndists)
        throw Rcpp::exception("Tried to create a plan object with more regions than ndists!");
    if (num_regions == 0)
        throw Rcpp::exception("Tried to create a plan with 0 regions");
    // Now move the data in the matrix

    bool const use_graph_space = sampling_space == SamplingSpace::GraphSpace;
    bool const use_forest_space = sampling_space == SamplingSpace::ForestSpace;
    bool const use_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    bool const use_lct_graph_space = sampling_space == SamplingSpace::LCTGraphSpace;

    if (!use_graph_space && !use_lct_graph_space) {
        if (rng_states.size() > (pool.getNumThreads() == 0 ? 1 : pool.getNumThreads())) {
            throw Rcpp::exception("RNG States vector is more than the number of threads!\n");
        }
    }

    if (verbosity >= 3) {
        Rcpp::Rcout << "Loading Partial Plans!" << std::endl;
    }

    // trick to give each thread a unique id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};

    auto tree_ptr = std::make_unique<UniformValidSplitter>(map_params.map_graph);

    int const n_threads = get_num_threads(pool); 
    std::vector<USTSampler> ust_sampler_buffers(n_threads,
                                                USTSampler(map_params, splitting_schedule));
    std::vector<PlanMultigraph> plan_multigraph_buffers(n_threads, PlanMultigraph(map_params));
    std::vector<Graph> region_graph_buffers(n_threads, Graph(num_regions));

    RcppThread::ProgressBar bar(nsims, 1);
    pool.parallelFor(0, nsims, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }


        USTSampler &ust_sampler = ust_sampler_buffers[thread_id];
        PlanMultigraph &plan_multigraph = plan_multigraph_buffers[thread_id];
        Graph &region_graph = region_graph_buffers[thread_id];
        // create the plan attributes for this specific plan
        PlanVector plan_region_ids(flattened_all_plans, V * i, V * (i + 1));
        RegionSizes plan_sizes(flattened_all_region_sizes, ndists * i, ndists * (i + 1));
        IntPlanAttribute plan_pops(flattened_all_region_pops, ndists * i, ndists * (i + 1));
        IntPlanAttribute plan_region_order_added(flattened_all_region_order_added, ndists * i,
                                                 ndists * (i + 1));
        PlanEdgeBits plan_forest_edge_bits(
            flattened_all_forest_edge_bits,
            i * num_forest_edge_bit_words_per_plan,
            (i + 1) * num_forest_edge_bit_words_per_plan
        );

        // copy the sizes from matrix into the vector
        std::copy(region_sizes_mat.begin() + i * num_regions,
                  region_sizes_mat.begin() + (i + 1) * num_regions, plan_sizes.begin());

        // create the plans
        if (use_graph_space) {
            plan_ptr_vec[i] =
                std::make_unique<GraphPlan>(num_regions, map_params.pop, plan_region_ids,
                                            plan_sizes, plan_pops, plan_region_order_added, plan_forest_edge_bits);
        } else if (use_lct_graph_space) {
            plan_ptr_vec[i] =
                std::make_unique<LCTGraphPlan>(num_regions, map_params.pop, plan_region_ids,
                                               plan_sizes, plan_pops, plan_region_order_added, plan_forest_edge_bits);
        } else if (use_forest_space) {
            plan_ptr_vec[i] = std::make_unique<ForestPlan>(
                ndists, num_regions, map_params.pop, plan_region_ids, plan_sizes, plan_pops,
                plan_region_order_added, plan_forest_edge_bits, 
                ust_sampler_buffers[thread_id],
                rng_states[thread_id]);
        } else if (use_linking_edge_space) {
            plan_ptr_vec[i] = std::make_unique<LinkingEdgePlan>(
                ndists, num_regions, map_params.pop, plan_region_ids, plan_sizes, plan_pops,
                plan_region_order_added, plan_forest_edge_bits, 
                *tree_ptr, ust_sampler, plan_multigraph, region_graph,
                rng_states[thread_id]);
        } else {
            throw Rcpp::exception("This plan type not supported!\n");
        }
        if (verbosity >= 3) {
            ++bar;
        }
    });

    pool.wait();
}

Rcpp::IntegerMatrix PlanEnsemble::get_R_plans_matrix() {
    // make the plans matrix
    Rcpp::IntegerMatrix plan_mat(V, nsims);
    // copy data over
    std::copy(flattened_all_plans.begin(), flattened_all_plans.end(), plan_mat.begin());
    // now add 1 to everything
    std::transform(plan_mat.begin(), plan_mat.end(), plan_mat.begin(),
                   [](int x) { return x + 1; });
    return plan_mat;
}

Rcpp::IntegerMatrix PlanEnsemble::get_R_sizes_matrix(RcppThread::ThreadPool &pool) {
    int const num_regions = plan_ptr_vec[0]->num_regions;
    // make the sizes matrix
    Rcpp::IntegerMatrix sizes_mat(num_regions, nsims);
    // to avoid wasting space if not ndists we don't copy all
    if (num_regions < ndists) {
        // copy over the non-zero sizes for each plan
        pool.parallelFor(0, nsims, [&](int i) {
            std::copy(plan_ptr_vec[i]->region_sizes.begin(),
                      plan_ptr_vec[i]->region_sizes.begin() + num_regions,
                      sizes_mat.column(i).begin() // Start of column in Rcpp::IntegerMatrix
            );
        });
        pool.wait();
    } else {
        // else we can just copy the entire vector
        std::copy(flattened_all_region_sizes.begin(), flattened_all_region_sizes.end(),
                  sizes_mat.begin());
    }

    return sizes_mat;
}

int PlanEnsemble::count_unique_plans(RcppThread::ThreadPool &pool) const {
    int const num_regions = plan_ptr_vec[0]->num_regions;
    int const num_threads = get_num_threads(pool);
    std::vector<std::unordered_set<std::string>> plan_count_maps_vec(num_threads);

    // Trick to give each thread a unique id
    // We need the extra steps to avoid the problem where
    // only some threads persist from previous calls but the counter resets
    // resulting in multiple threads with the same id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};

    std::vector<std::vector<RegionID>> reindexed_plan_vecs(num_threads,
                                                           std::vector<RegionID>(V));
    std::vector<std::vector<int>> reindex_vecs(num_threads, std::vector<int>(num_regions));

    pool.parallelFor(0, nsims, [&](int i) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }

        if (thread_id < 0 || thread_id >= num_threads) {
            RcppThread::Rcerr << "Thread id out of range: " << thread_id << " / " << num_threads
                              << "\n";
            throw std::runtime_error("Thread id out of range");
        }

        // reset the vector indices
        std::fill(reindex_vecs[thread_id].begin(), reindex_vecs[thread_id].end(), -1);

        int current_region_relabel_counter = 0;

        for (size_t v = 0; v < V; v++) {
            auto v_region = plan_ptr_vec[i]->region_ids[v];

            // check if this region has been relabelled yet
            if (reindex_vecs[thread_id][v_region] < 0) {
                // if not then we haven't set a relabel for this region
                reindex_vecs[thread_id][v_region] = current_region_relabel_counter;
                ++current_region_relabel_counter;
            }

            if (reindex_vecs[thread_id][v_region] < 0 ||
                reindex_vecs[thread_id][v_region] >= num_regions) {
                // REprintf("%d Regions - Vertex %d - Reindexed Region %d to %d \n",
                //     num_regions, v, v_region, reindex_vec[v_region]
                // );
                // REprintf("Reindex Vec is: ");
                // for (int j = 0; j < reindex_vec.size(); j++)
                // {
                //     REprintf("|%d -> %d| ", j, reindex_vecs[thread_id][j]);
                // }
                // REprintf("\nPlan Vectors:");
                // for (int u = 0; u < V; u++)
                // {
                //      REprintf("|%d -> %d| ",
                //         plan_ptr_vec[i]->region_ids[u],
                //         reindex_vecs[thread_id][plan_ptr_vec[i]->region_ids[u]]
                //     );
                // }
                // REprintf("\n");
                // throw Rcpp::exception("!\n");
                std::ostringstream oss;

                oss << num_regions << " Regions - Vertex " << v << " - Reindexed Region "
                    << v_region << " to " << reindex_vecs[thread_id][v_region] << '\n';

                oss << "Reindex Vec is: ";
                for (int j = 0; j < static_cast<int>(reindex_vecs[thread_id].size()); ++j) {
                    oss << '|' << j << " -> " << reindex_vecs[thread_id][j] << "| ";
                }

                oss << "\nPlan Vectors:";
                for (int u = 0; u < V; ++u) {
                    int rid = plan_ptr_vec[i]->region_ids[u];
                    int mapped =
                        (rid >= 0 && rid < static_cast<int>(reindex_vecs[thread_id].size()))
                            ? reindex_vecs[thread_id][rid]
                            : -1; // sentinel if out of range
                    oss << '|' << rid << " -> " << mapped << "| ";
                }
                oss << '\n';

                // Safe in workers: buffers and flushes on main thread at pool.wait()/join()
                Rcpp::Rcout << oss.str();

                // Throw a plain C++ exception from the worker (safe). Catch on main thread and
                // convert to R error if desired.
                throw std::range_error("Out\n");
            }

            // now relabel
            reindexed_plan_vecs[thread_id][v] = reindex_vecs[thread_id][v_region];
        }

        // now hash the plan
        std::ostringstream oss;

        for (int row = 0; row < V; ++row) {
            oss << reindexed_plan_vecs[thread_id][row];
            if (row < V - 1) {
                oss << ",";
            }
        }
        auto key = oss.str();
        plan_count_maps_vec[thread_id].insert(key);
    });
    pool.wait();

    // now combine into one map
    std::unordered_set<std::string> pattern_counts;
    for (size_t i = 0; i < plan_count_maps_vec.size(); i++) {
        for (const auto &plans_str : plan_count_maps_vec[i]) {
            pattern_counts.insert(plans_str);
        }
    }

    return pattern_counts.size();
}

Rcpp::IntegerMatrix PlanEnsemble::get_region_pops_matrix(RcppThread::ThreadPool &pool) {
    int const num_regions = plan_ptr_vec[0]->num_regions;
    // make the sizes matrix
    Rcpp::IntegerMatrix pops_mat(num_regions, nsims);
    // to avoid wasting space if not ndists we don't copy all
    if (num_regions < ndists) {
        // copy over the non-zero sizes for each plan
        pool.parallelFor(0, nsims, [&](int i) {
            std::copy(plan_ptr_vec[i]->region_pops.begin(),
                      plan_ptr_vec[i]->region_pops.begin() + num_regions,
                      pops_mat.column(i).begin() // Start of column in Rcpp::IntegerMatrix
            );
        });
        pool.wait();
    } else {
        // else we can just copy the entire vector
        std::copy(flattened_all_region_pops.begin(), flattened_all_region_pops.end(),
                  pops_mat.begin());
    }

    return pops_mat;
}

PlanEnsemble get_plan_ensemble(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity) {
    if (num_regions == 1) {
        return PlanEnsemble(map_params, 
            std::accumulate(map_params.pop.begin(), map_params.pop.end(), 0),
            nsims, sampling_space, pool,
                            verbosity);
    } else {
        return PlanEnsemble(map_params, splitting_schedule, num_regions, nsims, sampling_space,
                            plans_mat, region_sizes_mat, rng_states, pool, verbosity);
    }
}

std::unique_ptr<PlanEnsemble> get_plan_ensemble_ptr(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity) {
    if (num_regions == 1) {
        return std::make_unique<PlanEnsemble>(map_params, 
            std::accumulate(map_params.pop.begin(), map_params.pop.end(), 0), 
            nsims, sampling_space, pool, verbosity);
    } else {
        return std::make_unique<PlanEnsemble>(map_params, splitting_schedule, num_regions,
                                              nsims, sampling_space, plans_mat,
                                              region_sizes_mat, rng_states, pool, verbosity);
    }
}


void PlanEnsemble::check_all_plans_valid(
    MapParams const &map_params,
    std::string_view where
) {
    bool const check_forest = sampling_space != SamplingSpace::GraphSpace;
    std::ostringstream oss;

    auto fail = [&](std::string const &msg) {
        std::ostringstream full;
        full << "PlanEnsemble validity check failed.\n";
        full << "Where: " << where << "\n";
        full << msg << "\n";
        throw std::runtime_error(full.str());
    };

    if (nsims < 0) {
        fail("nsims is negative.");
    }

    if (V != map_params.V) {
        std::ostringstream msg;
        msg << "PlanEnsemble V does not match map_params.V. "
            << "V=" << V << ", map_params.V=" << map_params.V;
        fail(msg.str());
    }

    if (ndists != map_params.ndists) {
        std::ostringstream msg;
        msg << "PlanEnsemble ndists does not match map_params.ndists. "
            << "ndists=" << ndists << ", map_params.ndists=" << map_params.ndists;
        fail(msg.str());
    }

    if (static_cast<int>(plan_ptr_vec.size()) != nsims) {
        std::ostringstream msg;
        msg << "plan_ptr_vec.size() != nsims. "
            << "plan_ptr_vec.size()=" << plan_ptr_vec.size()
            << ", nsims=" << nsims;
        fail(msg.str());
    }

    if (static_cast<int>(flattened_all_plans.size()) != nsims * V) {
        std::ostringstream msg;
        msg << "flattened_all_plans has wrong size. "
            << "size=" << flattened_all_plans.size()
            << ", expected=" << nsims * V;
        fail(msg.str());
    }

    if (static_cast<int>(flattened_all_region_sizes.size()) != nsims * ndists) {
        std::ostringstream msg;
        msg << "flattened_all_region_sizes has wrong size. "
            << "size=" << flattened_all_region_sizes.size()
            << ", expected=" << nsims * ndists;
        fail(msg.str());
    }

    if (static_cast<int>(flattened_all_region_pops.size()) != nsims * ndists) {
        std::ostringstream msg;
        msg << "flattened_all_region_pops has wrong size. "
            << "size=" << flattened_all_region_pops.size()
            << ", expected=" << nsims * ndists;
        fail(msg.str());
    }

    if (static_cast<int>(flattened_all_region_order_added.size()) != nsims * ndists) {
        std::ostringstream msg;
        msg << "flattened_all_region_order_added has wrong size. "
            << "size=" << flattened_all_region_order_added.size()
            << ", expected=" << nsims * ndists;
        fail(msg.str());
    }

    if (num_forest_edge_bit_words_per_plan < 0) {
        fail("num_forest_edge_bit_words_per_plan is negative.");
    }

    if (num_forest_edge_bit_words_per_plan > 0) {
        std::size_t const expected_forest_words =
            static_cast<std::size_t>(nsims) *
            static_cast<std::size_t>(num_forest_edge_bit_words_per_plan);

        if (flattened_all_forest_edge_bits.size() != expected_forest_words) {
            std::ostringstream msg;
            msg << "flattened_all_forest_edge_bits has wrong size. "
                << "size=" << flattened_all_forest_edge_bits.size()
                << ", expected=" << expected_forest_words;
            fail(msg.str());
        }
    }

    if (static_cast<int>(map_params.pop.size()) != V) {
        std::ostringstream msg;
        msg << "map_params.pop has wrong size. "
            << "pop.n_elem=" << map_params.pop.size()
            << ", V=" << V;
        fail(msg.str());
    }

    int map_total_pop = 0;
    for (int v = 0; v < V; ++v) {
        map_total_pop += static_cast<int>(map_params.pop[v]);
    }

    for (int i = 0; i < nsims; ++i) {
        Plan *plan = plan_ptr_vec[i].get();

        if (plan == nullptr) {
            std::ostringstream msg;
            msg << "plan_ptr_vec[" << i << "] is null.";
            fail(msg.str());
        }

        auto fail_plan = [&](std::string const &msg) {
            std::ostringstream full;
            full << "Invalid plan in PlanEnsemble.\n";
            full << "Where: " << where << "\n";
            full << "Plan index: " << i << "\n";
            full << msg << "\n";
            full << plan->debug_string(true);
            throw std::runtime_error(full.str());
        };

        if (static_cast<int>(plan->region_ids.size()) != V) {
            std::ostringstream msg;
            msg << "region_ids has wrong size. "
                << "size=" << plan->region_ids.size()
                << ", expected V=" << V;
            fail_plan(msg.str());
        }

        if (static_cast<int>(plan->region_sizes.size()) != ndists) {
            std::ostringstream msg;
            msg << "region_sizes has wrong size. "
                << "size=" << plan->region_sizes.size()
                << ", expected ndists=" << ndists;
            fail_plan(msg.str());
        }

        if (static_cast<int>(plan->region_pops.size()) != ndists) {
            std::ostringstream msg;
            msg << "region_pops has wrong size. "
                << "size=" << plan->region_pops.size()
                << ", expected ndists=" << ndists;
            fail_plan(msg.str());
        }

        if (static_cast<int>(plan->region_added_order.size()) != ndists) {
            std::ostringstream msg;
            msg << "region_added_order has wrong size. "
                << "size=" << plan->region_added_order.size()
                << ", expected ndists=" << ndists;
            fail_plan(msg.str());
        }

        if (plan->num_regions < 1 || plan->num_regions > ndists) {
            std::ostringstream msg;
            msg << "num_regions is invalid. "
                << "num_regions=" << plan->num_regions
                << ", ndists=" << ndists;
            fail_plan(msg.str());
        }

        std::vector<int> counted_region_vertices(ndists, 0);
        std::vector<int> counted_region_pops(ndists, 0);

        for (int v = 0; v < V; ++v) {
            int const region_id = static_cast<int>(plan->region_ids[v]);

            if (region_id < 0 || region_id >= plan->num_regions) {
                std::ostringstream msg;
                msg << "Invalid region id. "
                    << "vertex=" << v
                    << ", region_id=" << region_id
                    << ", num_regions=" << plan->num_regions;
                fail_plan(msg.str());
            }

            ++counted_region_vertices[region_id];
            counted_region_pops[region_id] += static_cast<int>(map_params.pop[v]);
        }

        int total_counted_vertices = 0;
        int total_counted_pop = 0;

        for (int region_id = 0; region_id < plan->num_regions; ++region_id) {
            int const stored_size = static_cast<int>(plan->region_sizes[region_id]);
            int const stored_pop = static_cast<int>(plan->region_pops[region_id]);

            total_counted_vertices += counted_region_vertices[region_id];
            total_counted_pop += counted_region_pops[region_id];

            if (stored_pop != counted_region_pops[region_id]) {
                std::ostringstream msg;
                msg << "Stored region pop does not match counted pop. "
                    << "region_id=" << region_id
                    << ", stored_pop=" << stored_pop
                    << ", counted_pop=" << counted_region_pops[region_id];
                fail_plan(msg.str());
            }

            if (stored_size <= 0) {
                std::ostringstream msg;
                msg << "Active region has nonpositive size. "
                    << "region_id=" << region_id
                    << ", stored_size=" << stored_size;
                fail_plan(msg.str());
            }
        }

        if (total_counted_vertices != V) {
            std::ostringstream msg;
            msg << "Total counted vertices does not equal V. "
                << "total_counted_vertices=" << total_counted_vertices
                << ", V=" << V;
            fail_plan(msg.str());
        }

        if (total_counted_pop != map_total_pop) {
            std::ostringstream msg;
            msg << "Total counted pop does not equal map total pop. "
                << "total_counted_pop=" << total_counted_pop
                << ", map_total_pop=" << map_total_pop;
            fail_plan(msg.str());
        }

        for (int region_id = plan->num_regions; region_id < ndists; ++region_id) {
            int const stored_size = static_cast<int>(plan->region_sizes[region_id]);
            int const stored_pop = static_cast<int>(plan->region_pops[region_id]);

            if (stored_size != 0) {
                std::ostringstream msg;
                msg << "Inactive region has nonzero size. "
                    << "region_id=" << region_id
                    << ", stored_size=" << stored_size
                    << ", num_regions=" << plan->num_regions;
                fail_plan(msg.str());
            }

            if (stored_pop != 0) {
                std::ostringstream msg;
                msg << "Inactive region has nonzero pop. "
                    << "region_id=" << region_id
                    << ", stored_pop=" << stored_pop
                    << ", num_regions=" << plan->num_regions;
                fail_plan(msg.str());
            }
        }

        if (check_forest && num_forest_edge_bit_words_per_plan > 0) {
            std::string const forest_msg =
                std::string(where) +
                ", PlanEnsemble::check_all_plans_valid, plan i = " +
                std::to_string(i);

            plan->check_forest_integrity(
                map_params.graph_edge_index,
                forest_msg
            );
        }
    }
}

// Reorders all the plans in the vector by order a region was split
//
// Takes a vector of plans and uses the vector of dummy plans to reorder
// each of the plans by the order a region was split.
//
//
// @title Reorders all the plans in the vector by order a region was split
//
// @param pool A threadpool for multithreading
// @param plans_vec A vector of plans
// @param dummy_plans_vec A vector of dummy plans
//
// @details Modifications
//    - Each plan in the `plans_vec` object is reordered by when the region was split
//    - Each plan is a shallow copy of the plans in `plans_vec`
//
// @noRd
// @keywords internal
void reorder_all_plans(RcppThread::ThreadPool &pool,
                       std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
                       std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec) {

    int M = (int)plan_ptrs_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&](int i) {
        // reorder every plan
        plan_ptrs_vec.at(i)->reorder_plan_by_oldest_split(*dummy_plan_ptrs_vec.at(i));
    });

    // Wait for all the threads to finish
    pool.wait();

    return;
}

std::vector<std::unique_ptr<TreeSplitter>>
get_tree_splitter_ptrs(MapParams const &map_params, SplittingMethodType const splitting_method,
                       SamplingSpace const sampling_space,
                       Rcpp::List const &control, int const nsims, int const num_threads) {
    int V = map_params.V;
    double target = map_params.target;

    // create the pointer
    std::vector<std::unique_ptr<TreeSplitter>> tree_splitters_ptr_vec;
    tree_splitters_ptr_vec.reserve(nsims);

    if (splitting_method == SplittingMethodType::NaiveTopK) {
        // set splitting k to -1
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), num_threads,
                        [V] { return std::make_unique<NaiveTopKSplitter>(V, 1); });
    } else if (splitting_method == SplittingMethodType::UnifValid) {
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), num_threads,
                        [&map_graph = map_params.map_graph] { return std::make_unique<UniformValidSplitter>(map_graph); });
    } else if (splitting_method == SplittingMethodType::ExpBiggerAbsDev) {
        double alpha = Rcpp::as<double>(control["splitting_alpha"]);
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), num_threads,
                        [&map_graph = map_params.map_graph, alpha, target] {
                            return std::make_unique<ExpoWeightedSplitter>(map_graph, alpha, target);
                        });
    } else if (splitting_method == SplittingMethodType::ExpSmallerAbsDev) {
        double alpha = Rcpp::as<double>(control["splitting_alpha"]);
        std::generate_n(
            std::back_inserter(tree_splitters_ptr_vec), num_threads, [&map_graph = map_params.map_graph, alpha, target] {
                return std::make_unique<ExpoWeightedSmallerDevSplitter>(map_graph, alpha, target);
            });
    } else if (splitting_method == SplittingMethodType::Constraint) {
        int const ndists = map_params.ndists;
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), num_threads, [&map_graph = map_params.map_graph, V, ndists] {
            return std::make_unique<ConstraintSplitter>(map_graph, V, ndists);
        });
    } else if (splitting_method == SplittingMethodType::Experimental) {
        double epsilon = Rcpp::as<double>(control["splitting_epsilon"]);
        std::generate_n(std::back_inserter(tree_splitters_ptr_vec), num_threads,
                        [&map_graph = map_params.map_graph, epsilon, target] {
                            return std::make_unique<ExperimentalSplitter>(map_graph, epsilon, target);
                        });
    } else {
        throw Rcpp::exception("Invalid Splitting Method!");
    }

    return tree_splitters_ptr_vec;
}

std::vector<bool> vector_tree_to_edge_vector(
    GraphEdgeIndex const &edge_index,
    Tree const &tree
){
    // make a boolean vector of length edge_vec
    std::vector<bool> edge_vec(edge_index.num_edges, false);
    // Now we walk through the tree and mark edges in edge_vec as true
    for (int v = 0; v < edge_index.V; v++)
    {
        for (int const u: tree[v])
        {
            // make this edge as true 
            edge_vec[edge_index.get_edge_id(v, u)] = true;
        }
    }

    return edge_vec;
}

Rcpp::List graph_edge_index_to_list(
    GraphEdgeIndex const &edge_index
) {
    Rcpp::List edge_list(edge_index.num_edges);

    for (int edge_id = 0; edge_id < edge_index.num_edges; ++edge_id) {
        auto const [u, v] = edge_index.get_edge_endpoints(
            static_cast<EdgeID>(edge_id)
        );

        edge_list[edge_id] = Rcpp::IntegerVector::create(
            static_cast<int>(u) + 1,
            static_cast<int>(v) + 1
        );
    }

    return edge_list;
}

SMCDiagnostics::SMCDiagnostics(SamplingSpace const sampling_space,
                               SplittingMethodType const splitting_method_type,
                               SplittingSizeScheduleType const splitting_schedule_type,
                               std::vector<bool> const &merge_split_step_vec, int const V,
                               int const nsims, int const ndists, int const total_seats,
                               int const initial_num_regions, int const total_smc_steps,
                               int const total_ms_steps, 
                               bool const estimated_unbiased_normalizing_constant,
                               int const diagnostic_level,
                               bool const splitting_all_the_way, bool const split_district_only)
    : diagnostic_level(diagnostic_level), total_steps(total_smc_steps + total_ms_steps),
      log_wgt_stddevs(total_smc_steps), acceptance_rates(total_steps),
      nunique_parents(total_smc_steps), nunique_plans(total_steps), n_eff(total_smc_steps),
      num_merge_split_attempts_vec(total_ms_steps),
      cut_k_values(sampling_space == SamplingSpace::GraphSpace ? total_steps : 0),
      tries_before_extra_particle(estimated_unbiased_normalizing_constant ? total_smc_steps : 0),
      smc_step_parameter_estimation_times(total_smc_steps),
      smc_split_times(total_smc_steps),
      smc_weight_times(total_smc_steps),
      ms_step_parameter_estimation_times(total_ms_steps),
      ms_step_times(total_ms_steps),
      wilson_call_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      md_selection_times(perf_config::track_granular_times ? total_smc_steps : 0),
      plan_updating_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      hard_constraint_split_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      total_plan_smc_split_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_smc_steps : 0
      ),
      get_valid_pairs_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      plan_scores_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      region_scores_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      log_tau_times(perf_config::track_granular_times ? total_smc_steps + total_ms_steps : 0),
      retro_splitting_prob_times(perf_config::track_granular_times ? total_smc_steps : 0),
      total_plan_smc_weight_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_smc_steps : 0
      ),
      selecting_merge_pair_times(perf_config::track_granular_times ? total_ms_steps : 0),
      eff_boundary_times(perf_config::track_granular_times ? total_ms_steps : 0),
      total_plan_mcmc_times(
        perf_config::track_granular_times ? nsims : 0,
        perf_config::track_granular_times ? total_ms_steps : 0
      )
       {
    // Level 1 Diagnostics. Not too big relative to plan size
    log_incremental_weights_mat = arma::dmat(
        nsims, total_smc_steps, arma::fill::none); // entry [i][s] is the log unnormalized
                                                   // weight of particle i AFTER split s
    draw_tries_mat =
        Rcpp::IntegerMatrix(nsims, total_steps); // Entry [i][s] is the number of tries it took
                                                 // to form particle i on split s
    parent_index_mat = Rcpp::IntegerMatrix(
        nsims,
        total_smc_steps); // Entry [i][s] is the index of the parent of particle i at split s
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of
    // successful merge splits performed for plan i on merge split round s
    merge_split_successes_mat =
        Rcpp::IntegerMatrix(total_ms_steps > 0 ? nsims : 1, total_ms_steps);
    // counts the size of the trees
    tree_sizes_mat = Rcpp::IntegerMatrix(total_seats, total_steps);
    successful_tree_sizes_mat = Rcpp::IntegerMatrix(total_seats, total_steps);

    // Level 2
    parent_unsuccessful_tries_mat = Rcpp::IntegerMatrix(nsims, total_smc_steps);

    bool diagnostic_mode = diagnostic_level == 1;
    // level 3
    all_steps_plan_region_ids_list.reserve(diagnostic_mode ? total_steps : 0);
    all_steps_forests_adj_list.resize(
        (diagnostic_mode && sampling_space != SamplingSpace::GraphSpace) ? total_steps : 0);
    all_steps_linking_edge_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::LinkingEdgeSpace) ? total_steps
                                                                               : 0);
    all_steps_valid_region_sizes_to_split.resize(diagnostic_mode ? total_smc_steps : 0);
    all_steps_valid_split_region_sizes.resize(diagnostic_mode ? total_smc_steps : 0);

    // Store size at every step but last one if needed
    int plan_dval_list_size = (diagnostic_mode & !split_district_only) ? total_steps - 1 : 0;
    if (!splitting_all_the_way)
        plan_dval_list_size++;

    region_sizes_mat_list.reserve(plan_dval_list_size);

    // If diagnostic mode track vertex region ids from every round
    if (diagnostic_mode) {
        // The number of regions starts at 1
        int curr_num_regions = initial_num_regions;
        for (size_t i = 0; i < total_steps; i++) {
            all_steps_plan_region_ids_list.emplace_back(V, nsims);
            // Create V by nsims matrix for the plan
            // This is a vector where every entry is a V by nsims Rcpp::IntegerMatrix

            // increase number of regions by 1 if that step is an smc one
            if (!merge_split_step_vec.at(i))
                curr_num_regions++;

            // If not doing district only splits, and its not the final one or
            // we're only doing partial plans then make size matrix
            if (!split_district_only && (i < total_steps - 1 || !splitting_all_the_way)) {
                // This is number of regions by nsims
                region_sizes_mat_list.emplace_back(curr_num_regions, nsims);
            }
        }
    }
}

void SMCDiagnostics::add_full_step_diagnostics(
    int const total_steps, bool const splitting_all_the_way, int const step_num,
    int const merge_split_step_num, int const smc_step_num, bool const is_smc_step,
    SamplingSpace const sampling_space, RcppThread::ThreadPool &pool,
    PlanEnsemble &plan_ensemble, PlanEnsemble &new_plans_ensemble,
    SplittingSchedule const &splitting_schedule) {
    // if(diagnostic_mode){ // record if in diagnostic mode and generalized splits
    //  reorder the plans by oldest split if either we'vxe done any merge split or
    //  its generalized region splits

    bool const split_district_only =
        splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;
    int const nsims = plan_ensemble.nsims;

    // if smc step update splitting step info
    if (is_smc_step) {
        int current_num_regions = plan_ensemble.plan_ptr_vec[0]->num_regions;
        // save the acceptable split sizes
        for (int region_size = 1;
             region_size <= splitting_schedule.total_seats - current_num_regions + 2;
             region_size++) {
            if (splitting_schedule.valid_split_region_sizes[region_size]) {
                all_steps_valid_split_region_sizes[smc_step_num].push_back(region_size);
            }
            if (splitting_schedule.valid_region_sizes_to_split[region_size]) {

                all_steps_valid_region_sizes_to_split[smc_step_num].push_back(region_size);
                ;
            }
        }
    }

    if (merge_split_step_num > 0 || !split_district_only) {
        reorder_all_plans(pool, plan_ensemble.plan_ptr_vec, new_plans_ensemble.plan_ptr_vec);
    }

    // Copy the vertex plan matrix
    all_steps_plan_region_ids_list.at(step_num) = plan_ensemble.get_R_plans_matrix();

    // store the
    if (!(sampling_space == SamplingSpace::GraphSpace)) {
        all_steps_forests_adj_list.at(step_num).reserve(nsims);
        for (size_t i = 0; i < nsims; i++) {
            // add the forests from each plan at this step
            all_steps_forests_adj_list.at(step_num).push_back(
                plan_ensemble.plan_ptr_vec[i]->get_forest_adj());
        }
        if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
            for (size_t i = 0; i < nsims; i++) {
                // add the forests from each plan at this step
                all_steps_linking_edge_list.at(step_num).push_back(
                    plan_ensemble.plan_ptr_vec[i]->get_linking_edges());
            }
        }
    }

    // Copy the sizes if neccesary
    if (!split_district_only && (step_num < total_steps - 1 || !splitting_all_the_way)) {
        region_sizes_mat_list.at(step_num) = plan_ensemble.get_R_sizes_matrix(pool);
    }

    return;
}

void SMCDiagnostics::add_diagnostics_to_out_list(Rcpp::List &out) {
    // make parent index 1 indexed in place
    std::transform(parent_index_mat.begin(), parent_index_mat.end(), parent_index_mat.begin(),
                   [](int x) { return x + 1; });

    // add granular time info 
    Rcpp::List granular_timing;
    if constexpr (perf_config::track_granular_times){
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = true,
            Rcpp::_["wilson_call_times"] = wilson_call_times,
            Rcpp::_["multidistrict_selection_times"] = md_selection_times,
            Rcpp::_["plan_updating_times"] = plan_updating_times,
            Rcpp::_["hard_constraint_split_times"] = hard_constraint_split_times,
            Rcpp::_["total_plan_smc_split_times"] = total_plan_smc_split_times,
            Rcpp::_["getting_valid_pairs_times"] = get_valid_pairs_times,
            Rcpp::_["computing_plan_scores_times"] = plan_scores_times,
            Rcpp::_["computing_region_scores_times"] = region_scores_times,
            Rcpp::_["computing_spanning_tree_count_times"] = log_tau_times,   
            Rcpp::_["computing_retro_splitting_prob_times"] = retro_splitting_prob_times,     
            Rcpp::_["total_plan_smc_weight_times"] = total_plan_smc_weight_times,
            Rcpp::_["selecting_merge_pair"] = selecting_merge_pair_times,
            Rcpp::_["eff_boundary_times"] = eff_boundary_times,
            Rcpp::_["total_plan_mcmc_times"] = total_plan_mcmc_times
        );

    }else{
        granular_timing = Rcpp::List::create(
            Rcpp::_["granular_time_tracked"] = false
        );
    }

    // optional add special timing 
    if constexpr (perf_config::special_timing){
        Rcpp::List special_timing_list = Rcpp::List::create(
        );

        granular_timing["special_timing_list"] = special_timing_list;
    }
    
    out["acceptance_rates"] = acceptance_rates;
    out["draw_tries_mat"] = draw_tries_mat;
    out["parent_index"] = parent_index_mat;
    out["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat;
    out["step_n_eff"] = n_eff;
    out["nunique_parent_indices"] = nunique_parents;
    out["nunique_plans"] = nunique_plans;
    out["tree_sizes"] = tree_sizes_mat;
    out["successful_tree_sizes"] = successful_tree_sizes_mat;
    out["log_weight_stddev"] = log_wgt_stddevs;
    out["cut_k_vals"] = cut_k_values;
    out["log_incremental_weights_mat"] = log_incremental_weights_mat;
    out["ms_step_counts"] = num_merge_split_attempts_vec;
    out["merge_split_success_mat"] = merge_split_successes_mat;
    out["tries_before_extra_particle"] = tries_before_extra_particle;
    out["smc_step_parameter_estimation_times"] = smc_step_parameter_estimation_times;
    out["smc_split_times"] = smc_split_times;
    out["smc_weight_times"] = smc_weight_times;
    out["ms_step_parameter_estimation_times"] = ms_step_parameter_estimation_times;
    out["ms_step_times"] = ms_step_times;
    out["granular_times"] = granular_timing;
    out["region_ids_mat_list"] = all_steps_plan_region_ids_list;
    out["region_seats_mat_list"] = region_sizes_mat_list;
    out["forest_adjs_list"] = all_steps_forests_adj_list;
    out["linking_edges_list"] = all_steps_linking_edge_list;
    out["valid_split_region_sizes_list"] = all_steps_valid_split_region_sizes;
    out["valid_region_sizes_to_split_list"] = all_steps_valid_region_sizes_to_split;


    return;
}


