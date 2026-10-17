/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2024/10
 * Purpose: Helper functions for all redist algorithm types
 ********************************************************/

#include <iterator>
#include <numeric>

#include "redist_alg_helpers.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

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



int check_ensemble_dims(MapParams const &map_params, int const nsims,
                        SamplingSpace const sampling_space) {
    if (nsims <= 0) {
        throw Rcpp::exception("Tried to create a plan ensemble with no plans!\n");
    }

    constexpr std::int64_t max_elements =
        static_cast<std::int64_t>(std::numeric_limits<int>::max());

    auto check_buffer = [&](int const per_plan, char const *what) {
        if (per_plan <= 0) return; // buffer is unused for this sampling space
        std::int64_t const total =
            static_cast<std::int64_t>(per_plan) * static_cast<std::int64_t>(nsims);
        if (total > max_elements) {
            std::ostringstream oss;
            oss << "Cannot create an ensemble of " << nsims << " plans: the " << what
                << " buffer would need " << total << " elements, more than the "
                << max_elements << " an int offset can address.\n"
                << "Reduce the number of plans to at most " << (max_elements / per_plan)
                << ".\n";
            throw Rcpp::exception(oss.str().c_str());
        }
    };

    check_buffer(map_params.V, "plan region id");
    check_buffer(map_params.ndists, "per region");
    if (sampling_space == SamplingSpace::ForestSpace ||
        sampling_space == SamplingSpace::LinkingEdgeSpace) {
        check_buffer(map_params.num_edge_bit_words, "forest edge bit");
    }

    return nsims;
}

// creates plan ensemble of blank plans
PlanEnsemble::PlanEnsemble(MapParams const &map_params, int const total_pop, int const nsims,
                           SamplingSpace const sampling_space, RcppThread::ThreadPool &pool,
                           int const verbosity)
    // `nsims` is the first member, so this check runs before any buffer below is sized
    : nsims(check_ensemble_dims(map_params, nsims, sampling_space)), V(map_params.V),
      ndists(map_params.ndists),
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
            throw std::runtime_error(
                "PlanEnsemble: can't create blank plans in the " +
                sampling_space_to_str(sampling_space) +
                " sampling space. Supported spaces are graph, LCT graph, spanning forest, "
                "and linking edge.\n");
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
                           int const verbosity,
                           InitialAugmentedSamples const &initial_augmented_samples)
    // `nsims` is the first member, so this check runs before any buffer below is sized
    : nsims(check_ensemble_dims(map_params, nsims, sampling_space)), V(map_params.V),
      ndists(map_params.ndists),
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
    int const num_threads = get_num_threads(pool);
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
        if (rng_states.size() < static_cast<std::size_t>(num_threads)) {
            std::ostringstream oss;
            oss << "Not enough RNG states.\n";
            oss << "rng_states.size()="
                << rng_states.size() << "\n";
            oss << "n_threads="
                << num_threads << "\n";

            throw std::runtime_error(oss.str());
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

        // the saved forest and linking edges for this plan if there are any
        EdgeBitWord const *initial_forest_words =
            initial_augmented_samples.forest_words.empty()
                ? nullptr
                : initial_augmented_samples.forest_words.data() +
                      static_cast<std::size_t>(i) * num_forest_edge_bit_words_per_plan;
        std::vector<LinkingEdge> const *initial_linking_edges =
            initial_augmented_samples.linking_edges.empty()
                ? nullptr
                : &initial_augmented_samples.linking_edges[i];

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
                rng_states[thread_id], initial_forest_words);
        } else if (use_linking_edge_space) {
            plan_ptr_vec[i] = std::make_unique<LinkingEdgePlan>(
                ndists, num_regions, map_params.pop, plan_region_ids, plan_sizes, plan_pops,
                plan_region_order_added, plan_forest_edge_bits,
                *tree_ptr, ust_sampler, plan_multigraph, region_graph,
                rng_states[thread_id], initial_forest_words, initial_linking_edges);
        } else {
            throw std::runtime_error(
                "PlanEnsemble: can't create partial plans in the " +
                sampling_space_to_str(sampling_space) +
                " sampling space. Supported spaces are graph, LCT graph, spanning forest, "
                "and linking edge.\n");
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
    // widen to R's 4-byte integers and shift to 1-indexed in a single pass
    std::transform(flattened_all_plans.begin(), flattened_all_plans.end(), plan_mat.begin(),
                   [](RegionID x) { return static_cast<int>(x) + 1; });
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

std::vector<PlanTallyEntry> PlanEnsemble::get_unique_plan_tally(
    RcppThread::ThreadPool &pool) const {
    if (nsims == 0) {
        return {};
    }

    int const num_regions = plan_ptr_vec[0]->num_regions;
    int const num_threads = get_num_threads(pool);

    // A plan is just a run of `V` region ids in the flattened buffer, and a
    // raw pointer into it already indexes the way the shared helpers expect.
    auto const plan_ptr = [&](int const plan_idx) {
        return flattened_all_plans.data() +
               static_cast<std::size_t>(plan_idx) * static_cast<std::size_t>(V);
    };

    std::vector<std::vector<int>> label_scratches(num_threads,
                                                  std::vector<int>(num_regions));

    // Stores only {hash, plan index}: 16 bytes per plan.
    std::vector<std::pair<std::uint64_t, int>> hashed_plans(nsims);

    // Trick to give each thread a unique id
    // We need the extra steps to avoid the problem where
    // only some threads persist from previous calls but the counter resets
    // resulting in multiple threads with the same id
    static std::atomic<int> global_generation_counter{0};
    int const generation = global_generation_counter.fetch_add(1, std::memory_order_relaxed);
    std::atomic<int> thread_id_counter{0};

    // First pass: compute each plan's canonical hash.
    pool.parallelFor(0, nsims, [&](int const plan_idx) {
        static thread_local int thread_generation_counter = -1;
        static thread_local int thread_id = -1;

        // check if the thread id was generated this function call
        if (thread_generation_counter != generation) {
            // if not then give it a new id
            thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
            thread_generation_counter = generation;
        }

        hashed_plans[plan_idx] = {
            canonical_plan_hash(plan_ptr(plan_idx), V, num_regions,
                                label_scratches[thread_id]),
            plan_idx
        };
    });

    pool.wait();

    // Put plans with matching hashes next to each other.
    std::sort(hashed_plans.begin(), hashed_plans.end());

    std::vector<int> plan1_to_plan2(num_regions);
    std::vector<int> plan2_to_plan1(num_regions);

    std::vector<PlanTallyEntry> tally;
    tally.reserve(nsims);

    // Slots in `tally` for the group of plans sharing the current hash.
    // Usually holds one entry; it grows only when distinct plans happen to
    // share a 64-bit hash, which is why matches are confirmed exactly below.
    std::vector<int> representative_slots;
    representative_slots.reserve(2);

    std::size_t group_begin = 0;

    while (group_begin < hashed_plans.size()) {
        std::size_t group_end = group_begin + 1;

        while (group_end < hashed_plans.size() &&
               hashed_plans[group_end].first == hashed_plans[group_begin].first) {
            ++group_end;
        }

        representative_slots.clear();

        for (std::size_t i = group_begin; i < group_end; ++i) {
            int const candidate_idx = hashed_plans[i].second;
            bool matched = false;

            for (int const slot : representative_slots) {
                if (plans_canonically_equal(plan_ptr(candidate_idx),
                                            plan_ptr(tally[slot].representative), V,
                                            num_regions, plan1_to_plan2,
                                            plan2_to_plan1)) {
                    ++tally[slot].count;
                    matched = true;
                    break;
                }
            }

            if (!matched) {
                if (!representative_slots.empty()) {
                    // genuine hash collision, chain it off the last one
                    tally[representative_slots.back()].next_same_hash =
                        static_cast<int>(tally.size());
                }
                representative_slots.push_back(static_cast<int>(tally.size()));
                tally.push_back(PlanTallyEntry{candidate_idx, 1, -1});
            }
        }

        group_begin = group_end;
    }

    return tally;
}


int PlanEnsemble::count_unique_plans(RcppThread::ThreadPool &pool) const {
    return static_cast<int>(get_unique_plan_tally(pool).size());
}

void PlanEnsemble::release_all_but_plan_ids() {
    // Swapping with a fresh empty vector is what actually hands the memory
    // back: `shrink_to_fit` is only a non-binding request, whereas the
    // temporary here takes ownership of the buffer and deallocates it when it
    // dies at the end of the statement.
    auto release = [](auto &buffer) { std::decay_t<decltype(buffer)>().swap(buffer); };

    // The plans are views into the buffers below, so they have to go first.
    release(plan_ptr_vec);
    release(flattened_all_region_sizes);
    release(flattened_all_region_pops);
    release(flattened_all_region_order_added);
    release(flattened_all_forest_edge_bits);
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

Rcpp::LogicalMatrix PlanEnsemble::get_R_forest_matrix(int const num_edges,
                                                      RcppThread::ThreadPool &pool) {
    if (num_forest_edge_bit_words_per_plan == 0) {
        throw Rcpp::exception(
            "Spanning forests are only stored in forest or linking edge space!\n");
    }
    if (compute_num_edge_bit_words(num_edges) != num_forest_edge_bit_words_per_plan ||
        flattened_all_forest_edge_bits.size() !=
            static_cast<std::size_t>(nsims) * num_forest_edge_bit_words_per_plan) {
        throw Rcpp::exception("Forest edge bits do not match the number of graph edges!\n");
    }

    Rcpp::LogicalMatrix forest_mat(num_edges, nsims);
    // take the raw pointer up front so the threads never touch the R API
    int *const forest_mat_ptr = LOGICAL(forest_mat);

    // plan i's forest is the i-th run of words in the flattened buffer, with
    // edge id e stored in bit e % 64 of word e / 64
    pool.parallelFor(0, nsims, [&](int i) {
        EdgeBitWord const *plan_words =
            flattened_all_forest_edge_bits.data() +
            static_cast<std::size_t>(i) * num_forest_edge_bit_words_per_plan;
        int *const plan_col = forest_mat_ptr + static_cast<std::size_t>(i) * num_edges;

        for (int edge_id = 0; edge_id < num_edges; ++edge_id) {
            plan_col[edge_id] = static_cast<int>(
                (plan_words[edge_id / EDGE_BITS_PER_WORD] >> (edge_id % EDGE_BITS_PER_WORD)) &
                EdgeBitWord{1});
        }
    });
    pool.wait();

    return forest_mat;
}

Rcpp::DataFrame PlanEnsemble::get_R_linking_edges(GraphEdgeIndex const &edge_index) {
    if (sampling_space != SamplingSpace::LinkingEdgeSpace) {
        throw Rcpp::exception("Linking edges are only stored in linking edge space!\n");
    }

    // Plans are not assumed to all have the same number of linking edges
    std::size_t total_linking_edges = 0;
    for (int i = 0; i < nsims; ++i) {
        total_linking_edges += plan_ptr_vec[i]->get_linking_edges_ref().size();
    }

    Rcpp::IntegerVector draw(total_linking_edges);
    Rcpp::IntegerVector edge_id(total_linking_edges);
    Rcpp::IntegerVector vertex1(total_linking_edges);
    Rcpp::IntegerVector vertex2(total_linking_edges);
    Rcpp::NumericVector log_prob(total_linking_edges);

    std::size_t row = 0;
    std::vector<std::pair<EdgeID, LinkingEdge const *>> plan_edges;
    for (int i = 0; i < nsims; ++i) {
        // sort by edge id so the output only depends on the linking edge set
        plan_edges.clear();
        for (auto const &linking_edge : plan_ptr_vec[i]->get_linking_edges_ref()) {
            plan_edges.emplace_back(
                edge_index.get_edge_id(linking_edge.vertex1, linking_edge.vertex2),
                &linking_edge);
        }
        std::sort(plan_edges.begin(), plan_edges.end(),
                  [](auto const &a, auto const &b) { return a.first < b.first; });

        for (auto const &[plan_edge_id, linking_edge] : plan_edges) {
            draw[row] = i + 1;
            edge_id[row] = static_cast<int>(plan_edge_id) + 1;
            vertex1[row] = linking_edge->vertex1;
            vertex2[row] = linking_edge->vertex2;
            log_prob[row] = linking_edge->valid_log_prob ? linking_edge->log_prob : NA_REAL;
            ++row;
        }
    }

    return Rcpp::DataFrame::create(Rcpp::_["draw"] = draw, Rcpp::_["edge_id"] = edge_id,
                                   Rcpp::_["vertex1"] = vertex1, Rcpp::_["vertex2"] = vertex2,
                                   Rcpp::_["log_prob"] = log_prob);
}

InitialAugmentedSamples get_initial_augmented_samples(Rcpp::List const &control,
                                                      MapParams const &map_params,
                                                      int const nsims) {
    InitialAugmentedSamples initial_augmented_samples;

    if (control.containsElementNamed("initial_spanning_forests")) {
        Rcpp::LogicalMatrix const forests = control["initial_spanning_forests"];
        int const num_edge_bit_words = map_params.num_edge_bit_words;
        initial_augmented_samples.forest_words.assign(
            static_cast<std::size_t>(nsims) * num_edge_bit_words, EdgeBitWord{0});

        // same packing as the ensemble, edge id e is bit e % 64 of word e / 64
        for (int i = 0; i < nsims; ++i) {
            EdgeBitWord *plan_words = initial_augmented_samples.forest_words.data() +
                                      static_cast<std::size_t>(i) * num_edge_bit_words;
            for (int edge_id = 0; edge_id < map_params.num_edges; ++edge_id) {
                if (forests(edge_id, i) == TRUE) {
                    plan_words[edge_id / EDGE_BITS_PER_WORD] |=
                        EdgeBitWord{1} << (edge_id % EDGE_BITS_PER_WORD);
                }
            }
        }
    }

    if (control.containsElementNamed("initial_linking_edges")) {
        Rcpp::List const linking_edges = Rcpp::as<Rcpp::List>(control["initial_linking_edges"]);
        Rcpp::IntegerVector const draw = linking_edges["draw"];
        Rcpp::IntegerVector const vertex1 = linking_edges["vertex1"];
        Rcpp::IntegerVector const vertex2 = linking_edges["vertex2"];

        initial_augmented_samples.linking_edges.resize(nsims);
        for (R_xlen_t k = 0; k < draw.size(); ++k) {
            // no log probability so it is marked invalid and recomputed when needed
            initial_augmented_samples.linking_edges.at(draw[k] - 1).emplace_back(vertex1[k],
                                                                                vertex2[k]);
        }
    }

    return initial_augmented_samples;
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
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity,
    InitialAugmentedSamples const &initial_augmented_samples) {
    if (num_regions == 1) {
        return std::make_unique<PlanEnsemble>(map_params,
            std::accumulate(map_params.pop.begin(), map_params.pop.end(), 0),
            nsims, sampling_space, pool, verbosity);
    } else {
        return std::make_unique<PlanEnsemble>(map_params, splitting_schedule, num_regions,
                                              nsims, sampling_space, plans_mat,
                                              region_sizes_mat, rng_states, pool, verbosity,
                                              initial_augmented_samples);
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
        throw std::runtime_error("get_tree_splitter_ptrs: no tree splitter for the " +
                                 splitting_method_to_str(splitting_method) +
                                 " splitting method.\n");
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

Rcpp::IntegerMatrix graph_edge_index_to_matrix(
    GraphEdgeIndex const &edge_index
) {
    Rcpp::IntegerMatrix edge_mat(edge_index.num_edges, 2);

    for (int edge_id = 0; edge_id < edge_index.num_edges; ++edge_id) {
        auto const [u, v] = edge_index.get_edge_endpoints(
            static_cast<EdgeID>(edge_id)
        );

        edge_mat(edge_id, 0) = static_cast<int>(u) + 1;
        edge_mat(edge_id, 1) = static_cast<int>(v) + 1;
    }

    return edge_mat;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix graph_edge_list(Rcpp::List const &adj_list) {
    Graph const g = list_to_graph(adj_list);
    return graph_edge_index_to_matrix(GraphEdgeIndex(g, count_undirected_edges(g)));
}

