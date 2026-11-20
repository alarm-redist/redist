#pragma once
#ifndef SMC_ALG_HELPERS_H
#define SMC_ALG_HELPERS_H


#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <Rcpp.h>

#include "advanced_types.h"
#include "base_plan_type.h"



class TreeSplitter;
class SplittingSchedule;
class RNGState;

namespace RcppThread {
    class ThreadPool;
}


/********************************************************
 * Canonical hashing and equality for plans.
 *
 * Two plans describe the same partition when one's region ids are a
 * permutation of the other's, so comparing them means comparing their
 * canonical labellings: the region of vertex 0 becomes label 0, the region of
 * the next vertex in a different region becomes 1, and so on.
 *
 * These are templated on how a plan is indexed, so the same code serves an
 * `Rcpp::IntegerMatrix` column, a `PlanVector` out of a `PlanEnsemble`, or a
 * plain `std::vector`. Canonical labelling also makes the result independent
 * of whether the source ids were 0- or 1-indexed.
 ********************************************************/

// Marks a region that has not been seen yet in a relabelling scratch vector.
inline constexpr int UNSEEN_REGION = -1;

/*
 * Final avalanche step of splitmix64.
 *
 * The accumulator below is a cheap multiply-xor chain that mixes the low bits
 * well but leaves structure in the high ones. Running the result through this
 * spreads every input bit across the whole word, which keeps the bucket
 * distribution flat.
 */
inline std::uint64_t mix64_hash(std::uint64_t x) noexcept {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

/*
 * Hash of a plan's canonical labelling.
 *
 * `label_scratch` needs room for `num_regions` entries and is reset here, so
 * a caller can hand the same buffer back on every call and the whole routine
 * allocates nothing. That is the point: hashing via a per plan string
 * allocates on every call, which serializes threads on the allocator.
 */
template <typename PlanID>
std::uint64_t canonical_plan_hash(PlanID const &region_ids, int const V,
                                  int const num_regions,
                                  std::vector<int> &label_scratch) {
    label_scratch.assign(num_regions, UNSEEN_REGION);

    std::uint64_t hash = 0xcbf29ce484222325ULL; // FNV-1a offset basis
    int next_label = 0;

    for (int v = 0; v < V; ++v) {
        int const region = static_cast<int>(region_ids[v]);

        int label = label_scratch[region];
        if (label == UNSEEN_REGION) {
            label = next_label++;
            label_scratch[region] = label;
        }

        hash ^= static_cast<std::uint64_t>(label) + 1ULL;
        hash *= 0x100000001b3ULL; // FNV-1a prime
    }

    return mix64_hash(hash);
}

/*
 * Whether two plans are the same partition, ignoring how the regions happen
 * to be numbered.
 *
 * Walks both at once building the bijection between their region ids, bailing
 * the moment the two disagree. Both scratch vectors need room for
 * `num_regions` entries and are reset here.
 *
 * Used to confirm a hash match rather than trust it, so counts stay exact
 * even in the vanishingly unlikely event of a collision.
 */
template <typename PlanIDA, typename PlanIDB>
bool plans_canonically_equal(PlanIDA const &a, PlanIDB const &b, int const V,
                             int const num_regions,
                             std::vector<int> &a_to_b_scratch,
                             std::vector<int> &b_to_a_scratch) {
    a_to_b_scratch.assign(num_regions, UNSEEN_REGION);
    b_to_a_scratch.assign(num_regions, UNSEEN_REGION);

    for (int v = 0; v < V; ++v) {
        int const region_a = static_cast<int>(a[v]);
        int const region_b = static_cast<int>(b[v]);

        int const mapped_a = a_to_b_scratch[region_a];
        int const mapped_b = b_to_a_scratch[region_b];

        if (mapped_a == UNSEEN_REGION && mapped_b == UNSEEN_REGION) {
            // neither side has been seen before, so pair them up
            a_to_b_scratch[region_a] = region_b;
            b_to_a_scratch[region_b] = region_a;
        } else if (mapped_a != region_b || mapped_b != region_a) {
            // contradicts a pairing established earlier
            return false;
        }
    }

    return true;
}


/*
 * One distinct plan found by `tally_unique_plans`.
 *
 * `next_same_hash` chains entries whose canonical hashes collide. Keeping the
 * link inside the entry means a bucket costs nothing beyond the map node: no
 * per bucket container is allocated, and chains are essentially always length
 * one anyway.
 */
struct PlanTallyEntry {
    int representative;  // first particle carrying this plan
    int count;           // how many particles carry it
    int next_same_hash;  // next entry with the same hash, or -1
};

/*
 * Tally how many particles share each distinct plan.
 *
 * `plan_at(i)` must return something indexable by vertex for particle `i`,
 * which lets this serve a plans matrix, a `PlanEnsemble`, or anything else
 * without copying the plans first.
 *
 * Exact rather than hash-only: a hash match is confirmed with
 * `plans_canonically_equal` before two particles are merged, so a collision
 * costs an extra comparison rather than silently undercounting.
 */
template <typename PlanAt>
std::vector<PlanTallyEntry> tally_unique_plans(PlanAt const &plan_at, int const nsims,
                                               int const V, int const num_regions) {
    std::vector<PlanTallyEntry> tally;
    tally.reserve(nsims);

    // canonical hash -> index of the first tally entry with that hash
    std::unordered_map<std::uint64_t, int> heads;
    heads.reserve(nsims * 2);

    std::vector<int> label_scratch(num_regions);
    std::vector<int> a_to_b(num_regions);
    std::vector<int> b_to_a(num_regions);

    for (int i = 0; i < nsims; ++i) {
        auto const plan = plan_at(i);
        std::uint64_t const hash = canonical_plan_hash(plan, V, num_regions, label_scratch);

        auto const inserted = heads.emplace(hash, static_cast<int>(tally.size()));

        if (inserted.second) {
            // first plan with this hash
            tally.push_back(PlanTallyEntry{i, 1, -1});
            continue;
        }

        // walk the chain looking for the same partition
        int slot = inserted.first->second;
        int last = -1;
        bool matched = false;

        while (slot != -1) {
            auto const other = plan_at(tally[slot].representative);
            if (plans_canonically_equal(plan, other, V, num_regions, a_to_b, b_to_a)) {
                ++tally[slot].count;
                matched = true;
                break;
            }
            last = slot;
            slot = tally[slot].next_same_hash;
        }

        if (!matched) {
            // a genuine hash collision between different partitions
            tally[last].next_same_hash = static_cast<int>(tally.size());
            tally.push_back(PlanTallyEntry{i, 1, -1});
        }
    }

    return tally;
}


Rcpp::List maximum_input_sizes();

/*
 * Convert zero-indxed R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const Rcpp::List &l);

/*
 * Creates a reindexing vector for the plan with two regions merged
 * does this by making region2_id map to region1_id. Nothing else
 * changes
 */
void set_merged_region_reindex_vec(int const num_regions, std::vector<int> &region_reindex_vec,
                                   int const region1_id, int const region2_id);

/*
 *  Reorders all the plans in the vector by order a region was split
 *
 *  Takes a vector of plans and uses the vector of dummy plans to reorder
 *  each of the plans by the order a region was split.
 *
 *
 *  @title Reorders all the plans in the vector by order a region was split
 *
 *  @param pool A threadpool for multithreading
 *  @param plan_ptrs_vec A vector of pointers to plans
 *  @param dummy_plans_vec A vector of pointers to dummy plans
 *
 *  @details Modifications
 *     - Each plan in the `plans_vec` object is reordered by when the region was split
 *     - Each plan is a shallow copy of the plans in `plans_vec`
 *
 */
void reorder_all_plans(RcppThread::ThreadPool &pool,
                       std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
                       std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec);

std::vector<std::unique_ptr<TreeSplitter>>
get_tree_splitter_ptrs(MapParams const &map_params, SplittingMethodType const splitting_method,
                       SamplingSpace const sampling_space,
                       Rcpp::List const &control, int const nsims, int const num_threads);

/*
 *  Checks that an ensemble of `nsims` plans can actually be addressed, throwing
 *  if not, and returns `nsims` unchanged so it can be called from a member
 *  initializer list.
 *
 *  Every plan attribute lives in one flat buffer indexed by `int` offsets (see
 *  `PlanAttribute`), so each buffer's element count has to fit in an `int`. The
 *  products are formed in 64 bits here because computing them the natural way,
 *  as `V * nsims`, is exactly the signed overflow this guards against - that
 *  would be undefined behaviour and would silently under-allocate rather than
 *  fail. R cannot hold a matrix of more than INT_MAX elements either, so a run
 *  past this limit could never return its plans in any case.
 */
int check_ensemble_dims(MapParams const &map_params, int const nsims,
                        SamplingSpace const sampling_space);

// Saved spanning forests and linking edges (from `return_augmented_samples`)
// used to initialize partial plans instead of drawing them at random. They
// are not checked. Empty members mean draw at random.
struct InitialAugmentedSamples {
    // num_edge_bit_words packed words per plan, plan i's words starting at
    // i * num_edge_bit_words
    std::vector<EdgeBitWord> forest_words;
    // the linking edges of each plan
    std::vector<std::vector<LinkingEdge>> linking_edges;
};

// Reads the saved spanning forests (a num_edges by nsims logical matrix) and
// linking edges (a long data frame with 1-indexed `draw` and 0-indexed
// `vertex1` and `vertex2` columns) out of the `control` list elements
// `initial_spanning_forests` and `initial_linking_edges` if present.
InitialAugmentedSamples get_initial_augmented_samples(Rcpp::List const &control,
                                                      MapParams const &map_params,
                                                      int const nsims);

// lightweight container for plans
class PlanEnsemble {

  public:
    // constructor for empty plans
    PlanEnsemble(MapParams const &map_params, int const total_pop, int const nsims,
                 SamplingSpace const sampling_space, RcppThread::ThreadPool &pool,
                 int const verbosity = 3);
    // constructor for non-empty starting plans
    PlanEnsemble(MapParams const &map_params, SplittingSchedule const &splitting_schedule,
                 int const num_regions, int const nsims, SamplingSpace const sampling_space,
                 Rcpp::IntegerMatrix const &plans_mat,
                 Rcpp::IntegerMatrix const &region_sizes_mat, std::vector<RNGState> &rng_states,
                 RcppThread::ThreadPool &pool, int const verbosity = 3,
                 InitialAugmentedSamples const &initial_augmented_samples = {});

    int const nsims;
    int const V;
    int const ndists;
    int const total_seats;
    SamplingSpace const sampling_space;
    std::vector<RegionID> flattened_all_plans;
    std::vector<RegionID> flattened_all_region_sizes;
    std::vector<int> flattened_all_region_pops;
    std::vector<int> flattened_all_region_order_added;
    // Empty unless sampling space is ForestSpace or LinkingEdgeSpace.
    int const num_forest_edge_bit_words_per_plan;
    std::vector<EdgeBitWord> flattened_all_forest_edge_bits;
    std::vector<std::unique_ptr<Plan>> plan_ptr_vec;

    // exports current plans to 1-indexed Rcpp matrix
    Rcpp::IntegerMatrix get_R_plans_matrix();
    // export current region sizes to Rcpp matrix
    Rcpp::IntegerMatrix get_R_sizes_matrix(RcppThread::ThreadPool &pool);
    // export current region populations to Rcpp matrix
    Rcpp::IntegerMatrix get_region_pops_matrix(RcppThread::ThreadPool &pool);
    // export current spanning forests to a num_edges by nsims logical matrix
    // where entry (e, i) is whether edge id e is in plan i's forest. Only
    // valid in forest or linking edge space.
    Rcpp::LogicalMatrix get_R_forest_matrix(int const num_edges, RcppThread::ThreadPool &pool);
    // export current linking edges in long format as a data frame with one row
    // per linking edge and columns `draw` (1-indexed plan), `edge_id`
    // (1-indexed canonical edge id), `vertex1` and `vertex2` (the 0-indexed
    // endpoints as stored in the plan) and `log_prob` (the cached log
    // probability the edge is chosen when splitting the merged regions, NA if
    // stale). Rows are grouped by draw and sorted by edge id within a draw.
    // Plans may have different numbers of linking edges. Only valid in
    // linking edge space.
    Rcpp::DataFrame get_R_linking_edges(GraphEdgeIndex const &edge_index);
    // One entry per distinct plan in the ensemble: a representative particle,
    // how many particles carry that plan, and the hash collision chain.
    std::vector<PlanTallyEntry> get_unique_plan_tally(RcppThread::ThreadPool &pool) const;
    // counts the number of unique plans in the ensemble
    int count_unique_plans(RcppThread::ThreadPool &pool) const;

    // Frees every buffer except the plan region ids, leaving the ensemble
    // valid only for `get_R_plans_matrix`. The R plan matrix needs 4 bytes per
    // unit where the sampler needs 1, so it is the largest allocation of the
    // whole run; this keeps the rest of the ensemble from being resident while
    // it is made.
    void release_all_but_plan_ids();

    // debugging methods
    // checks all plans are valid. 
    void check_all_plans_valid(
        MapParams const &map_params,
        std::string_view where
    );
};

PlanEnsemble get_plan_ensemble(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity);

std::unique_ptr<PlanEnsemble> get_plan_ensemble_ptr(
    MapParams const &map_params, SplittingSchedule const &splitting_schedule,
    int const num_regions, int const nsims, SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, Rcpp::IntegerMatrix const &region_sizes_mat,
    std::vector<RNGState> &rng_states, RcppThread::ThreadPool &pool, int const verbosity,
    InitialAugmentedSamples const &initial_augmented_samples = {});


// converts trees to compact edge list form 
// This only supports undirected trees 
std::vector<bool> vector_tree_to_edge_vector(
    GraphEdgeIndex const &edge_index,
    Tree const &tree
);

// Converts a graph edge index to an R num_edges by 2 matrix where
// row e (1-indexed e + 1) is the 1-indexed vertex pair (u, v), u < v,
// for edge id e
Rcpp::IntegerMatrix graph_edge_index_to_matrix(
    GraphEdgeIndex const &edge_index
);

// Returns the canonical edge list of a 0-indexed adjacency list as a
// num_edges by 2 1-indexed matrix. Row e is the edge with edge id e - 1.
Rcpp::IntegerMatrix graph_edge_list(Rcpp::List const &adj_list);


#endif
