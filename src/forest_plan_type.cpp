#include <RcppArmadillo.h>
#include "forest_plan_type.h"

#include "scoring.h"
#include "splitting_schedule_types.h"
#include "tree_splitting.h"
#include "wilson.h"

// constructor for partial plan (more than 1 region)
ForestPlan::ForestPlan(int const ndists, int const num_regions, const std::vector<unsigned int> &pop,
                       PlanVector &this_plan_region_ids, RegionSizes &this_plan_region_sizes,
                       IntPlanAttribute &this_plan_region_pops,
                       IntPlanAttribute &this_plan_order_added, 
                       PlanEdgeBits &this_plan_forest_edge_bits,
                       USTSampler &ust_sampler,
                       RNGState &rng_state, const Rcpp::List &initial_forest_adj_list)
    : Plan(num_regions, pop, this_plan_region_ids, this_plan_region_sizes,
           this_plan_region_pops, this_plan_order_added, this_plan_forest_edge_bits) {

    if (initial_forest_adj_list.size() > 1) {
        throw Rcpp::exception("Input Forest list not supported right now\n");
    } else {
        // else just build a forest at random
        for (size_t region_id = 0; region_id < num_regions; region_id++) {
            auto result = draw_tree_on_region(ust_sampler, region_id,
                                              rng_state, 1000000);

            if (!result.first) {
                std::ostringstream oss;

                oss << "Failed to draw tree on region " << region_id
                    << " after 1000000 attempts.\n";
                oss << debug_string(true);

                throw std::runtime_error(oss.str());
            }
        }
    }

    if constexpr(perf_config::object_integrity_checking){
        check_forest_equality(
            ust_sampler.get_vertex_tree(),
            forest_edges.get_graph_tree(ust_sampler.map_params.graph_edge_index),
            ust_sampler.map_params.graph_edge_index,
            "IN Partial Forest Plan Constructor, checking forest_graph vs forest edges (through get_graph_tree)"
        );
    }

}

Tree ForestPlan::get_forest_adj() { throw std::runtime_error("get_forest_adj not supported right now!\n") ; }

// IT IS VERY IMPORTANT THAT FOR SMC split_region1_id is the id of the multidistrict
// The idea is any other split regions have not actually been updated yet 
void ForestPlan::update_vertex_and_plan_specific_info_from_cut(
    TreeSplitter const &tree_splitter, USTSampler &ust_sampler, EdgeCut const cut_edge,
    const int split_region1_id, const int split_region2_id, bool const add_region) {

    update_plan_ids_and_forest_from_cut(tree_splitter, 
        ust_sampler, cut_edge,
        split_region1_id, split_region2_id, add_region
    );

    return;
}

/*
 * @title Get Valid Adjacent Regions to Effective Tree Boundary Lengths Map
 *
 * Returns a unordered_map mapping valid pairs of adjacent regions to
 * the length of the effect tree boundary between them in the forest
 *
 * Finds all valid pairs of adjacent regions (meaning their sizes are valid
 * to merge) and returns a hash map mapping the pairs of valid adjacent
 * regions to the length of the effective tree boundary between them in `g`.
 * Recall this is the probability that we chose the actual split edge among
 * the tree made by merging each boundary edge
 *
 *
 * @param g A graph (adjacency list) passed by reference
 * @param plan A plan object
 * @param check_adj_to_regions A vector tracking whether or not we should
 * check for edges adjacent to vertices in a region of a particular size. For
 * example, `check_adj_to_regions[i] == true` means we will attempt to find
 * edges adjacent to any vertex in a region of size i. This vector is 1-indexed
 * meaning we don't need to subtract region size by 1 when accessing.
 * @param valid_merge_pairs A 2D `ndists+1` by `ndists+1` boolean matrix that
 * uses region size to check whether or not two regions are considered a valid
 * merge that can be counted in the map. For example `valid_merge_pairs[i][j]`
 * being true means that any regions where the sizes are (i,j) are considered
 * valid to merge. 2D matrix is 1 indexed (don't need to subtract region size)
 *
 * @details No modifications to inputs made
 * @return A hash map mapping (std::pair<int,int> to double) that maps a pair of
 * region ids in the form (smaller id, bigger id) to the effective tree boundary length
 * between them in `g`.
 */
std::vector<std::tuple<RegionID, RegionID, double>> compute_log_tree_eff_boundary_lens(
    const ForestPlan &plan, EdgeBitset const &forest_edges,
    PlanMultigraph &plan_multigraph,
    const SplittingSchedule &splitting_schedule, USTSampler &ust_sampler,
    TreeSplitter &edge_splitter, ScoringFunction const &scoring_function) {

    if constexpr (perf_config::object_integrity_checking){
        plan.check_forest_integrity(
            ust_sampler.map_params.graph_edge_index,
            "IN compute_log_tree_eff_boundary_lens, checking forest_edges integrity"
        );
    }

    // now make the output vector
    std::vector<std::tuple<RegionID, RegionID, double>> region_pairs_tuple_vec;
    region_pairs_tuple_vec.reserve(plan_multigraph.pair_map.num_hashed_pairs);

    for (auto const key_val_pair : plan_multigraph.pair_map.get_all_values()) {
        RegionID const region1_id = key_val_pair.first.first;
        RegionID const region2_id = key_val_pair.first.second;
        PairHashData const &pair_data = key_val_pair.second;    
        int const region1_size = plan.region_sizes[region1_id];
        int const region2_size = plan.region_sizes[region2_id];
        /*
        * Skip pairs that have already been marked as not needing their
        * effective boundary computed.
        */
        if (!pair_data.count_pair) {
            continue;
        }
        // skip if we don't check adjacent pairs for these sizes
        if (
            !splitting_schedule.check_adj_to_regions[region1_size] &&
            !splitting_schedule.check_adj_to_regions[region2_size]
        ) {
            continue;
        }

        int const merged_region_size = region1_size + region2_size;

        auto const cut_size_bounds = splitting_schedule
                .all_regions_min_and_max_possible_cut_sizes[
                    merged_region_size
                ];

        int const min_possible_cut_size = cut_size_bounds.first;
        int const max_possible_cut_size = cut_size_bounds.second;

        /*
        * If the regions are administratively adjacent, only within-county
        * map edges contribute. Otherwise across-county edges are allowed.
        *
        * When counties are off, all map edges are countable.
        */
        bool const count_edges_across = !plan_multigraph.counties_on || !pair_data.admin_adjacent;

       /*
        * The representative edge is stored in canonical region-ID order:
        *
        *   representative_vertex1 belongs to region1_id
        *   representative_vertex2 belongs to region2_id
        *
        * where region1_id < region2_id.
        */
        int const region1_root = pair_data.representative_vertex1;
        int const region2_root = pair_data.representative_vertex2;
        // now compute the log eff becoundary length
        double const log_eff_boundary_len = edge_splitter.get_log_eff_boundary_len_for_adj_region_pair(
                plan_multigraph.map_params,
                scoring_function,
                forest_edges,
                ust_sampler.pops_below_vertex,
                region1_root,
                region2_root,
                plan,
                min_possible_cut_size,
                max_possible_cut_size,
                splitting_schedule
                    .all_regions_smaller_cut_sizes_to_try[
                        merged_region_size
                    ],
                count_edges_across
            );

        region_pairs_tuple_vec.emplace_back(
            region1_id,
            region2_id,
            log_eff_boundary_len
        );
    }

    return region_pairs_tuple_vec;
}

std::vector<std::tuple<RegionID, RegionID, double>>
ForestPlan::get_valid_adj_regions_and_eff_log_boundary_lens(
    PlanMultigraph &plan_multigraph, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function, bool const is_final_split, USTSampler &ust_sampler,
    TreeSplitter &tree_splitter) const {
    // build the multigraph
    plan_multigraph.build_plan_multigraph(region_ids, num_regions);
    // remove invalid hierarchical merges if needed
    plan_multigraph.remove_invalid_hierarchical_merge_pairs(*this);
    // remove invalid size pairs
    plan_multigraph.remove_invalid_size_pairs(*this, splitting_schedule);
    // remove invalid hard constraint pairs
    plan_multigraph.remove_invalid_hard_constraint_pairs(*this, scoring_function,
                                                         is_final_split);

    // get pairs and log tree effective boundary
    auto region_pairs_tuple_vec = compute_log_tree_eff_boundary_lens(
        *this, forest_edges, plan_multigraph, splitting_schedule, ust_sampler, tree_splitter,
        scoring_function);

    return region_pairs_tuple_vec;
}

double ForestPlan::get_log_eff_boundary_len(PlanMultigraph &plan_multigraph,
                                            const SplittingSchedule &splitting_schedule,
                                            USTSampler &ust_sampler,
                                            TreeSplitter &tree_splitter,
                                            ScoringFunction const &scoring_function,
                                            const int region1_id, int const region2_id) const {

    if constexpr (perf_config::object_integrity_checking){
        check_forest_integrity(
            ust_sampler.map_params.graph_edge_index,
            "IN get_log_eff_boundary_len, checking forest_edges integrity"
        );
    }

    int const merged_region_size = region_sizes[region1_id] + region_sizes[region2_id];

    auto cut_size_bounds =
        splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    auto const search_result = plan_multigraph.pair_map.get_value(
            region1_id,
            region2_id
        );

    PairHashData const pair_data = search_result.second;
    int const region1_root = pair_data.representative_vertex1;
    int const region2_root = pair_data.representative_vertex2;

    // figure out if we are allowed to count across county boundaries
    bool count_edges_across;
    if (plan_multigraph.counties_on) {
        // only count across if not admin adjacent
        count_edges_across = !search_result.second.admin_adjacent;
    } else {
        count_edges_across = true;
    }

    double log_eff_boundary = tree_splitter.get_log_eff_boundary_len_for_adj_region_pair(
                        plan_multigraph.map_params, scoring_function, forest_edges,
                        ust_sampler.pops_below_vertex, region1_root,
                        region2_root, *this, min_possible_cut_size, max_possible_cut_size,
                        splitting_schedule
                            .all_regions_smaller_cut_sizes_to_try[merged_region_size],
                        count_edges_across);

    return log_eff_boundary;
}