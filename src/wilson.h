#ifndef WILSON_H
#define WILSON_H

#include "tree_op.h"
#include "base_plan_type.h"
#include "redist_types.h"
#include "tree_op.h"
#include "redist_constants.h"
#include <RcppArmadillo.h>


class Plan;
class ScoringFunction;
class TreeSplitter;


struct USTDrawResult {
    bool successful;
    int num_vertices;
    int root;
};

// Class for wrapping wilson code in 
class USTSampler {

  private:

    // private method which calls `sample_sub_ust` and assumes that visited and ignore have been properly set up
    // first entry is whether tree was successfully drawn
    // second is the root of the tree 
    std::pair<bool, int> draw_ust(double const lower, double const upper,
      RNGState &rng_state);

    std::pair<bool, EdgeCut>
    try_to_sample_splittable_tree(Plan const &plan, int const split_region1,
                                  int const split_region2, int const root,
                                  ScoringFunction const &scoring_function, RNGState &rng_state,
                                  TreeSplitter &tree_splitter, int const region_populations,
                                  int const region_size, bool const save_selection_prob);

  public:
    USTSampler(MapParams const &map_params, SplittingSchedule const &splitting_schedule)
        :
        ust(init_tree(map_params.V)),
        pops_below_vertex(map_params.V, 0),
          visited(map_params.V), ignore(map_params.V), stack(map_params.V + 1),
          county_tree(init_tree(map_params.num_counties)),
          county_stack(map_params.num_counties + 1),
          dummy_county_tree_queue(map_params.V),
          county_pop(map_params.num_counties, arma::fill::zeros),
          county_members(map_params.num_counties, std::vector<int>{}),
          c_visited(map_params.num_counties, true), cty_pop_below(map_params.num_counties, 0),
          vertex_queue(map_params.V), map_params(map_params),
          splitting_schedule(splitting_schedule) {
            // reserve the max capacity now 
            for (size_t v = 0; v < map_params.V; v++)
            {
              ust[v].reserve(map_params.g[v].size());
            }
             
          };

    Tree ust;
    std::vector<int> pops_below_vertex;
    std::vector<bool> visited, ignore;
    TreePopStack stack;
    Tree county_tree;
    TreePopStack county_stack;
    DummyTreeQueue dummy_county_tree_queue;
    arma::uvec county_pop;
    std::vector<std::vector<int>> county_members;
    std::vector<bool> c_visited;
    std::vector<int> cty_pop_below;
    std::vector<std::array<int, 3>> county_path;
    std::vector<int> path;
    CircularQueue<std::pair<int, int>> vertex_queue;
    MapParams const &map_params;
    SplittingSchedule const &splitting_schedule;

    // just used to draw a tree on a generic subgraph.
    bool draw_tree_on_subgraph(RNGState &rng_state, bool membership_vector);

    // Attempts to draw a tree on a region
    // defaults to map_params.lower * min_possible_cut_size and
    // map_params.upper * max_possible_cut_size as the bounds if 
    // use_custom_bounds = false
    USTDrawResult attempt_to_draw_tree_on_region(RNGState &rng_state, Plan const &plan,
                                        const int region_to_draw_tree_on, 
                                        bool const use_custom_bounds = false,
                                        double const custom_sample_sub_ust_lower = 0,
                                        double const custom_sample_sub_ust_upper = 0);

    // Attempts to draw a tree on a region formed by merging the two regions
    USTDrawResult attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
                                               const int region1_to_draw_tree_on,
                                               const int region2_to_draw_tree_on, 
                                                bool const use_custom_bounds = false,
                                                double const custom_sample_sub_ust_lower = 0,
                                                double const custom_sample_sub_ust_upper = 0);

    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_split(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const region_to_split,
        int const new_region_id, bool const save_selection_prob);


    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_mergesplit(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const merge_region1,
        int const merge_region2, bool const save_selection_prob);

    // checks that all the vertices in the tree are valid and 
    // its a directed tree 
    void check_tree_integrity(
      Tree const &a_ust,
      std::string_view where,
      int root,
      int expected_tree_vertices,
      bool check_vertex_count
    );
};




#endif