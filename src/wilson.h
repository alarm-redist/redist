#ifndef WILSON_H
#define WILSON_H

#include "tree_op.h"
#include "advanced_types.h"

class Plan;
class ScoringFunction;
class TreeSplitter;
class SplittingSchedule;


struct USTDrawResult {
    bool successful;
    int num_vertices;
    int root;
};

// used for timing manual code 
struct WilsonTimes {
  double input_prep_time = 0.0;
  double sub_ust_call_time = 0.0;
};


// Simple struct to hold sratch objects 
// used by Wilson code on  graphs 
class WilsonGraphScratch {

  public:
    WilsonGraphScratch(int const V)
        : remaining(0),
          smallest_v_seen(-1),
          dummy_county_tree_queue(V + 1),
          next_vertex(V, -1) {}


    int remaining; // The number of vertices remaining 
    int smallest_v_seen; // The 
    DummyTreeQueue dummy_county_tree_queue;
    std::vector<int> next_vertex;

};


class WilsonMultiGraphScratch {
  public:
    explicit WilsonMultiGraphScratch(
        int const num_admin_units
    )
        : num_admin_units(num_admin_units),
          total_pop(0),
          smallest_county_seen(-1),
          c_remaining(0),
          county_stack(num_admin_units + 1),
          county_pop(num_admin_units, 0),
          c_visited(num_admin_units, true),
          admin_ignore(num_admin_units, true),
          cty_pop_below(num_admin_units, 0),
          next_admin_edge(
              num_admin_units,
              AdminEdge{-1, -1, -1}
          ),
          admin_roots(num_admin_units, -1) {}

    int num_admin_units;
    int total_pop; // tracks total population
    int smallest_county_seen;
    int c_remaining; // tracks counties left to split
    TreePopStack county_stack; // county 
    std::vector<unsigned int> county_pop; // county
    std::vector<bool> c_visited; // county
    std::vector<bool> admin_ignore;
    std::vector<int> cty_pop_below; // county

    /*
     * Multigraph edges selected along the active walk.
     */
    std::vector<AdminEdge> next_admin_edge;

    // Maps each admin unit to a vertex in said admin unit
    // This is used for filling in a tree after acceptance
    std::vector<int> admin_roots;
};

// Class for wrapping wilson code in 
class USTSampler {

  private:
    std::vector<bool> visited;

    // This takes a function `is_active` which given a vertex v
    // returns true if it should be included and false if not 
    // ie false value means ignore
    // Using this it preps all the inputs for a fresh ust call
    template <typename IsActive>
    int prep_fresh_ust_call(
        IsActive const &is_active
    );


    // Draws a completely fresh ust. 
    template <bool SkipUnsplittableTrees, typename IsActive>
    USTDrawResult draw_fresh_ust(
        double const lower,
        double const upper,
        RNGState &rng_state,
        IsActive const &is_active
    );      

    // Finds and 
    std::pair<bool, EdgeCut>
    try_to_find_and_erase_splittable_edge(Plan const &plan, int const split_region1,
                                  int const split_region2, int const root,
                                  ScoringFunction const &scoring_function, RNGState &rng_state,
                                  TreeSplitter &tree_splitter, int const region_populations,
                                  int const region_size, bool const save_selection_prob);

  public:
    USTSampler(MapParams const &map_params, SplittingSchedule const &splitting_schedule)
        :
        // ust(map_params.map_graph.get_flat_empty_tree()),
        ust(init_tree(map_params.V)),
        // wilson_submap(map_params.map_graph),
        pops_below_vertex(map_params.V, 0),
          visited(map_params.V), ignore(map_params.V), stack(map_params.V + 1),
          county_tree(init_tree(map_params.num_counties)),
          g_scratch(map_params.V),
          mg_scratch(map_params.num_counties),
          vertex_queue(map_params.V), map_params(map_params),
          splitting_schedule(splitting_schedule) {
            // reserve the max capacity now 
            for (size_t v = 0; v < map_params.V; v++)
            {
              ust[v].reserve(map_params.g[v].size());
            }
             
          };

    // FlatGraph ust;
    Tree ust;
    // FlatGraph wilson_submap; // subgraph of g restricted to only the vertices we care about 
    std::vector<int> pops_below_vertex;
    std::vector<bool> ignore; // TODO make visited private and create one for tree splitter
    TreePopStack stack; // graph
    Tree county_tree; // county 
    WilsonGraphScratch g_scratch;
    WilsonMultiGraphScratch mg_scratch;
    CircularQueue<std::pair<int, int>> vertex_queue; // not used in sample ust
    MapParams const &map_params;
    SplittingSchedule const &splitting_schedule;

    // just used to draw a tree on a generic subgraph.
    // has toggable option to turn on or off skipping drawing 
    // a tree when an admin unit can't be split
    // Skipping is usually much faster but makes it impossible to test the
    // code as now not all trees have an equal probability of being sampled
    std::pair<bool, int> draw_tree_on_subgraph(
      RNGState &rng_state, std::vector<bool> const &vertices_to_ignore,
      bool const skip_unsplittable_subtrees, 
      const double lower, const double upper,
      WilsonTimes &wilson_times
    );

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


    // Fills in hierarchical subtrees that were skipped 
    // Will throw an error if for any subunit no trees can be filled after 
    // max_tries 
    void fill_in_skipped_subtrees(RNGState &rng_state, int max_tries = 1000);




    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_split(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const region_to_split,
        int const new_region_id, bool const save_selection_prob);


    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_mergesplit(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const merge_region1,
        int const merge_region2, bool const save_selection_prob);

    // Tree get_vertex_tree ()const {return ust.to_vertex_graph();}
    Tree get_vertex_tree ()const {return ust;}

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