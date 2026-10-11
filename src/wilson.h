#ifndef WILSON_H
#define WILSON_H

#include <string_view>

#include "tree_op.h"
#include "advanced_types.h"

class Plan;
class ScoringFunction;
class TreeSplitter;
class SplittingSchedule;
struct CutSizeOptions;


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


// The number of admin units the Wilson scratch objects must be able to hold.
//
// When plans are guaranteed hierarchically valid the admin units are exactly
// the counties. Otherwise a county's active portion can break into several
// disconnected pieces, each of which becomes its own fake county, so in the
// worst case every vertex is its own admin unit.
inline int admin_unit_capacity(MapParams const &map_params) {
    return map_params.plans_may_be_non_hierarchical
        ? map_params.V
        : map_params.num_counties;
}


class WilsonMultiGraphScratch {
  public:
    explicit WilsonMultiGraphScratch(
        int const admin_capacity
    )
        : admin_capacity(admin_capacity),
          num_admin_units(0),
          total_pop(0),
          smallest_county_seen(-1),
          c_remaining(0),
          county_stack(admin_capacity + 1),
          county_pop(admin_capacity, 0),
          c_visited(admin_capacity, true),
          admin_ignore(admin_capacity, true),
          cty_pop_below(admin_capacity, 0),
          next_admin_edge(
              admin_capacity,
              AdminEdge{-1, -1, -1}
          ),
          admin_roots(admin_capacity, -1),
          deterministic_counties() {
            deterministic_counties.reserve(admin_capacity);
          }

    // Allocated length of every per-admin-unit vector below. Never changes,
    // and matches the size of the multigraph the admin walk runs on.
    int const admin_capacity;
    // Number of admin units actually live in the current draw. Always
    // <= admin_capacity, and equal to it when plans are hierarchically valid.
    // Reset by `prep_fresh_ust_call` on every draw.
    //
    // THIS, not admin_capacity, bounds every loop over admin units. Entries
    // at or past this index hold stale values from earlier draws.
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
    // vector of the counties where a deterministic tree was filled in
    std::vector<int> deterministic_counties;
};

// Class for wrapping wilson code in 
class USTSampler {

  private:
    std::vector<bool> visited;
    std::vector<bool> ignore; // TODO make visited private and create one for tree splitter

    // This takes a function `is_active` which given a vertex v
    // returns true if it should be included and false if not
    // ie false value means ignore
    // Using this it preps all the inputs for a fresh ust call
    //
    // `counties` is the 1-indexed admin unit label of every map vertex and
    // `num_admin_units` is how many of those units exist. They are passed in
    // rather than read off `map_params` so a draw can run over fake counties.
    // This is the only writer of `mg_scratch.num_admin_units`.
    template <typename IsActive>
    int prep_fresh_ust_call(
        IsActive const &is_active,
        std::vector<unsigned int> const &counties,
        int const num_admin_units
    );


    // Refines the counties into fake counties on the active subgraph, where
    // each fake county is one connected piece of (active vertices) intersect
    // (a county). Needed when a plan is not hierarchically valid, since then
    // a county's active portion can be disconnected and the within county
    // Wilson walk could never reach the far pieces.
    //
    // Fills `fake_counties` and `fake_cg` and returns how many fake counties
    // there are. Does not touch `mg_scratch.num_admin_units`; the caller
    // hands the returned count to `draw_fresh_ust` like any other partition.
    template <typename IsActive>
    int build_fake_counties(IsActive const &is_active);


    // Draws a completely fresh ust over the administrative partition given by
    // (`counties`, `admin_multigraph`, `num_admin_units`), which must all
    // describe the same partition.
    template <bool SkipUnsplittableTrees, typename IsActive>
    USTDrawResult draw_fresh_ust(
        double const lower,
        double const upper,
        RNGState &rng_state,
        IsActive const &is_active,
        std::vector<unsigned int> const &counties,
        Multigraph const &admin_multigraph,
        int const num_admin_units
    );

    // Finds and 
    std::pair<bool, EdgeCut>
    try_to_find_and_erase_splittable_edge(Plan const &plan, int const split_region1,
                                  int const split_region2, int const root,
                                  ScoringFunction const &scoring_function, RNGState &rng_state,
                                  TreeSplitter &tree_splitter, int const region_populations,
                                  int const region_size, CutSizeOptions const &cut_options,
                                  bool const save_selection_prob);

  public:
    USTSampler(MapParams const &map_params, SplittingSchedule const &splitting_schedule)
        :
          visited(map_params.V), ignore(map_params.V),
          assume_hierarchical(!map_params.plans_may_be_non_hierarchical),
          // ust(map_params.map_graph.get_flat_empty_tree()),
          ust(init_tree(map_params.V)),
          // wilson_submap(map_params.map_graph),
          pops_below_vertex(map_params.V, 0),
          stack(map_params.V + 1),
          county_tree(init_tree(admin_unit_capacity(map_params))),
          g_scratch(map_params.V),
          mg_scratch(admin_unit_capacity(map_params)),
          vertex_queue(map_params.V), map_params(map_params),
          splitting_schedule(splitting_schedule),
          // Only allocated when plans can be non-hierarchical; empty otherwise
          // so the hierarchical path pays nothing for them.
          fake_counties(
              map_params.plans_may_be_non_hierarchical ? map_params.V : 0, 0
          ),
          fake_cg(
              map_params.plans_may_be_non_hierarchical ? map_params.V : 0
          ) {
            // reserve the max capacity now 
            for (size_t v = 0; v < map_params.V; v++)
            {
              ust[v].reserve(map_params.g[v].size());
            }
             
          };

    bool const assume_hierarchical; // Whether or not we can assume the inputted plans will always be hierarchical
    // FlatGraph ust;
    Tree ust;
    // FlatGraph wilson_submap; // subgraph of g restricted to only the vertices we care about 
    std::vector<int> pops_below_vertex;
    TreePopStack stack; // graph
    Tree county_tree; // county 
    WilsonGraphScratch g_scratch;
    WilsonMultiGraphScratch mg_scratch;
    CircularQueue<std::pair<int, int>> vertex_queue; // not used in sample ust
    MapParams const &map_params;
    SplittingSchedule const &splitting_schedule;
    // These are only needed if we're not guaranteed hierarchical plans, and
    // are left empty otherwise. `fake_counties` gives every map vertex its
    // fake county label and `fake_cg` is the multigraph over those units,
    // sized to `mg_scratch.admin_capacity` with only the first
    // `mg_scratch.num_admin_units` rows live in any given draw.
    //
    // Fake county labels are 1-indexed to match `map_params.counties`, so a
    // label of 0 marks a vertex outside the active subgraph.
    std::vector<unsigned int> fake_counties;
    Multigraph fake_cg;


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
    bool fill_in_skipped_subtrees(
      EdgeBitset &packed_forest_edges,
      RNGState &rng_state, 
      int const max_tries_multiple = 1000
    );

    // `cut_options` are the allowed ways to split `region_to_split`, from the
    // splitting schedule
    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_split(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const region_to_split,
        int const new_region_id, CutSizeOptions const &cut_options,
        bool const save_selection_prob);


    // `cut_options` are the allowed ways to split the merged region, from the
    // splitting schedule
    std::pair<bool, EdgeCut> attempt_to_find_valid_tree_mergesplit(
        RNGState &rng_state, ScoringFunction const &scoring_function,
        TreeSplitter &tree_splitter, Plan const &plan, int const merge_region1,
        int const merge_region2, CutSizeOptions const &cut_options,
        bool const save_selection_prob);

    // Tree get_vertex_tree ()const {return ust.to_vertex_graph();}
    Tree get_vertex_tree ()const {return ust;}

    /*
     * The same tree with every edge stored at both endpoints.
     *
     * `ust` is directed, holding each edge only under its parent, whereas a
     * forest recovered from the packed edge bitset is undirected. Comparing
     * the two directly fails at the first leaf, so anything checking the
     * sampler's tree against a packed forest wants this instead.
     */
    Tree get_undirected_vertex_tree() const {
        Tree undirected(ust.size());

        for (std::size_t v = 0; v < ust.size(); ++v) {
            for (int const child : ust[v]) {
                undirected[v].push_back(static_cast<int>(child));
                undirected[child].push_back(static_cast<int>(v));
            }
        }

        return undirected;
    }

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