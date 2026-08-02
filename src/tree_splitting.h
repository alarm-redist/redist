#pragma once
#ifndef TREE_SPLITTING_H
#define TREE_SPLITTING_H

#include <RcppArmadillo.h>
#include <queue>
#include <utility>
#include <vector>

#include "advanced_types.h"
#include "tree_op.h"

class ScoringFunction;
class Plan;
class RNGState;


// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut
class TreeSplitter {

  public:
    // Default Constructor
    TreeSplitter(FlatGraph const &map_graph) :
    forest_graph(map_graph.get_flat_empty_tree()), 
    visited(map_graph.size()),
    no_valid_edges_vertices(map_graph.size()),
    stack(map_graph.size()+1)
    // balanced_edge_cuts()
    {
        // balanced_edge_cuts.reserve(
        //     static_cast<int>(std::ceil(.1 * map_graph.size()))
        // );
    };
    TreeSplitter(int const V) : 
        forest_graph(), 
        visited(0),
        no_valid_edges_vertices(V),
        stack(V + 1) { 
        // balanced_edge_cuts.reserve(
        //     static_cast<int>(std::ceil(.1 * V))
        // );
    };


    virtual ~TreeSplitter() = default;

    FlatGraph forest_graph; // used for computing get_log_retroactive_splitting_prob_for_joined_flattree
    std::vector<bool> visited; // used in retroactive prob so not needed for naive k
    mutable std::vector<bool> no_valid_edges_vertices; // used in finding balanced edge cuts
    mutable TreePopStack stack; // used in splitting so needed for all 
    // mutable std::vector<EdgeCut> balanced_edge_cuts;

    // Returns a vector of all the valid edges in the tree
    std::vector<EdgeCut> get_all_valid_pop_edge_cuts_in_directed_tree(
        MapParams const &map_params, Tree const &ust, // FlatGraph const &ust, 
        const int root, 
        std::vector<int> &pops_below_vertex, 
        int const region_population, int const region_size, const int min_potential_cut_size,
        const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try) const;

    // Takes a spanning tree and returns the edge to cut if successful
    virtual std::pair<bool, EdgeCut> attempt_to_find_edge_to_cut(
        const MapParams &map_params, ScoringFunction const &scoring_function,
        RNGState &rng_state, Plan const &plan, int const split_region1, int const split_region2,
        Tree const &ust, // FlatGraph const &ust, 
        const int root, 
        std::vector<int> &pops_below_vertex, 
        int const region_population, int const region_size, const int min_potential_cut_size,
        const int max_potential_cut_size, std::vector<int> const &smaller_cut_sizes_to_try,
        bool save_selection_prob = false);

    // Get probability a specific edge was cut in the tree made by joining
    // the trees in the two regions where the forest is stored as a packed forest
    // This is called by linking edge plans since its not worth it to copy 
    // over to a vertex forest 
    virtual double get_log_retroactive_splitting_prob_for_joined_packed_tree(
        MapParams const &map_params, ScoringFunction const &scoring_function,
        EdgeBitset const &forest_edges, 
        std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
        Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try);

    // Get probability a specific edge was cut in the tree made by joining
    // the trees in the two regions where the forest is stored as a vertex of vertex forest
    // This is primarily called by forest space plans since its cheaper to copy to vertex
    // forest once for the indexing gains 
    virtual double get_log_retroactive_splitting_prob_for_joined_flattree(
        MapParams const &map_params, ScoringFunction const &scoring_function,
        FlatGraph const &forest_graph, 
        std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
        Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try);


    virtual double get_log_retroactive_splitting_prob_from_valid_pop_cut_list(
        std::vector<EdgeCut> &valid_edges, EdgeCut const actual_cut_edge
    );

    // Takes a vector of valid edge cuts and returns the log probability
    // the one an index idx would have been chosen
    virtual double get_log_selection_prob(std::vector<EdgeCut> &valid_edges, int idx) const;

    // used to update the k parameter for top k splitter
    virtual void update_single_int_param(int int_param) {
        throw Rcpp::exception("Update single int param not implemented!\n");
    };

    // used to get the k parameter for top k splitter
    virtual int get_single_int_param() const {
        throw Rcpp::exception("Update single int param not implemented!\n");
        return -1;
    };

    // returns edge cut and log probability it was chosen
    // for tree splitters where the selection probability is solely a function of the 
    // balanced tree cuts. It can't depend on anything else
    // For that just make a custom version of 
    // get_log_retroactive_splitting_prob_for_joined_flattree
    virtual std::pair<bool, EdgeCut> select_edge_to_cut(
        RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
        bool save_selection_prob
    ) const;

    virtual double compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const {
        throw Rcpp::exception("Not implemented for this class!");
    };
};

// Splitting method that just tries to pick one of the top k edges unif
class NaiveTopKSplitter : public TreeSplitter {

  public:
    // Constructor for NaiveTopKSplitter
    NaiveTopKSplitter(int const V, int k_param) 
    : TreeSplitter(V), k_param(k_param) {}

    // Attributes specific to NaiveTopKSplitter
    int k_param; // Top k valuex

    // how to update the k param
    void update_single_int_param(int int_param) override;

    int get_single_int_param() const override { return k_param; };

    std::pair<bool, EdgeCut> select_edge_to_cut(RNGState &rng_state,
                                                std::vector<EdgeCut> &valid_edges,
                                                bool save_selection_prob) const override;

    double get_log_selection_prob(std::vector<EdgeCut> &valid_edges, int idx) const override {
        throw Rcpp::exception("No log selection prob implemented for naive k!\n");
    }
};

// Splitting method that just tries to pick one of the top k edges unif
class UniformValidSplitter : public TreeSplitter {

  public:
    // implementation of the pure virtual function
    UniformValidSplitter(FlatGraph const &map_graph) : 
    TreeSplitter(map_graph) {};

    std::pair<bool, EdgeCut> select_edge_to_cut(RNGState &rng_state,
                                                std::vector<EdgeCut> &valid_edges,
                                                bool save_selection_prob) const override;

    // since uniform log prob is just -log(# of candidates)
    double get_log_selection_prob(std::vector<EdgeCut> &valid_edges, int idx) const override {
        return -std::log(valid_edges.size());
    }
};

// Splitting method that picks edge w/ prob ∝ exp(-alpha*bigger dev)
class ExpoWeightedSplitter : public TreeSplitter {

  public:
    ExpoWeightedSplitter(FlatGraph const &map_graph,
        double alpha, double target)
        : TreeSplitter(map_graph), alpha(alpha), target(target) {
        if (alpha < 0.0)
            throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;
    double target;

    double compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const override;
};

class ExpoWeightedSmallerDevSplitter : public TreeSplitter {

  public:
    ExpoWeightedSmallerDevSplitter(FlatGraph const &map_graph,
        double alpha, double target)
        : TreeSplitter(map_graph), alpha(alpha), target(target) {
        if (alpha < 0.0)
            throw Rcpp::exception("Alpha must be greater than zero!");
    }

    double alpha;
    double target;

    virtual double compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const override;
};

// Splitting method that just tries to pick one of the top k edges unif
class PopTemperSplitter : public TreeSplitter {

  public:
    // implementation of the pure virtual function
    PopTemperSplitter(FlatGraph const &map_graph,
        double const target, double const pop_temper, int const ndists)
        : TreeSplitter(map_graph), target(target), pop_temper(pop_temper), ndists(ndists) {};

    double const target;
    double const pop_temper;
    int const ndists;

    // since uniform log prob is just -log(# of candidates)
    double compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const override;
};

class ExperimentalSplitter : public TreeSplitter {

  public:
    ExperimentalSplitter(FlatGraph const &map_graph,
        double epsilon, double target)
        : TreeSplitter(map_graph), epsilon(epsilon), target(target) {
        if (epsilon < 0.0)
            throw Rcpp::exception("Epsilon must be greater than zero!");
    }

    double epsilon;
    double target;

    std::pair<bool, EdgeCut> select_edge_to_cut(RNGState &rng_state,
                                                std::vector<EdgeCut> &valid_edges,
                                                bool save_selection_prob) const override;

    double get_log_selection_prob(std::vector<EdgeCut> &valid_edges, int idx) const override;
};

class ConstraintSplitter : public TreeSplitter {

  private:
    std::vector<RegionID> underlying_plans_vec;
    std::vector<RegionID> underlying_sizes_vec;
    std::vector<int> underlying_pops_vec;

  public:
    ConstraintSplitter(FlatGraph const &map_graph, int const V, int const ndists)
        : TreeSplitter(map_graph), underlying_plans_vec(V, 0), underlying_sizes_vec(ndists, 0),
          underlying_pops_vec(ndists, 0), V(V), ndists(ndists),
          region_ids(underlying_plans_vec, 0, V), region_sizes(underlying_sizes_vec, 0, ndists),
          region_pops(underlying_pops_vec, 0, ndists), vertex_queue(V), dummy_forest(V) {
        if (ndists < 0.0)
            throw Rcpp::exception("Ndists must be greater than zero!");
    }

    int const V;
    int const ndists;
    PlanVector region_ids;
    RegionSizes region_sizes;
    IntPlanAttribute region_pops;
    CircularQueue<std::pair<int, int>> vertex_queue;
    Tree dummy_forest;

    // Takes a spanning tree and returns the edge to cut if successful
    std::pair<bool, EdgeCut> attempt_to_find_edge_to_cut(
        const MapParams &map_params, ScoringFunction const &scoring_function,
        RNGState &rng_state, Plan const &plan, int const split_region1, int const split_region2,
        Tree const &ust, // FlatGraph const &ust, 
        const int root, 
        std::vector<int> &pops_below_vertex, 
        int const region_population, int const region_size, const int min_potential_cut_size,
        const int max_potential_cut_size, std::vector<int> const &smaller_cut_sizes_to_try,
        bool save_selection_prob = false) override;

    double get_log_retroactive_splitting_prob_for_joined_packed_tree(
        MapParams const &map_params, ScoringFunction const &scoring_function,
        EdgeBitset const &forest_edges,
        std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
        Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try) override;


    double get_log_retroactive_splitting_prob_for_joined_flattree(
        MapParams const &map_params, ScoringFunction const &scoring_function,
        FlatGraph const &forest_graph, 
        std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
        Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
        std::vector<int> const &smaller_cut_sizes_to_try) override;

    
    double custom_get_log_retroactive_splitting_prob_from_valid_pop_cut_list(
            std::vector<EdgeCut> &valid_edges, EdgeCut const actual_cut_edge,
            MapParams const &map_params, ScoringFunction const &scoring_function,
            Plan const &plan
    );

};



#endif
