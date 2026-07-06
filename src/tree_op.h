#ifndef TREE_OP_H
#define TREE_OP_H

#include "redist_types.h"
#include "smc_base.h"
#include <RcppArmadillo.h>
#include <limits>
#include <queue>
#include <stack>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;



/*
 * Generate a random vertex (integer) among unvisited vertices
 * `lower` is a lower bound (inclusive) on the index of the first unvisited element
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining, int &lower,
         RNGState &rng_state);

/*
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
int rnbor(const Graph &g, int vtx, RNGState &rng_state);

/*
 * Initialize empty tree structure on graph with `V` vertices
 */
// TESTED
Tree init_tree(int V);

/*
 * Clear a tree
 */
void clear_tree(Tree &tree);

// print a tree
void print_tree(Tree const &ust);

/*
 * Count population below each node in tree and get parent
 */
// TESTED
int tree_pop(Tree &ust, int vtx, const arma::uvec &pop, std::vector<int> &pop_below,
             std::vector<int> &parent);

/*
 * Just Count population below each node in tree
 */
// TESTED
void get_tree_pops_below(const Tree &ust, const int root, TreePopStack &stack,
                         const arma::uvec &pop, std::vector<int> &pop_below);


// Essentially a packed forest, aka a spanning forest stored as a list of edges
// in a bitpacked boolean array form  
class EdgeBitset {
  private:
    PlanEdgeBits edge_bits;

    static int word_index(EdgeID edge_id) {
        return static_cast<int>(edge_id >> 6);
    }

    static int bit_index(EdgeID edge_id) {
        return static_cast<int>(edge_id & 63);
    }

    static EdgeBitWord bit_mask(EdgeID edge_id) {
        return EdgeBitWord{1} << bit_index(edge_id);
    }

  public:
    explicit EdgeBitset(PlanEdgeBits &this_plan_edge_bits)
        : edge_bits(this_plan_edge_bits) {
            clear_all();
        }

    std::size_t size_words() const noexcept {
        return edge_bits.size();
    }

    bool empty() const {return edge_bits.empty();}

    void clear_all() {
        std::fill(edge_bits.begin(), edge_bits.end(), EdgeBitWord{0});
    }

    void copy(EdgeBitset const &other) {
        edge_bits.copy(other.edge_bits);
    }

    bool test_edge_id(EdgeID edge_id) const {
        if constexpr (perf_config::unnecessary_input_checks){
            int const word = word_index(edge_id);
            if (word < 0 || word >= static_cast<int>(edge_bits.size())) {
                std::ostringstream oss;
                Rcpp::Rcerr << "EdgeBitset::test_edge_id out of range. "
                    << "edge_id=" << edge_id
                    << ", word=" << word
                    << ", edge_bits.size()=" << edge_bits.size();
                throw Rcpp::exception("test_edge_id");
                throw std::runtime_error(oss.str());
            }
        }

        return (edge_bits[word_index(edge_id)] & bit_mask(edge_id)) != 0;
    }

    void set_edge_id(EdgeID edge_id) {
        if constexpr (perf_config::unnecessary_input_checks){
            int const word = word_index(edge_id);
            if (word < 0 || word >= static_cast<int>(edge_bits.size())) {
                std::ostringstream oss;
                Rcpp::Rcerr << "EdgeBitset::test_edge_id out of range. "
                    << "edge_id=" << edge_id
                    << ", word=" << word
                    << ", edge_bits.size()=" << edge_bits.size();
                throw Rcpp::exception("set_edge_id");
                throw std::runtime_error(oss.str());
            }
        }
        edge_bits[word_index(edge_id)] |= bit_mask(edge_id);
    }

    void clear_edge_id(EdgeID edge_id) {
        if constexpr (perf_config::unnecessary_input_checks){
            int const word = word_index(edge_id);
            if (word < 0 || word >= static_cast<int>(edge_bits.size())) {
                std::ostringstream oss;
                Rcpp::Rcerr << "EdgeBitset::test_edge_id out of range. "
                    << "edge_id=" << edge_id
                    << ", word=" << word
                    << ", edge_bits.size()=" << edge_bits.size();
                throw Rcpp::exception("clear_edge_id");
                throw std::runtime_error(oss.str());
            }
        }
        edge_bits[word_index(edge_id)] &= ~bit_mask(edge_id);
    }

    void set_edge(int v, int u, GraphEdgeIndex const &edge_index) {
        set_edge_id(edge_index.get_edge_id(v, u));
    }

    void clear_edge(int v, int u, GraphEdgeIndex const &edge_index) {
        clear_edge_id(edge_index.get_edge_id(v, u));
    }

    bool test_edge(int v, int u, GraphEdgeIndex const &edge_index) const {
        return test_edge_id(edge_index.get_edge_id(v, u));
    }

    Tree get_graph_tree(GraphEdgeIndex const &edge_index) const {
        Tree full_tree(edge_index.V);

        for (EdgeID edge_id = 0;
            edge_id < static_cast<EdgeID>(edge_index.edges.size());
            ++edge_id) {

            if (!test_edge_id(edge_id)) {
                continue;
            }

            auto const endpoints = edge_index.edges[edge_id];

            int const u = static_cast<int>(endpoints.first);
            int const v = static_cast<int>(endpoints.second);

            full_tree[u].push_back(v);
            full_tree[v].push_back(u);
        }

        return full_tree;
    }

    void print(GraphEdgeIndex const &edge_index) const {
        REprintf("Packed Forest with %zu edges!", edge_index.edges.size());
        for (EdgeID edge_id = 0;
            edge_id < static_cast<EdgeID>(edge_index.edges.size());
            ++edge_id) {

            auto const endpoints = edge_index.edges[edge_id];

            int const u = static_cast<int>(endpoints.first);
            int const v = static_cast<int>(endpoints.second);

            Rcpp::Rcout << "(" << u << ", " << v << "): "
                        << (test_edge_id(edge_id) ? "INCLUDED" : "NOT INCLUDED")
                        << "\n";
        }
    }

    void print_full_tree(GraphEdgeIndex const &edge_index) const {
        Tree full_tree(edge_index.V);

        for (EdgeID edge_id = 0;
            edge_id < static_cast<EdgeID>(edge_index.edges.size());
            ++edge_id) {

            if (!test_edge_id(edge_id)) {
                continue;
            }

            auto const endpoints = edge_index.edges[edge_id];

            int const u = static_cast<int>(endpoints.first);
            int const v = static_cast<int>(endpoints.second);

            full_tree[u].push_back(v);
            full_tree[v].push_back(u);
        }

        print_tree(full_tree);
    }

    // Clears every currently stored forest edge whose two endpoints are both
    // assigned to region_id.
    //
    // This is the packed-forest equivalent of clearing the old tree inside
    // a region before writing a newly sampled tree for that region.
    // NOTE: In the future for more than 2 region splits you could
    // just pass in a vector of region ids 
    void clear_region_tree(
        PlanVector const &region_ids,
        RegionID const region_id,
        GraphEdgeIndex const &edge_index
    ); 

    void clear_region_trees(
        PlanVector const &region_ids,
        RegionID const region_id1, RegionID const region_id2,
        GraphEdgeIndex const &edge_index
    ); 

    // Takes a vector tree ust, clears it and replace it with the current 
    // packed forest in the plan 
    void fill_vector_tree(
        GraphEdgeIndex const &edge_index,
        Tree &ust
    ) const;

};




/*
 * Assign `new_region_num_id` to all descendants of `root` in `ust`
 */
void assign_region_id_from_tree(const Tree &ust, PlanVector &region_ids, int root,
                                const int new_region_num_id,
                                CircularQueue<std::pair<int, int>> &vertex_queue);

void assign_region_id_and_forest_from_tree(const Tree &ust, PlanVector &region_ids,
                                           Tree &forest_graph, EdgeBitset &forest_edges,
                                           int root,
                                           const int new_region_id,
                                           const GraphEdgeIndex &graph_edge_index,
                                           CircularQueue<std::pair<int, int>> &vertex_queue);

/*
 * Assign `new_region_num_id` to all descendants of `root` in `ust` where
 * `ust` is a directed spanning tree 
 */
void assign_region_id_and_forest_from_tree_NEW(const Tree &ust, PlanVector &region_ids,
                                           EdgeBitset &forest_edges, int root,
                                           const int new_region_id,
                                           const GraphEdgeIndex &graph_edge_index,
                                           CircularQueue<std::pair<int, int>> &vertex_queue);

/*
 *  Erases an edge from a tree
 *
 *  Erases the directed edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`)
 *  from the tree `ust`. The directed edge here means we have `child_vertex` being one of
 *  the values in `ust[parent_vertex]`.
 *
 *
 *  @param ust A directed spanning tree passed by reference
 *  @param cut_edge An `EdgeCut` object representing the edge cut
 *
 *  @details Modifications
 *     - The edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`)
 *     is removed from `ust`
 *
 *
 */
void erase_tree_edge(Tree &ust, EdgeCut cut_edge);






#endif
