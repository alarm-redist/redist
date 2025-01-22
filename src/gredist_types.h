#pragma once
#ifndef GREDIST_TYPES_H
#define GREDIST_TYPES_H


#define PRINT_LN Rcout << __func__ << "(), " << __FILE__ << ":" << __LINE__ << "\n";

#include <vector>
#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]


typedef std::vector<std::vector<int>> Tree;
typedef std::vector<std::vector<int>> Graph;
typedef std::vector<std::vector<std::vector<int>>> Multigraph;

/*
 * Initialize empty multigraph structure on graph with `V` vertices
 */
// TESTED
Multigraph init_multigraph(int V);


/*
 * Make a county graph from a precinct graph and list of counties
 * County graph is list of list of 3: <cty of nbor, index of vtx, index of nbor>
 */
// TESTED
Multigraph county_graph(const Graph &g, const arma::uvec &counties);


/*
 * Convert R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const Rcpp::List &l);

// Essentially just a useful container for map parameters 

class MapParams {
    public:
    // Constructor 
    MapParams(Rcpp::List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int ndists, double lower, double target, double upper);

    Graph g; // The graph as undirected adjacency list 
    arma::uvec counties; // county labels
    Multigraph cg; // county multigraph
    arma::uvec pop; // population of each vertex
    int V; // Number of vertices in the graph
    int ndists; // The number of districts a final plan should have
    double lower; // lower bound on district population
    double target; // target district population
    double upper; // upper bound on district population

};

// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut 
class EdgeCut {

public:
    // Default Constructor 
    EdgeCut()
        : tree_root(0), 
        cut_vertex(0), 
        cut_vertex_parent(0), 
        cut_below_region_size(0), 
        cut_below_pop(0), 
        cut_above_region_size(0), 
        cut_above_pop(0) {}
    
    // Constructor
    EdgeCut(const int tree_root, 
            const int cut_vertex, const int cut_vertex_parent,
            const int cut_below_region_size, const int cut_below_pop,
            const int cut_above_region_size, const int cut_above_pop)
        : tree_root(tree_root), 
          cut_vertex(cut_vertex), 
          cut_vertex_parent(cut_vertex_parent), 
          cut_below_region_size(cut_below_region_size), 
          cut_below_pop(cut_below_pop), 
          cut_above_region_size(cut_above_region_size), 
          cut_above_pop(cut_above_pop) {}
    
    // Attributes
    int tree_root; // The root of the tree
    int cut_vertex; // The vertex where we are cutting below it
    int cut_vertex_parent; // The parent of `cut_vertex` so we are cutting `(cut_vertex_parent, cut_vertex)`
    int cut_below_region_size; // The size of the region below made by cutting 
    int cut_below_pop; // The population of the region below made by cutting 
    int cut_above_region_size; // The size of the region above made by cutting 
    int cut_above_pop; // The population of the region above made by cutting 

    // Gets the information on the two regions formed from an edge cut by reference
    void get_split_regions_info(
        int &split_region1_tree_root, int &split_region1_dval, int &split_region1_pop,
        int &split_region2_tree_root, int &split_region2_dval, int &split_region2_pop
    );

    // Gets the signed (not absolute value) deviation of the two regions from the targets
    // first entry is below and second is above
    std::array<double, 2> compute_signed_pop_deviances(double target);

    // returns absolute population deviation
    std::array<double, 2> compute_abs_pop_deviances(double target);

    // Equality operator
    bool operator==(const EdgeCut& other) const {
        return tree_root == other.tree_root &&
               cut_vertex == other.cut_vertex &&
               cut_vertex_parent == other.cut_vertex_parent &&
               cut_below_region_size == other.cut_below_region_size &&
               cut_below_pop == other.cut_below_pop &&
               cut_above_region_size == other.cut_above_region_size &&
               cut_above_pop == other.cut_above_pop;
    }

    // Not-equal operator
    bool operator!=(const EdgeCut& other) const {
        return !(*this == other);
    }

    // Less-than operator
    bool operator<(const EdgeCut& other) const {
        return cut_vertex < other.cut_vertex;
    }

};


// enum for various methods of splitting a plan
enum class SplittingMethodType : unsigned char
{
    NaiveTopK, // picks 1 of top k edges even if invalid
    UnifValid, // picks uniform valid edge at random 
    ExpBiggerAbsDev // propto exp(-alpha*bigger abs dev of pair)
};

// loads a splitting type enum from a control string
SplittingMethodType get_splitting_type(std::string const &splitting_type_str);

// Get convinient string representation
std::string splitting_method_to_str(SplittingMethodType splitting_method);

enum class SplitRegionSizeType : unsigned char
{
    DistrictOnly,
    AnyValidSize,
    CustomSizes
};

// load from control spring 
SplitRegionSizeType get_splitting_size_regime(std::string const &splitting_size_regime_str);


class SpanningTree {
public:
    Tree ust;
    std::vector<int> parents, pop_below;
    std::vector<bool> visited, ignore;
    int V, tree_size;
};


#endif
