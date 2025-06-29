#pragma once
#ifndef GREDIST_TYPES_H
#define GREDIST_TYPES_H


#define PRINT_LN Rcout << __func__ << "(), " << __FILE__ << ":" << __LINE__ << "\n";

#include <vector>
#include <queue>
#include <cstdint>
#include <stdint.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

typedef uint16_t VertexID; // type for trees to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_VERTICES = static_cast<unsigned int>(
    std::numeric_limits<VertexID>::max()
); 

typedef std::uint_least8_t RegionID; // type for plan vectors to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_DISTRICTS = static_cast<unsigned int>(
    std::numeric_limits<RegionID>::max()
); 

typedef std::uint_least16_t CountyID; // type for county vectors to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_COUNTIES = static_cast<unsigned int>(
    std::numeric_limits<CountyID>::max()
); 

typedef std::uint_least16_t CountyRegion; // type for region_intersect_county lookup
constexpr unsigned int MAX_SUPPORTED_COUNTYREGION_VALUE = static_cast<unsigned int>(
    std::numeric_limits<CountyRegion>::max()
); 

//
typedef std::vector<RegionID> AllPlansVector;
typedef std::vector<RegionID> AllRegionSizesVector;
typedef std::vector<RegionID> RegionSizeVector;
typedef std::vector<std::vector<VertexID>> VertexGraph; // graphs on vertices
typedef std::vector<std::unordered_map<int, int>> RegionMultigraphCount;
typedef std::tuple<CountyRegion, RegionID, CountyID> CountyComponentVertex;
typedef std::vector<std::vector<CountyComponentVertex>> CountyComponentGraph;

// old types
typedef std::vector<std::vector<int>> Tree;
typedef std::vector<std::vector<int>> Graph;
typedef std::vector<std::vector<std::array<int, 3>>> Multigraph;



// This is a template class for storing plan attributes 
// in a way to minimize memory 
template <typename T> class PlanAttribute {
    private:
        int const offset_start; // starting index for this subset of vector 
        int const offset_end; // end index for this subset of vector 
        std::vector<T> &long_vec; // ref to underlying long vector 
    
    public:
        PlanAttribute(std::vector<T> &underyling_long_vec, int offset_start, int offset_end):
        offset_start(offset_start), offset_end(offset_end), long_vec(underyling_long_vec)
        {};

        // methods for accessing 
        // Const version for read-only access
        const T operator[](int index) const {
            return long_vec[offset_start + index];
        };
        // non constant for modification 
        T& operator [](int index) {
            return long_vec[offset_start + index];
        };

        // Non-const iterator accessors
        auto begin() { return long_vec.begin() + offset_start; }
        auto end() { return long_vec.begin() + offset_end; }

        // Const iterator accessors
        auto begin() const { return long_vec.begin() + offset_start; }
        auto end() const { return long_vec.begin() + offset_end; }

        // This copies one plans data from another 
        void copy(PlanAttribute const &other_attr){
            std::copy(
                other_attr.begin(),
                other_attr.end(),
                begin()
            );
        }

        std::size_t size() const noexcept{return offset_end - offset_start;};

};

typedef PlanAttribute<RegionID> PlanVector;
typedef PlanAttribute<RegionID> RegionSizes;
typedef PlanAttribute<int> IntPlanAttribute;


// uniquely maps pairs (x,y) of the form 
// 0 <= x < y < container_size
// to a value in (0, container_size choose 2)
// Currently used in the following scenarios
//  - Pairs of regions (region1_id, region2_id) with 0 <= region1_id < region2_id < num_regions
inline int index_from_ordered_pair(int x, int y, int container_size) {
    return ( (x * (2 * container_size - x - 1)) / 2 ) + (y - x - 1);
}

// indexing (i,j) in a matrix 
// assuming the length of a row is row_length
// i is the row index
// j is the column index so j < row_length
inline int mat_index_from_pair(int i, int j, int row_length){
    return i * row_length + j;
}


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

/*
 * Create a forest where each tree is a spanning tree on a county along
 with a vector of roots. Lets you traverse all counties in O(V) time and
 space
 */
std::pair<Tree,std::vector<int>> build_county_forest(
    const Graph &g, const arma::uvec &counties, int const num_counties
);

/*
 * Given a graph G and county assignments this creates the potentially disconnected graph
 * created when all edges across counties are removed from G. This guarantees that any
 * search started from a vertex in one county will never leave that county
 *  
 */
Graph build_restricted_county_graph(Graph const &g,  arma::uvec const &counties);

// Essentially just a useful container for map parameters 

class MapParams {
    public:
    // Constructor 
    MapParams(Rcpp::List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int const ndists, int const total_seats, std::vector<int> const &district_seat_sizes,
        double const lower, double const target, double const upper) :
        g(list_to_graph(adj_list)), counties(counties), num_counties(max(counties)),
        cg(county_graph(g, counties)),
        county_restricted_graph(num_counties > 1 ? build_restricted_county_graph(g, counties) : Graph(0)),
        // silly but call function twice so attributes can be constant 
        county_forest(build_county_forest(g, counties, num_counties).first), 
        county_forest_roots(build_county_forest(g, counties, num_counties).second),
        pop(pop),
        V(static_cast<int>(g.size())), ndists(ndists), total_seats(total_seats),
        lower(lower), target(target), upper(upper),
        smallest_district_size(*min_element(district_seat_sizes.begin(), district_seat_sizes.end())),
        largest_district_size(*max_element(district_seat_sizes.begin(), district_seat_sizes.end())),
        district_seat_sizes(district_seat_sizes),
        is_district(total_seats+1, false){
            // check the sizes are ok 
            if(ndists-1 > MAX_SUPPORTED_NUM_DISTRICTS || total_seats-1 > MAX_SUPPORTED_NUM_DISTRICTS){
                REprintf("The maximum number of districts supported right now is %u!\n",
                    MAX_SUPPORTED_NUM_DISTRICTS);
                throw Rcpp::exception("Too many districts! This number of districts isn't supported!\n");
            }
            if(num_counties-1 > MAX_SUPPORTED_NUM_COUNTIES){
                REprintf("The maximum number of counties supported right now is %u!\n",
                    MAX_SUPPORTED_NUM_COUNTIES);
                throw Rcpp::exception("Too many counties! This number of counties isn't supported!\n");
            }
            if(V-1 > MAX_SUPPORTED_NUM_VERTICES){
                REprintf("The maximum number of vertices supported right now is %u!\n",
                    MAX_SUPPORTED_NUM_VERTICES);
                throw Rcpp::exception("Too many vertices in the map! This number of vertices isn't supported!\n");
            }
            if(ndists > total_seats){
                throw Rcpp::exception("The number of distrcts must be less than or equal to the total number of seats!\n");
            }else if(ndists != total_seats){
                for (auto const &a_size: district_seat_sizes){
                    if(a_size < 0) throw Rcpp::exception("District Seat Sizes must be strictly positive!\n");
                    if(a_size >= total_seats)  throw Rcpp::exception("District Seat Sizes must be less than total seats!\n");
                    // mark this as a district size 
                    is_district[a_size] = true;
                }
                
            }else{
                is_district[1] = true;
            }
            
        };

    Graph const g; // The graph as undirected adjacency list 
    arma::uvec const counties; // county labels
    int const num_counties; // The number of distinct counties
    Multigraph const cg; // county multigraph
    Graph county_restricted_graph; // g but with all edges crossing counties removed 
    Tree const county_forest; // Spanning forest on the counties, ie each tree is a tree on a specific county
    std::vector<int> const county_forest_roots; // roots of each county tree, so [i] is root of tree on county[i+1]
    arma::uvec const pop; // population of each vertex
    int const V; // Number of vertices in the graph
    int const ndists; // The number of districts a final plan should have
    int const total_seats; // the total number of seats 
    double const lower; // lower bound on district population
    double const target; // target district population
    double const upper; // upper bound on district population
    int const smallest_district_size; // smallest district size
    int const largest_district_size; // largest district size
    std::vector<int> const &district_seat_sizes; // vector of all district seat sizes 
    std::vector<bool> is_district; // of length total_seats that says whether or not that size is a district

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
        cut_above_pop(0),
        log_prob(0) {}
    
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
          cut_above_pop(cut_above_pop),
          log_prob(0) {}
    
    // Attributes
    int tree_root; // The root of the tree
    int cut_vertex; // The vertex where we are cutting below it
    int cut_vertex_parent; // The parent of `cut_vertex` so we are cutting `(cut_vertex_parent, cut_vertex)`
    int cut_below_region_size; // The size of the region below made by cutting 
    int cut_below_pop; // The population of the region below made by cutting 
    int cut_above_region_size; // The size of the region above made by cutting 
    int cut_above_pop; // The population of the region above made by cutting 
    double log_prob; // Log Probability this edge was chosen to split in the tree 

    // Gets the information on the two regions formed from an edge cut by reference
    void get_split_regions_info(
        int &split_region1_tree_root, int &split_region1_dval, int &split_region1_pop,
        int &split_region2_tree_root, int &split_region2_dval, int &split_region2_pop
    ) const;

    // Gets the signed (not absolute value) deviation of the two regions from the targets
    // first entry is below and second is above
    std::array<double, 2> compute_signed_pop_deviances(double target) const;

    // returns absolute population deviation
    std::array<double, 2> compute_abs_pop_deviances(double target) const;

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

// enum for sampling spaces
enum class SamplingSpace : unsigned char
{
    GraphSpace, // Sampling on the space of graph partitions
    ForestSpace, // Sample on the space of spanning forests
    LinkingEdgeSpace // Sample on space of forests and linking edges
};

// loads a sampling spaces type enum from a control string
SamplingSpace get_sampling_space(std::string const &sampling_space_str);

// Get convinient string representation
std::string sampling_space_to_str(SamplingSpace sampling_space);

// enum for various methods of splitting a plan
enum class SplittingMethodType : unsigned char
{
    NaiveTopK, // picks 1 of top k edges even if invalid
    UnifValid, // picks uniform valid edge at random 
    ExpBiggerAbsDev, // propto exp(-alpha*bigger abs dev of pair)
    ExpSmallerAbsDev, // propto exp(-alpha*smaller abs dev of pair)
    Experimental // Just for testing
};

// loads a splitting type enum from a control string
SplittingMethodType get_splitting_type(std::string const &splitting_type_str);

// Get convinient string representation
std::string splitting_method_to_str(SplittingMethodType splitting_method);

enum class SplittingSizeScheduleType : unsigned char
{
    DistrictOnlySMD,
    AnyValidSizeSMD,
    DistrictOnlyMMD,
    AnyValidSizeMMD,
    OneCustomSize,
    PureMergeSplitSize,
    CustomSizes
};

// load from control spring 
SplittingSizeScheduleType get_splitting_size_regime(std::string const &splitting_size_regime_str);




#endif
