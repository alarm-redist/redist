#pragma once
#ifndef ADVANCED_REDIST_TYPES_H
#define ADVANCED_REDIST_TYPES_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "redist_constants.h"
#include "redist_types.h"

using Clock = std::chrono::steady_clock;

typedef uint16_t VertexID; // type for trees to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_VERTICES =
    static_cast<unsigned int>(std::numeric_limits<VertexID>::max());

typedef std::uint_least8_t RegionID; // type for plan vectors to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_DISTRICTS =
    static_cast<unsigned int>(std::numeric_limits<RegionID>::max());

typedef std::uint_least16_t CountyID; // type for county vectors to save space from normal int
constexpr unsigned int MAX_SUPPORTED_NUM_COUNTIES =
    static_cast<unsigned int>(std::numeric_limits<CountyID>::max());

typedef std::uint_least16_t CountyRegion; // type for region_intersect_county lookup
constexpr unsigned int MAX_SUPPORTED_COUNTYREGION_VALUE =
    static_cast<unsigned int>(std::numeric_limits<CountyRegion>::max());

//
typedef std::vector<RegionID> AllPlansVector;
typedef std::vector<RegionID> AllRegionSizesVector;
typedef std::vector<RegionID> RegionSizeVector;
typedef std::tuple<CountyRegion, RegionID, CountyID> CountyComponentVertex;
typedef std::vector<std::vector<CountyComponentVertex>> CountyComponentGraph;

// -----------------------------------------------------------------------------
// Packed forest edge storage
// The idea is we have a list of all the edges in the graph and represent
// a forest as a boolean vector where true means the edge is in the forest and 
// false means not 
// -----------------------------------------------------------------------------
typedef std::uint32_t EdgeID;       // identifies edge number e
typedef std::uint64_t EdgeBitWord;  // stores 64 edge booleans

constexpr int EDGE_BITS_PER_WORD = 64;

inline int compute_num_edge_bit_words(int num_graph_edges) {
    return (num_graph_edges + EDGE_BITS_PER_WORD - 1) / EDGE_BITS_PER_WORD;
}



// This is a template class for storing plan attributes
// in a way to minimize memory
template <typename T> class PlanAttribute {
  private:
    int const offset_start;   // starting index for this subset of vector
    int const offset_end;     // end index for this subset of vector
    std::vector<T> &long_vec; // ref to underlying long vector

  public:
    PlanAttribute(std::vector<T> &underyling_long_vec, int offset_start, int offset_end)
        : offset_start(offset_start), offset_end(offset_end), long_vec(underyling_long_vec) {};

    // Checks if empty 
    bool empty() const { return offset_start == offset_end; }
    // methods for accessing
    // Const version for read-only access
    const T operator[](int index) const { 
        if constexpr(perf_config::bounds_checking){
            if (index < 0 || index >= static_cast<int>(size())) {
                std::ostringstream oss;
                oss << "PlanAttribute::operator[] const out of range.\n";
                oss << "index=" << index << "\n";
                oss << "size=" << size() << "\n";
                oss << "offset_start=" << offset_start << "\n";
                oss << "offset_end=" << offset_end << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        return long_vec[offset_start + index]; 
    };
    // non constant for modification
    T &operator[](int index) { 
        if constexpr(perf_config::bounds_checking){
            if (index < 0 || index >= static_cast<int>(size())) {
                std::ostringstream oss;
                oss << "PlanAttribute::operator[] out of range.\n";
                oss << "index=" << index << "\n";
                oss << "size=" << size() << "\n";
                oss << "offset_start=" << offset_start << "\n";
                oss << "offset_end=" << offset_end << "\n";
                throw std::runtime_error(oss.str());
            }
        }
        return long_vec[offset_start + index]; 
    };

    // Non-const iterator accessors
    auto begin() { return long_vec.begin() + offset_start; }
    auto end() { return long_vec.begin() + offset_end; }

    // Const iterator accessors
    auto begin() const { return long_vec.begin() + offset_start; }
    auto end() const { return long_vec.begin() + offset_end; }

    // This copies one plans data from another
    void copy(PlanAttribute const &other_attr) {
        if constexpr (perf_config::bounds_checking){
            if (size() != other_attr.size()) {
                std::ostringstream oss;
                oss << "PlanAttribute::copy size mismatch. "
                    << "dest.size()=" << size()
                    << ", source.size()=" << other_attr.size();
                throw std::runtime_error(oss.str());
            }
        }
        std::copy(other_attr.begin(), other_attr.end(), begin());
    }

    std::size_t size() const noexcept { return offset_end - offset_start; };

    // debug functions 
    int debug_offset_start() const {
        return offset_start;
    }

    int debug_offset_end() const {
        return offset_end;
    }

    std::size_t debug_size() const {
        return static_cast<std::size_t>(offset_end - offset_start);
    }

    std::vector<T> const *debug_backing_vector_address() const {
        return &long_vec;
    }

    T const *debug_data_begin() const {
        return long_vec.data() + offset_start;
    }

    T const *debug_data_end() const {
        return long_vec.data() + offset_end;
    }

    bool debug_overlaps(PlanAttribute<T> const &other) const {
        if (&long_vec != &other.long_vec) {
            return false;
        }

        return offset_start < other.offset_end &&
            other.offset_start < offset_end;
    }

};

typedef PlanAttribute<RegionID> PlanVector;
typedef PlanAttribute<RegionID> RegionSizes;
typedef PlanAttribute<int> IntPlanAttribute;
typedef PlanAttribute<double> DoublePlanAttribute;
typedef PlanAttribute<EdgeBitWord> PlanEdgeBits;

// uniquely maps pairs (x,y) of the form
// 0 <= x < y < container_size
// to a value in (0, container_size choose 2)
// Currently used in the following scenarios
//  - Pairs of regions (region1_id, region2_id) with 0 <= region1_id < region2_id < num_regions
inline int index_from_ordered_pair(int x, int y, int container_size) {
    return ((x * (2 * container_size - x - 1)) / 2) + (y - x - 1);
}

// indexing (i,j) in a matrix
// assuming the length of a row is row_length
// i is the row index
// j is the column index so j < row_length
inline int mat_index_from_pair(int i, int j, int row_length) { return i * row_length + j; }

// Class for storing the graph as a long vector of edges 
// Allows you to take a vertex and get neighbors 
class GraphEdgeIndex {
  public:
  // incident edge stores both the vertex adjacent to v and the edge_id associated with this edge
    struct IncidentEdge {
        VertexID neighbor;
        EdgeID edge_id;
    };

    int const num_edges;
    int const V;
    int max_num_edges; // TODO make constant and initialized lated 

    GraphEdgeIndex() = delete;

    explicit GraphEdgeIndex(Graph const &g, int const num_edges);

    // Takes two vertices and returns their edge id
    EdgeID get_edge_id(int v, int u) const {
        if constexpr (perf_config::bounds_checking){
            check_vertex(v, "GraphEdgeIndex::get_edge_id received invalid v!");
            check_vertex(u, "GraphEdgeIndex::get_edge_id received invalid u!");
        }
        // 
        VertexID const target = static_cast<VertexID>(u);

        for (auto const &incident_edge : incident_edges[v]) {
            if (incident_edge.neighbor == target) {
                return incident_edge.edge_id;
            }
        }

        std::ostringstream oss;
        oss << "GraphEdgeIndex::get_edge_id called on non-edge. "
            << "v=" << v
            << ", u=" << u
            << ", V=" << incident_edges.size()
            << ". Neighbors of v: ";

        if (v >= 0 && v < static_cast<int>(incident_edges.size())) {
            for (auto const &incident_edge : incident_edges[v]) {
                oss << static_cast<int>(incident_edge.neighbor) << " ";
            }
        }

        throw std::runtime_error(oss.str());
    }

    bool has_edge(int v, int u) const {
        if constexpr (perf_config::bounds_checking){
            if (!is_valid_vertex(v) || !is_valid_vertex(u)) {
                return false;
            }
        }

        VertexID const target = static_cast<VertexID>(u);

        for (auto const &incident_edge : incident_edges[v]) {
            if (incident_edge.neighbor == target) {
                return true;
            }
        }

        return false;
    }

    // takes an edge id and return the pair associated with it
    std::pair<VertexID, VertexID> get_edge_endpoints(EdgeID edge_id) const {
        if constexpr (perf_config::bounds_checking){
            if (static_cast<std::size_t>(edge_id) >= edges.size()) {
                throw std::runtime_error("GraphEdgeIndex::get_edge_endpoints received invalid edge_id!");
            }
        }

        return edges[edge_id];
    }

    std::vector<std::pair<VertexID, VertexID>> edges;
    std::vector<std::vector<IncidentEdge>> incident_edges;

  private:
    bool is_valid_vertex(int v) const {
        return v >= 0 && v < static_cast<int>(incident_edges.size());
    }

    void check_vertex(int v, char const *message) const {
        if (!is_valid_vertex(v)) {
            throw std::runtime_error(message);
        }
    }

};


// Designed to allow for different tree splitting methods
// This allows us to seperate cutting the tree from finding the edge to cut
class EdgeCut {

  public:
    // Default Constructor
    EdgeCut()
        : tree_root(0), cut_vertex(0), cut_vertex_parent(0), cut_below_region_size(0),
          cut_below_pop(0), cut_above_region_size(0), cut_above_pop(0), log_prob(0) {}

    // Constructor
    EdgeCut(const int tree_root, const int cut_vertex, const int cut_vertex_parent,
            const int cut_below_region_size, const int cut_below_pop,
            const int cut_above_region_size, const int cut_above_pop)
        : tree_root(tree_root), cut_vertex(cut_vertex), cut_vertex_parent(cut_vertex_parent),
          cut_below_region_size(cut_below_region_size), cut_below_pop(cut_below_pop),
          cut_above_region_size(cut_above_region_size), cut_above_pop(cut_above_pop),
          log_prob(0) {}

    // Attributes
    int tree_root;         // The root of the tree
    int cut_vertex;        // The vertex where we are cutting below it
    int cut_vertex_parent; // The parent of `cut_vertex` so we are cutting `(cut_vertex_parent,
                           // cut_vertex)`
    int cut_below_region_size; // The size of the region below made by cutting
    int cut_below_pop;         // The population of the region below made by cutting
    int cut_above_region_size; // The size of the region above made by cutting
    int cut_above_pop;         // The population of the region above made by cutting
    double log_prob;           // Log Probability this edge was chosen to split in the tree

    // Gets the information on the two regions formed from an edge cut by reference
    void get_split_regions_info(int &split_region1_tree_root, int &split_region1_dval,
                                int &split_region1_pop, int &split_region2_tree_root,
                                int &split_region2_dval, int &split_region2_pop) const;

    // Gets the signed (not absolute value) deviation of the two regions from the targets
    // first entry is below and second is above
    std::array<double, 2> compute_signed_pop_deviances(double target) const;

    // returns absolute population deviation
    std::array<double, 2> compute_abs_pop_deviances(double target) const;

    // Equality operator
    bool operator==(const EdgeCut &other) const {
        return tree_root == other.tree_root && cut_vertex == other.cut_vertex &&
               cut_vertex_parent == other.cut_vertex_parent &&
               cut_below_region_size == other.cut_below_region_size &&
               cut_below_pop == other.cut_below_pop &&
               cut_above_region_size == other.cut_above_region_size &&
               cut_above_pop == other.cut_above_pop;
    }

    // Not-equal operator
    bool operator!=(const EdgeCut &other) const { return !(*this == other); }

    // Less-than operator
    bool operator<(const EdgeCut &other) const { return cut_vertex < other.cut_vertex; }
};

// Flat representation of graphs where its 
// fixed at runtime 
class FlatGraph {
  private:
    std::vector<int> offsets;
    std::vector<int> sizes;
    std::vector<int> data;

  public:
    struct NeighborRange {
            int const *ptr;
            int len;

            int const *begin() const { return ptr; }
            int const *end() const { return ptr + len; }

            int size() const { return len; }
            bool empty() const { return len == 0; }

            int operator[](int i) const { 
                if constexpr (perf_config::bounds_checking) {
                    if (i >= len || i < 0) {
                        throw std::out_of_range(
                            "FlatGraph::NeighborRange index out of range."
                        );
                    }
                }
                return ptr[i]; 
            }
    };

    FlatGraph() = default;

    explicit FlatGraph(Graph const &g) {
        std::size_t const V = g.size();

        offsets.resize(V + 1);
        sizes.assign(V, 0);
        
        std::size_t total_cap = 0;
        for (std::size_t v = 0; v < V; ++v) {
            // We set the start of v's possible neighbors at the current
            // total cap
            offsets[v] = total_cap;
            // see the degree of this vertex
            auto const v_capacity = g[v].size();
            // add this to the total capacity
            total_cap += v_capacity;
            // also store the size in g
            sizes[v] = v_capacity;
        }
        offsets[V] = total_cap;

        data.resize(total_cap);
        // Fill flat storage with the original graph neighbors.
        for (std::size_t v = 0; v < V; v++)
        {
            auto const start = offsets[v];

            for (int j = 0; j < sizes[v]; ++j) {
                data[start + j] = g[v][j];
            }
        }
    }

    // Returns a flat empty tree on the current flat 
    FlatGraph get_flat_empty_tree() const {
        FlatGraph out;

        out.offsets = offsets;
        out.sizes.assign(sizes.size(), 0);
        out.data.resize(data.size());

        return out;
    }

    int size() const {
        return static_cast<int>(sizes.size());
    }

    void clear_vertex(int v) {
        sizes[v] = 0;
    }

    void clear() {
        std::fill(sizes.begin(), sizes.end(), 0);
    }

    
    void add_directed_edge(int const parent, int const child) {
        if constexpr (perf_config::bounds_checking) {
            if (parent < 0 || parent >= size()) {
                std::ostringstream oss;
                oss << "FlatGraph::add_directed_edge invalid parent.\n";
                oss << "parent=" << parent << "\n";
                oss << "child=" << child << "\n";
                oss << "size=" << size() << "\n";
                throw std::runtime_error(oss.str());
            }

            if (child < 0 || child >= size()) {
                std::ostringstream oss;
                oss << "FlatGraph::add_directed_edge invalid child.\n";
                oss << "parent=" << parent << "\n";
                oss << "child=" << child << "\n";
                oss << "size=" << size() << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        
        auto const pos = sizes[parent];

        if constexpr (perf_config::bounds_checking) {
            auto const v = parent;
            std::size_t const cap = offsets[v + 1] - offsets[v];

            if (pos >= cap) {
                std::ostringstream oss;
                oss << "FlatGraph::add_directed_edge exceeded capacity.\n";
                oss << "parent=" << parent << "\n";
                oss << "child=" << child << "\n";
                oss << "size=" << pos << "\n";
                oss << "cap=" << cap << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        data[offsets[parent] + pos] = child;
        sizes[parent] = pos + 1;
    }

    // adds an undirected edge to the graph 
    void add_undirected_edge(int const v, int const u) {
        add_directed_edge(v, u);
        add_directed_edge(u, v);
    }

    int degree(int v) const {
        return sizes[v];
    }

    int neighbor(int const v, int const j) const {
        if constexpr (perf_config::bounds_checking) {
            if (v < 0 || v >= size()) {
                throw std::runtime_error("FlatGraph::neighbor invalid v.");
            }

            if (j < 0 || j >= sizes[v]) {
                throw std::runtime_error("FlatGraph::neighbor invalid j.");
            }
        }

        return data[offsets[v] + j];
    }

    NeighborRange neighbors(int const v) const {
        auto const start = offsets[v];
        return NeighborRange{data.data() + start, sizes[v]};
    }

    // converts to a vertex graph. Useful for exporting to R
    Graph to_vertex_graph() const {
        Graph out(size());

        for (std::size_t v = 0; v < size(); ++v) {
            out[v].reserve(sizes[v]);

            auto const start = offsets[v];
            auto const stop = start + sizes[v];

            for (std::size_t idx = start; idx < stop; ++idx) {
                out[v].push_back(data[idx]);
            }
        }

        return out;
    }

    // erases the DIRECTED edge, meaning it only 
    // erases the directed edge (parent, child)
    void erase_directed_edge(EdgeCut const &cut_edge);
};

// enum for sampling spaces
enum class SamplingSpace : unsigned char {
    GraphSpace,       // Sampling on the space of graph partitions
    ForestSpace,      // Sample on the space of spanning forests
    LinkingEdgeSpace, // Sample on space of forests and linking edges
    LCTGraphSpace     // Graph partition with a link-cut-tree spanning forest (cyclewalk)
};

// loads a sampling spaces type enum from a control string
SamplingSpace get_sampling_space(std::string const &sampling_space_str);

// Get convinient string representation
std::string sampling_space_to_str(SamplingSpace sampling_space);

// enum for various methods of splitting a plan
enum class SplittingMethodType : unsigned char {
    NaiveTopK,        // picks 1 of top k edges even if invalid
    UnifValid,        // picks uniform valid edge at random
    ExpBiggerAbsDev,  // propto exp(-alpha*bigger abs dev of pair)
    ExpSmallerAbsDev, // propto exp(-alpha*smaller abs dev of pair)
    Constraint,       // propto constraint score (unif if no constraints)
    Experimental      // Just for testing
};

// loads a splitting type enum from a control string
SplittingMethodType get_splitting_type(std::string const &splitting_type_str);

// Get convinient string representation
std::string splitting_method_to_str(SplittingMethodType splitting_method);

enum class SplittingSizeScheduleType : unsigned char {
    DistrictOnlySMD,
    AnyValidSizeSMD,
    DistrictOnlyMMD,
    AnyValidSizeMMD,
    OneCustomSize,
    PureMergeSplitSize,
    CustomSizes
};

// load from control spring
SplittingSizeScheduleType
get_splitting_size_regime(std::string const &splitting_size_regime_str);

// Essentially just a useful container for map and some algorithm parameters
class MapParams {
  public:
    // Constructor
    MapParams(Graph const &g, 
              const std::vector<unsigned int> &counties, 
              const std::vector<unsigned int> &pop,
              int const ndists, int const total_seats,
              std::vector<int> const &district_seat_sizes, 
              double const lower,
              double const target, double const upper,
              SamplingSpace const sampling_space);

    // Constructor for when we only need map information
    // so fake district info is ued 
    MapParams(Graph const &g, 
              const std::vector<unsigned int> &counties, 
              const std::vector<unsigned int> &pop,
              double const lower, double const upper
             );



    Graph const g;                       // The graph as undirected adjacency list
    FlatGraph map_graph;                 // The graph stored as a FlatGraph type
    int const num_edges;                 // number of undirected edges in g
    GraphEdgeIndex const graph_edge_index; // bitpacked boolean storage of graph
    int const num_edge_bit_words; // used for bitpacked stuff
    std::vector<unsigned int> const counties;           // county labels
    int const num_counties;              // The number of distinct counties
    Multigraph const cg;                 // county multigraph
    Graph const county_restricted_graph; // g but with all edges crossing counties removed
    FlatGraph county_restricted_flat_graph; // flat version
    std::vector<unsigned int> const pop; // population of each vertex
    int const V;                         // Number of vertices in the graph
    int const ndists;                    // The number of districts a final plan should have
    int const total_seats;               // the total number of seats
    double const lower;                  // lower bound on district population
    double const target;                 // target district population
    double const upper;                  // upper bound on district population
    int const smallest_district_size;    // smallest district size
    int const largest_district_size;     // largest district size
    std::vector<int> const district_seat_sizes; // vector of all district seat sizes
    std::vector<bool> const
        is_district;   // of length total_seats that says whether or not that size is a district
    bool const is_mmd; // Whether or not multimember districting
    SamplingSpace const sampling_space;
};



// Stack where maximum size is known at runtime

template <typename T> class FixedStack {
  private:
    std::vector<T> buffer;
    size_t capacity;
    size_t count = 0;

  public:
    explicit FixedStack(size_t max_size) : buffer(max_size), capacity(max_size) {}

    bool empty() const { return count == 0; }
    bool full() const { return count == capacity; }
    size_t size() const { return count; }
    size_t max_size() const { return capacity; }

    void push(const T &value) {
        if constexpr (perf_config::bounds_checking){
            if (count >= capacity) {
                throw std::runtime_error("FixedStack overflow.");
            }
        }
        buffer[count++] = value; // copy
    }

    void push(T &&value) {
        if constexpr (perf_config::bounds_checking){
            if (count >= capacity) {
                throw std::runtime_error("FixedStack overflow.");
            }
        }
        buffer[count++] = std::move(value); // move
    }

    T &top() { 
        if constexpr (perf_config::bounds_checking){
            if (count == 0) {
                throw std::runtime_error("FixedStack underflow.");
            }
        }
        return buffer[count - 1]; 
    }

    const T &top() const { 
        if constexpr (perf_config::bounds_checking){
            if (count == 0) {
                throw std::runtime_error("FixedStack underflow.");
            }
        }
        return buffer[count - 1]; }

    T pop() { 
        if constexpr (perf_config::bounds_checking){
            if (count == 0) {
                throw std::runtime_error("FixedStack underflow.");
            }
        }
        return std::move(buffer[--count]); 
    }

    void clear() { count = 0; }
};

// type for main splitting code
typedef FixedStack<std::tuple<int, int, bool>> TreePopStack;

template <typename T> class CircularQueue {
  private:
    std::vector<T> buffer;
    size_t head = 0; // Points to the front item
    size_t tail = 0; // Points to the next write position
    size_t const capacity;
    size_t size = 0;

  public:
    explicit CircularQueue(size_t max_size) : buffer(max_size), capacity(max_size) {}

    bool empty() const { return size == 0; }
    bool full() const { return size == capacity; }
    size_t get_size() const { return size; }
    size_t max_size() const { return capacity; }

    void push(const T &value) {
        if constexpr (perf_config::bounds_checking){
            if (size >= capacity) {
                throw std::runtime_error("CircularQueue overflow.");
            }
        }

        buffer[tail] = value;
        tail = (tail + 1) % capacity;
        ++size;
    }

    void push(T &&value) {
        if constexpr (perf_config::bounds_checking){
            if (size >= capacity) {
                throw std::runtime_error("CircularQueue overflow.");
            }
        }
        buffer[tail] = std::move(value);
        tail = (tail + 1) % capacity;
        ++size;
    }

    T pop() {
        if constexpr (perf_config::bounds_checking){
            if (size == 0) {
                throw std::runtime_error("CircularQueue underlow.");
            }
        }

        T value = std::move(buffer[head]);
        head = (head + 1) % capacity;
        --size;
        return value;
    }

    T &front() { return buffer[head]; }

    const T &front() const { return buffer[head]; }

    T &back() { return buffer[(tail + capacity - 1) % capacity]; }

    const T &back() const { return buffer[(tail + capacity - 1) % capacity]; }

    void clear() { head = tail = size = 0; }
};

typedef CircularQueue<int> DummyTreeQueue;





#endif
