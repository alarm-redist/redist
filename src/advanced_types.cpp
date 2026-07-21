/********************************************************
 * Author: Philip O'Sullivan
 * Institution: Harvard University
 * Date Created: 2026/07
 * Purpose: Contains classes and types used in the new (Redist 5.0)
 * code. Seperated out from RedistTypes to reduce compile time since a 
 * lot of legacy code must link to redist_types
 ********************************************************/

#include "advanced_types.h"

#include <cmath>
#include <unordered_set>

namespace{



/*
 * Initialize empty multigraph structure on graph with `V` vertices
 */
// TESTED
Multigraph init_multigraph(int V) {
    Multigraph g;
    for (int i = 0; i < V; i++) {
        std::vector<std::array<int, 3>> el;
        g.push_back(el);
    }
    return g;
}

/*
 * Make a county graph from a precinct graph and list of counties
 * County graph is list of list of 3: <cty of nbor, index of vtx, index of nbor>
 */
// TESTED
Multigraph county_graph(const Graph &g, const std::vector<unsigned int> &counties) {
    int n_county = *std::max_element(counties.begin(), counties.end());
    Multigraph cg = init_multigraph(n_county);

    int V = g.size();
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        int county = counties.at(i) - 1;
        for (int j = 0; j < length; j++) {
            int nbor_cty = counties.at(nbors[j]) - 1;
            if (county == nbor_cty)
                continue;
            std::array<int, 3> el = {nbor_cty, i, nbors[j]};
            cg.at(county).push_back(el);
        }
    }

    return cg;
}

/*
 * Given a graph G and county assignments this creates the potentially disconnected graph
 * created when all edges across counties are removed from G. This guarantees that any
 * search started from a vertex in one county will never leave that county
 *
 */
Graph build_restricted_county_graph(Graph const &g, std::vector<unsigned int> const &counties) {
    Graph county_graph(g.size());
    // iterate through g and only add edges in the same county
    for (int v = 0; v < g.size(); v++) {
        // iterate over v's neighbors
        for (const auto u : g[v]) {
            // if same county add the edge
            if (counties[v] == counties[u]) {
                county_graph[v].push_back(u);
            }
        }
    }
    return county_graph;
}


// Counts the number of undirected edges in a graph. 
// It also checks the graph is actually symmetric 
int count_undirected_edges(Graph const &g) {
    int edge_count = 0;
    int total_edge_count = 0;

    for (int v = 0; v < static_cast<int>(g.size()); ++v) {
        for (int u : g[v]) {
            if (u < 0 || u >= static_cast<int>(g.size())) {
                throw std::runtime_error("GraphEdgeIndex found invalid neighbor index!");
            }
            total_edge_count++; // this counts total edges

            if (v < u) {
                ++edge_count;
            }
        }
    }

    if (edge_count > static_cast<std::size_t>(std::numeric_limits<EdgeID>::max())) {
        throw std::runtime_error("Too many graph edges for EdgeID!\n");
    }

    if (edge_count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Too many graph edges for int!\n");
    }

    if (total_edge_count != edge_count * 2){
        throw std::runtime_error("The inputted adjacency list is not symmetric!\n");
    }

    return edge_count;
}

}

void EdgeCut::get_split_regions_info(int &split_region1_tree_root, int &split_region1_dval,
                                     int &split_region1_pop, int &split_region2_tree_root,
                                     int &split_region2_dval, int &split_region2_pop) const {
    // Always make region 1 the smaller one by size (allowing for ties)

    if (cut_below_region_size <= cut_above_region_size) {
        // if true then cut below is smalelr so make that region 1
        split_region1_tree_root = cut_vertex;
        split_region1_dval = cut_below_region_size;
        split_region1_pop = cut_below_pop;

        split_region2_tree_root = tree_root;
        split_region2_dval = cut_above_region_size;
        split_region2_pop = cut_above_pop;
    } else {
        // if false then cut above is smalelr so make that region 1
        split_region2_tree_root = cut_vertex;
        split_region2_dval = cut_below_region_size;
        split_region2_pop = cut_below_pop;

        split_region1_tree_root = tree_root;
        split_region1_dval = cut_above_region_size;
        split_region1_pop = cut_above_pop;
    }
};

std::array<double, 2> EdgeCut::compute_signed_pop_deviances(double target) const {
    // get the target populations for the regions
    double cut_below_target = target * cut_below_region_size;
    double cut_above_target = target * cut_above_region_size;
    // get the deviation
    double below_dev =
        (static_cast<double>(cut_below_pop) - cut_below_target) / cut_below_target;
    double above_dev =
        (static_cast<double>(cut_above_pop) - cut_above_target) / cut_above_target;

    return std::array<double, 2>{below_dev, above_dev};
}

std::array<double, 2> EdgeCut::compute_abs_pop_deviances(double target) const {
    // get the raw unsigned deviations
    std::array<double, 2> unsigned_devs = compute_signed_pop_deviances(target);
    // take the absolute value
    std::array<double, 2> signed_devs = {std::fabs(unsigned_devs.at(0)),
                                         std::fabs(unsigned_devs.at(1))};

    return signed_devs;
}

// loads a sampling spaces type enum from a control string
SamplingSpace get_sampling_space(std::string const &sampling_space_str) {
    // find the type or throw an error
    if (sampling_space_str == "graph_plan") {
        return SamplingSpace::GraphSpace;
    } else if (sampling_space_str == "spanning_forest") {
        return SamplingSpace::ForestSpace;
    } else if (sampling_space_str == "linking_edge") {
        return SamplingSpace::LinkingEdgeSpace;
    } else if (sampling_space_str == "lct_graph") {
        return SamplingSpace::LCTGraphSpace;
    } else {
        throw std::invalid_argument(
            "Invalid sampling space '" + sampling_space_str +
            "'. Expected one of: graph_plan, spanning_forest, "
            "linking_edge, or lct_graph."
        );
    }
}

// Get convinient string representation
std::string sampling_space_to_str(SamplingSpace sampling_space) {
    if (sampling_space == SamplingSpace::GraphSpace) {
        return "Graph";
    } else if (sampling_space == SamplingSpace::ForestSpace) {
        return "Forest";
    } else if (sampling_space == SamplingSpace::LinkingEdgeSpace) {
        return "Linking Edge";
    } else if (sampling_space == SamplingSpace::LCTGraphSpace) {
        return "LCT Graph";
    } else {
        throw std::logic_error(
            "sampling_space_to_str received an unknown SamplingSpace value: " +
            std::to_string(static_cast<unsigned int>(sampling_space))
        );
    }
}

SplittingMethodType get_splitting_type(std::string const &splitting_type_str) {
    // find the type or throw an error
    if (splitting_type_str == "top_k") {
        return SplittingMethodType::NaiveTopK;
    } else if (splitting_type_str == "unif_valid") {
        return SplittingMethodType::UnifValid;
    } else if (splitting_type_str == "exp_abs_dev") {
        return SplittingMethodType::ExpBiggerAbsDev;
    } else if (splitting_type_str == "expo_smaller_abs_dev") {
        return SplittingMethodType::ExpSmallerAbsDev;
    } else if (splitting_type_str == "experimental") {
        return SplittingMethodType::Experimental;
    } else if (splitting_type_str == "constraint") {
        return SplittingMethodType::Constraint;
    } else {
        throw std::invalid_argument("Invalid splitting_type_str");
    }
}

std::string splitting_method_to_str(SplittingMethodType splitting_method) {
    if (splitting_method == SplittingMethodType::NaiveTopK) {
        return "Naive Top K Splitter";
    } else if (splitting_method == SplittingMethodType::UnifValid) {
        return "Uniform Valid Edge Splitter";
    } else if (splitting_method == SplittingMethodType::ExpBiggerAbsDev) {
        return "Exponentially Weighted Absolute Bigger Deviance Splitter";
    } else if (splitting_method == SplittingMethodType::ExpSmallerAbsDev) {
        return "Exponentially Weighted Absolute Smaller Deviance Splitter";
    } else if (splitting_method == SplittingMethodType::Experimental) {
        return "Experimental Splitter";
    } else if (splitting_method == SplittingMethodType::Constraint) {
        return "Constraint Splitter";
    } else {
        throw std::logic_error(
            "splitting_method_to_str received an unknown SplittingMethodType value: " +
            std::to_string(static_cast<unsigned int>(splitting_method))
        );
    }
}

SplittingSizeScheduleType
get_splitting_size_regime(std::string const &splitting_size_regime_str) {
    // find the type or throw an error
    if (splitting_size_regime_str == "split_district_only") {
        return SplittingSizeScheduleType::DistrictOnlySMD;
    } else if (splitting_size_regime_str == "any_valid_sizes") {
        return SplittingSizeScheduleType::AnyValidSizeSMD;
    } else if (splitting_size_regime_str == "split_district_only_mmd") {
        return SplittingSizeScheduleType::DistrictOnlyMMD;
    } else if (splitting_size_regime_str == "any_valid_sizes_mmd") {
        return SplittingSizeScheduleType::AnyValidSizeMMD;
    } else if (splitting_size_regime_str == "one_custom_size") {
        return SplittingSizeScheduleType::OneCustomSize;
    } else if (splitting_size_regime_str == "pure_ms_size") {
        return SplittingSizeScheduleType::PureMergeSplitSize;
    } else if (splitting_size_regime_str == "custom") {
        return SplittingSizeScheduleType::CustomSizes;
    } else {
        throw std::invalid_argument(
            "Invalid splitting size regime: '" +
            splitting_size_regime_str + "'"
        );
    }
};


GraphEdgeIndex::GraphEdgeIndex(Graph const &g, int const num_edges)
    :  num_edges(num_edges), V(g.size()), incident_edges(g.size()) {

    if (g.size() > MAX_SUPPORTED_NUM_VERTICES) {
        throw std::invalid_argument("Too many vertices for VertexID in GraphEdgeIndex!");
    }

    std::unordered_set<std::uint64_t> seen_edges;
    seen_edges.reserve(static_cast<std::size_t>(num_edges) * 2);

    int const V = static_cast<int>(g.size());
    max_num_edges = 1;

    for (int v = 0; v < V; ++v) {
        auto num_v_edges = g[v].size();
        incident_edges[v].reserve(num_v_edges);

        // update the max number of edges if this is larger
        if (num_v_edges > max_num_edges){
            max_num_edges = num_v_edges;
        }
        

        std::unordered_set<int> seen_neighbors_for_v;
        seen_neighbors_for_v.reserve(g[v].size());

        for (int u : g[v]) {
            if (u < 0 || u >= V) {
                throw std::runtime_error("GraphEdgeIndex found invalid neighbor index!");
            }

            if (!seen_neighbors_for_v.insert(u).second) {
                std::ostringstream oss;
                oss << "Duplicate neighbor in graph adjacency list: "
                    << "vertex " << v << " has neighbor " << u << " more than once.";
                
                throw std::runtime_error(oss.str());
            }

            if (v < u) {
                std::uint64_t const key =
                    (static_cast<std::uint64_t>(v) << 32) |
                    static_cast<std::uint32_t>(u);

                if (!seen_edges.insert(key).second) {
                    std::ostringstream oss;
                    oss << "Duplicate undirected edge in GraphEdgeIndex: ("
                        << v << ", " << u << ")";
                    
                    throw std::runtime_error(oss.str());
                }

                if (edges.size() > std::numeric_limits<EdgeID>::max()) {
                    throw std::runtime_error("Too many graph edges for EdgeID!");
                }

                EdgeID const eid = static_cast<EdgeID>(edges.size());

                edges.push_back({
                    static_cast<VertexID>(v),
                    static_cast<VertexID>(u)
                });

                incident_edges[v].push_back({
                    static_cast<VertexID>(u),
                    eid
                });

                incident_edges[u].push_back({
                    static_cast<VertexID>(v),
                    eid
                });
            }
        }
    }

    if (static_cast<int>(edges.size()) != num_edges) {
        std::ostringstream oss;
        oss << "GraphEdgeIndex edge count mismatch. "
            << "expected " << num_edges
            << ", built " << edges.size();
        
        throw std::runtime_error(oss.str());
    }
}

void FlatGraph::push_back(int v, int v_nbor_to_add) {
        if constexpr (perf_config::bounds_checking) {
            if (v < 0 || v >= size()) {
                throw std::runtime_error("FlatGraph::push_back invalid v");
            }
        }

        auto const pos = sizes[v];

        if constexpr (perf_config::bounds_checking) {
            VertexID const cap = offsets[v + 1] - offsets[v];

            if (pos >= cap) {
                std::ostringstream oss;
                oss << "FlatGraph::push_back exceeded capacity. "
                    << "v=" << v
                    << ", size=" << pos
                    << ", cap=" << cap
                    << ", u=" << v_nbor_to_add;
                throw std::runtime_error(oss.str());
            }
        }

        data[offsets[v] + pos] = static_cast<VertexID>(v_nbor_to_add);
        sizes[v] = pos + 1;
}

// adds a directed edge 
void FlatGraph::add_directed_edge(int const parent, int const child){
    push_back(parent, child);
}

// adds an undirected edge to the graph 
void FlatGraph::add_undirected_edge(int const v, int const u) {
    push_back(v, u);
    push_back(u, v);
}

void FlatGraph::erase_directed_edge(EdgeCut const &cut_edge) {
    int const parent = cut_edge.cut_vertex_parent;
    int const child = cut_edge.cut_vertex;

    if constexpr (perf_config::supposedly_safe_input_checks) {
        int const V = size();

        if (parent < 0 || parent >= V) {
            std::ostringstream oss;
            oss << "FlatGraph::erase_tree_edge invalid cut_vertex_parent.\n";
            oss << "cut_vertex_parent=" << parent << "\n";
            oss << "cut_vertex=" << child << "\n";
            oss << "tree_root=" << cut_edge.tree_root << "\n";
            oss << "FlatGraph size=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (child < 0 || child >= V) {
            std::ostringstream oss;
            oss << "FlatGraph::erase_tree_edge invalid cut_vertex.\n";
            oss << "cut_vertex_parent=" << parent << "\n";
            oss << "cut_vertex=" << child << "\n";
            oss << "tree_root=" << cut_edge.tree_root << "\n";
            oss << "FlatGraph size=" << V << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    VertexID const start = offsets[parent];
    VertexID const stop = start + sizes[parent];

    VertexID pos = stop;

    for (VertexID idx = start; idx < stop; ++idx) {
        if (data[idx] == static_cast<VertexID>(child)) {
            pos = idx;
            break;
        }
    }

    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (pos == stop) {
            std::ostringstream oss;
            oss << "FlatGraph::erase_tree_edge could not find cut edge.\n";
            oss << "cut_vertex_parent=" << parent << "\n";
            oss << "cut_vertex=" << child << "\n";
            oss << "tree_root=" << cut_edge.tree_root << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    // Shift later children left by one, preserving order.
    for (VertexID idx = pos + 1; idx < stop; ++idx) {
        data[idx - 1] = data[idx];
    }

    --sizes[parent];
}

// Constructor
MapParams::MapParams(Graph const &g, 
    const std::vector<unsigned int> &counties, 
    const std::vector<unsigned int> &pop,
            int const ndists, int const total_seats,
            std::vector<int> const &district_seat_sizes, 
            double const lower,
            double const target, double const upper,
            SamplingSpace const sampling_space)
    : g(g), map_graph(g),
    num_edges(count_undirected_edges(g)),
    graph_edge_index(g, num_edges), num_edge_bit_words(compute_num_edge_bit_words(num_edges)),
        counties(counties.begin(), counties.end()), 
        num_counties(*std::max_element(counties.begin(), counties.end())), 
        cg(county_graph(g, counties)),
        county_restricted_graph(num_counties > 1 ? build_restricted_county_graph(g, counties)
                                                : Graph(0)),
        pop(pop), V(static_cast<int>(g.size())), ndists(ndists), total_seats(total_seats),
        lower(lower), target(target), upper(upper),
        smallest_district_size(
            *std::min_element(district_seat_sizes.begin(), district_seat_sizes.end())),
        largest_district_size(
            *std::max_element(district_seat_sizes.begin(), district_seat_sizes.end())),
        district_seat_sizes(district_seat_sizes),
        is_district([ndists, total_seats, &district_seat_sizes]() {
            if (district_seat_sizes.size() == 0) {
                throw std::runtime_error("District Seat Sizes vector must be non-empty!\n");
            }
            // vector where index i is true iff i seats is a district
            std::vector<bool> is_district_vec(total_seats + 1, false);
            if (ndists > total_seats) {
                throw std::runtime_error("The number of distrcts must be less than or equal to "
                                    "the total number of seats!\n");
            } else if (ndists != total_seats) {
                for (auto const &a_size : district_seat_sizes) {
                    if (a_size < 0)
                        throw std::runtime_error(
                            "District Seat Sizes must be strictly positive!\n");
                    if (a_size >= total_seats)
                        throw std::runtime_error(
                            "District Seat Sizes must be less than total seats!\n");
                    // mark this as a district size
                    is_district_vec[a_size] = true;
                }

            } else {
                is_district_vec[1] = true;
            }
            return is_district_vec;
        }()),
        is_mmd(ndists != total_seats),
        sampling_space(sampling_space) {
    // check the sizes are ok
    if (ndists - 1 > MAX_SUPPORTED_NUM_DISTRICTS ||
        total_seats - 1 > MAX_SUPPORTED_NUM_DISTRICTS) {
        throw std::overflow_error(
            "Too many districts. The maximum supported number of districts is " +
            std::to_string(MAX_SUPPORTED_NUM_DISTRICTS) + "."
        );
    }
    if (num_counties - 1 > MAX_SUPPORTED_NUM_COUNTIES) {
        throw std::overflow_error(
            "Too many counties. The maximum supported number of counties is " +
            std::to_string(MAX_SUPPORTED_NUM_COUNTIES) + "."
        );
    }
    if (V - 1 > MAX_SUPPORTED_NUM_VERTICES) {
        throw std::overflow_error(
            "Too many vertices. The maximum supported number of vertices is " +
            std::to_string(MAX_SUPPORTED_NUM_VERTICES) + "."
        );
    }
    if (V <= 1) {
        throw std::underflow_error(
            "Too few vertices. There must be at least 2 vertices."
        );
    }
};