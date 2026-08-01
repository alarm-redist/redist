#include <cstdio>
#include <limits>
#include <RcppArmadillo.h>
#include "wilson.h"

#include "splitting_schedule_types.h"
#include "random.h"
#include "base_plan_type.h"
#include "tree_splitting.h"

constexpr bool DEBUG_WILSON_VERBOSE = false; // Compile-time constant

namespace{



// debug function 
template <typename T>
std::string vec_to_string(std::vector<T> const &x, std::string_view name, int const len) {
    std::ostringstream oss;

    oss << name << " = c(";
    for (std::size_t i = 0; i < len; ++i) {
        oss << x[i];
        if (i + 1 < len) {
            oss << ", ";
        }
    }
    oss << ")\n";

    return oss.str();
};

/*
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
inline int rnbor(Graph const &g, int const vtx, RNGState &rng_state) {
    auto const &nbors = g[vtx];
    uint32_t const n_nbors = static_cast<uint32_t>(nbors.size());
    return nbors[rng_state.r_int(n_nbors)];
}

// flat graph version
inline int rnbor(FlatGraph const &g, int const vtx, RNGState &rng_state) {
    uint32_t const deg = static_cast<uint32_t>(g.degree(vtx));
    return g.neighbor(vtx, rng_state.r_int(deg));
}

// finds a vertex not visited yet by starting from
// lower until one not visited is hit 
int find_unvisited_vertex(
    int const V,
    const std::vector<bool> &visited, 
    int &lower) {
    for (int i = lower; i < V; i++) {
        if (!visited[i]){
            lower = i;
            return i;
        }
    }
    return -1;
}

int find_unvisited_county_vertex(
    std::vector<int> const &county_vertices,
    std::vector<bool> const &visited,
    std::size_t &lower_index
) {
    while (
        lower_index < county_vertices.size() &&
        visited[county_vertices[lower_index]]
    ) {
        ++lower_index;
    }

    return lower_index < county_vertices.size()
        ? county_vertices[lower_index]
        : -1;
}

// First sees if the start vertex has any unvisited neighors 
// then finds the smallest unvisited vertex starting from lower
int find_unvisited_walk_vertex(
    Graph const &admin_restricted_graph,
    int const V,
    const std::vector<bool> &visited, 
    int const last_start,
    int &lower) {
    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (
            static_cast<int>(admin_restricted_graph.size()) != V ||
            static_cast<int>(visited.size()) != V
        ) {
            throw std::runtime_error(
                "find_unvisited_walk_vertex: size mismatch."
            );
        }

        if (last_start < -1 || last_start >= V) {
            throw std::runtime_error(
                "find_unvisited_walk_vertex: invalid last_start."
            );
        }
    }

    // first see if last_start has any unvisited neighbors 
    // since admin restricted we don't need to worry about finding 
    // a neighbor across an admin boundary
    for (auto const u: admin_restricted_graph[last_start]){
        if (!visited[u]) return u;
    }

    for (int i = lower; i < V; i++) {
        if (!visited[i]){
            lower = i;
            return i;
        }
    }
    return -1;
}


/*
 * Add a deterministic tree on the active portion of one county.
 *
 * The root has already been inserted into the map-level tree by the
 * administrative Wilson walk. Ignored vertices are already marked visited,
 * so the DFS automatically stays inside the active county portion.
 *
 * Returns the number of newly added vertices. The root is not counted.
 */
static int add_county_to_tree_dfs(
    Graph const &county_restricted_g,
    int const county_id,
    int const root,
    std::vector<bool> &visited,
    DummyTreeQueue &queue,
    Tree &ust
) {
    int const V = static_cast<int>(county_restricted_g.size());

    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (root < 0 ||
            root >= static_cast<int>(visited.size())) {
            throw std::runtime_error(
                "add_county_to_tree_dfs: invalid county root."
            );
        }

        // if (
        //     static_cast<int>(
        //         map_params.counties[root]
        //     ) - 1 != county_id
        // ) {
        //     throw std::runtime_error(
        //         "add_county_to_tree_dfs: root belongs to wrong county."
        //     );
        // }

        if (!visited[root]) {
            throw std::runtime_error(
                "add_county_to_tree_dfs: county root is not already visited."
            );
        }
    }

    if constexpr (perf_config::bounds_checking) {
        if (root < 0 || root >= V) {
            std::ostringstream oss;
            oss << "add_county_to_tree_dfs received invalid root.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "root=" << root << "\n";
            oss << "V=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (!visited[root]) {
            std::ostringstream oss;
            oss << "add_county_to_tree_dfs root is not already visited.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "root=" << root << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    queue.clear();
    queue.push(root);

    int added = 0;

    while (!queue.empty()) {
        int const v = queue.pop();

        if constexpr (perf_config::bounds_checking) {
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "add_county_to_tree_dfs popped invalid vertex.\n";
                oss << "county_id=" << county_id << "\n";
                oss << "v=" << v << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        for (int const u : county_restricted_g[v]) {
            if constexpr (perf_config::bounds_checking) {
                if (u < 0 || u >= V) {
                    std::ostringstream oss;
                    oss << "add_county_to_tree_dfs saw invalid neighbor.\n";
                    oss << "county_id=" << county_id << "\n";
                    oss << "v=" << v << "\n";
                    oss << "u=" << u << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            /*
             * This skips:
             *
             * 1. the county root;
             * 2. vertices already added by this DFS;
             * 3. inactive vertices, because preparation marked them visited.
             */
            if (visited[u]) {
                continue;
            }

            ust[v].push_back(u);
            visited[u] = true;
            queue.push(u);

            ++added;
        }
    }

    return added;
}



// simple struct for returning sample UST information without doing 
// pass by reference for stuff like root
struct SampleSubUSTResult {
    int code;
    int root;
};


void OLD_TO_UPDATE_get_tree_pops_below(const Tree &ust, const int root, TreePopStack &stack,
                         const std::vector<unsigned int> &pop, std::vector<int> &pop_below) {
    stack.clear();
    // add the root
    // we don't care about parent here
    stack.push({root, 0, false});

    while (!stack.empty()) {
        auto [vtx, parent, is_revisiting] = stack.pop();

        // if visiting again then that means we've seen all the children
        if (is_revisiting) {
            int total_pop = pop[vtx];
            for (int child : ust[vtx]) {
                total_pop += pop_below[child];
            }
            pop_below[vtx] = total_pop;
        } else {
            // else push again and push the children
            stack.push({vtx, 0, true});
            for (const auto &child : ust[vtx]) {
                stack.push({child, 0, false});
            }
        }
    }

    return;
}


/*
 * Perform one Wilson loop-erased random walk and add the resulting path to
 * an undirected tree.
 *
 */
int add_walk_from(
    FlatGraph const &graph,
    // Graph const &graph,
    Tree &tree,
    int const start_vertex,
    int const max_tries,
    std::vector<int> &next_vertex,
    std::vector<bool> &visited,
    std::vector<bool> const &ignore,
    RNGState &rng_state
) {
    int const V = graph.size();

    /*
     * Optional defensive checks.
     *
     * These should either be behind your perf_config checks or removed in the
     * hottest build. They are useful while validating the refactor.
     */
    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (start_vertex < 0 || start_vertex >= V) {
            throw std::runtime_error("add_walk_from: invalid start_vertex.");
        }

        if (static_cast<int>(tree.size()) != V) {
            throw std::runtime_error("add_walk_from: tree has wrong size.");
        }

        if (static_cast<int>(visited.size()) != V) {
            throw std::runtime_error("add_walk_from: visited has wrong size.");
        }
        if (visited[start_vertex]) {
            throw std::runtime_error(
                "add_walk_from: start_vertex is already in the tree."
            );
        }
        if (static_cast<int>(next_vertex.size()) != V) {
            throw std::runtime_error(
                "add_walk_from: next_vertex has wrong size."
            );
        }

        if (static_cast<int>(ignore.size()) != V) {
            throw std::runtime_error(
                "add_walk_from: ignore has wrong size."
            );
        }
        if (max_tries <= 0) {
            throw std::runtime_error("add_walk_from: max_tries must be positive.");
        }
    }


    int curr_vertex = start_vertex;
    bool hit_tree = false;

    for (int tries = 0; tries < max_tries; ++tries) {
        int const curr_nbor =
            rnbor(graph, curr_vertex, rng_state);

        if (ignore[curr_nbor]) {
            continue;
        }

        /*
         * Record the latest exit from curr_vertex.
         * Revisiting curr_vertex overwrites its previous exit,
         * automatically removing loops from the final chain.
         */
        next_vertex[curr_vertex] = curr_nbor;

        if (visited[curr_nbor]) {
            hit_tree = true;
            break;
        }

        curr_vertex = curr_nbor;
    }

    if (!hit_tree) {
        return 0;
    }

    int added = 0;
    curr_vertex = start_vertex;

    if constexpr (DEBUG_WILSON_VERBOSE){
        std::cerr << "Adding following path to tree:";
    }
    /*
     * Follow the final successor chain to the existing tree.
     */
    while (!visited[curr_vertex]) {
        if constexpr (
            perf_config::supposedly_safe_input_checks
        ) {
            if (added >= V) {
                throw std::runtime_error(
                    "add_walk_from: successor chain contains "
                    "a cycle."
                );
            }
        }
        int const parent =
            next_vertex[curr_vertex];

        if constexpr (perf_config::bounds_checking) {
            if (parent < 0 || parent >= V) {
                throw std::runtime_error(
                    "add_walk_from: invalid successor vertex."
                );
            }
        }
        /*
         * Directed away from the existing root.
         */
        tree[parent].push_back(curr_vertex);
        if constexpr (DEBUG_WILSON_VERBOSE){
            std::cerr << "(" << parent << ","  << curr_vertex << "), ";
        }

        visited[curr_vertex] = true;
        curr_vertex = parent;
        ++added;
    }

    if constexpr (DEBUG_WILSON_VERBOSE){
        std::cerr << std::endl;
    }

    return added;
}


/*
 * Perform one Wilson random walk on an administrative-unit multigraph and
 * add its loop-erased successor chain to both:
 *
 *   1. admin_tree: the directed administrative-unit tree;
 *   2. map_tree:   the directed tree on original map vertices.
 *
 * Multigraph edge format:
 *
 *   edge[0] = neighboring administrative unit
 *   edge[1] = map vertex in the current administrative unit
 *   edge[2] = map vertex in the neighboring administrative unit
 *
 * Tree orientation:
 *
 * Both trees point away from the original Wilson root. If the random walk
 * selects:
 *
 *     child_admin -> parent_admin
 *
 * then the committed tree edges are:
 *
 *     admin_tree[parent_admin].push_back(child_admin);
 *     map_tree[parent_map_vertex].push_back(child_map_vertex);
 *
 * Return value:
 *
 *   > 0 = number of newly added administrative units
 *     0 = failed to hit the existing administrative tree within max_tries
 */
int add_walk_from_admin(
    Multigraph const &multigraph,
    Tree &admin_tree,
    Tree &map_tree,
    int const start_admin,
    int const max_tries,
    std::vector<AdminEdge> &next_admin_edge,
    std::vector<int> &admin_roots,
    std::vector<bool> &admin_visited,
    std::vector<bool> &map_visited,
    std::vector<bool> const &ignore,
    RNGState &rng_state
) {
    int const num_admin =
        static_cast<int>(multigraph.size());

    int const map_V =
        static_cast<int>(map_tree.size());

    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (start_admin < 0 ||
            start_admin >= num_admin) {
            throw std::runtime_error(
                "add_walk_from_admin: invalid start_admin."
            );
        }

        if (static_cast<int>(admin_tree.size()) !=
            num_admin) {
            throw std::runtime_error(
                "add_walk_from_admin: admin_tree has wrong size."
            );
        }

        if (static_cast<int>(next_admin_edge.size()) !=
            num_admin) {
            throw std::runtime_error(
                "add_walk_from_admin: next_admin_edge has wrong size."
            );
        }

        if (static_cast<int>(admin_visited.size()) !=
            num_admin) {
            throw std::runtime_error(
                "add_walk_from_admin: admin_visited has wrong size."
            );
        }

        if (static_cast<int>(map_visited.size()) !=
            map_V) {
            throw std::runtime_error(
                "add_walk_from_admin: map_visited has wrong size."
            );
        }

        if (static_cast<int>(ignore.size()) !=
            map_V) {
            throw std::runtime_error(
                "add_walk_from_admin: ignore has wrong size."
            );
        }

        if (admin_visited[start_admin]) {
            throw std::runtime_error(
                "add_walk_from_admin: start_admin is already in tree."
            );
        }
        if (
            static_cast<int>(admin_roots.size()) !=
            num_admin
        ) {
            throw std::runtime_error(
                "add_walk_from_admin: admin_roots has wrong size."
            );
        }
        if (max_tries <= 0) {
            throw std::runtime_error(
                "add_walk_from_admin: max_tries must be positive."
            );
        }
    }

    int current_admin = start_admin;
    bool hit_admin_tree = false;

    /*
     * Perform the raw random walk.
     *
     * Every time current_admin is revisited, its previous successor edge is
     * overwritten. Consequently, loops are removed implicitly from the final
     * successor chain.
     */
    for (int tries = 0;
         tries < max_tries;
         ++tries) {
        // Get all the edges from this admin unit
        auto const &incident_edges =
            multigraph[current_admin];

        // optional check for if there are no active edges to walk along
        if constexpr (perf_config::supposedly_safe_input_checks) {
            int num_active_edges = 0;

            for (AdminEdge const &edge : incident_edges) {
                if (
                    !ignore[edge.current_map_vertex] &&
                    !ignore[edge.neighbor_map_vertex]
                ) {
                    ++num_active_edges;
                }
            }

            if (num_active_edges == 0) {
                std::ostringstream oss;
                oss << "add_walk_from_admin: active administrative unit has "
                    << "no active boundary edges.\n";
                oss << "current_admin=" << current_admin << "\n";
                oss << "start_admin=" << start_admin << "\n";
                oss << "incident_edges.size()="
                    << incident_edges.size() << "\n";

                throw std::runtime_error(oss.str());
            }
        }

        if (incident_edges.empty()) {
            /*
             * This administrative unit cannot reach the growing tree in this
             * multigraph.
             */
            return 0;
        }

        // randomly sample an edge 
        int const selected_index =
                rng_state.r_int(
                    static_cast<std::uint32_t>(
                        incident_edges.size()
                    )
                );

        // get the edge sampled 
        AdminEdge const &selected_edge =
            incident_edges[selected_index];

        int const proposed_admin = selected_edge.neighbor_admin;

        if constexpr (perf_config::bounds_checking) {
            if (proposed_admin < 0 ||
                proposed_admin >= num_admin) {
                throw std::runtime_error(
                    "add_walk_from_admin: invalid neighboring admin."
                );
            }

            if (selected_edge.current_map_vertex < 0 ||
                selected_edge.current_map_vertex >= map_V ||
                selected_edge.neighbor_map_vertex < 0 ||
                selected_edge.neighbor_map_vertex >= map_V) {
                throw std::runtime_error(
                    "add_walk_from_admin: invalid underlying map edge."
                );
            }
        }

        /*
         * The stored multigraph is currently defined over the complete map.
         * Reject boundary edges whose original-map endpoints are outside the
         * active region.
         * 
         * We do not need to seperately check that we are ignoring this county 
         * as a county is marked to ignore iff all vertices in it are marked to 
         * ignore so its safe.
         *
         * A rejected edge does not move the random walk.
         */
        if (
            ignore[
                selected_edge.neighbor_map_vertex
            ] ||
            ignore[
                selected_edge.current_map_vertex
            ]
        ) {
            continue;
        }

        /*
         * Record the most recent accepted exit from current_admin.
         *
         * If current_admin is revisited later, this assignment overwrites the
         * previous exit and thereby breaks the old loop.
         */
        next_admin_edge[current_admin] = selected_edge;

        if (admin_visited[proposed_admin]) {
            hit_admin_tree = true;
            break;
        }

        current_admin = proposed_admin;
    }

    if (!hit_admin_tree) {
        return 0;
    }

    /*
     * Follow the final successor chain from start_admin to the preexisting
     * administrative tree.
     */
    int added = 0;
    current_admin = start_admin;

    while (!admin_visited[current_admin]) {
        /*
         * There can be at most num_admin unvisited vertices in a valid chain.
         * This check catches a malformed successor cycle rather than hanging.
         */
        if constexpr (perf_config::supposedly_safe_input_checks) {
            if (added >= num_admin) {
                throw std::runtime_error(
                    "add_walk_from_admin: successor chain contains a cycle."
                );
            }
        }

        AdminEdge const &selected_edge =
            next_admin_edge[current_admin];

        int const parent_admin = selected_edge.neighbor_admin;

        int const child_map_vertex = selected_edge.current_map_vertex;

        int const parent_map_vertex = selected_edge.neighbor_map_vertex;

        if constexpr (perf_config::bounds_checking) {
            if (parent_admin < 0 ||
                parent_admin >= num_admin) {
                throw std::runtime_error(
                    "add_walk_from_admin: invalid successor admin."
                );
            }

            if (child_map_vertex < 0 ||
                child_map_vertex >= map_V ||
                parent_map_vertex < 0 ||
                parent_map_vertex >= map_V) {
                throw std::runtime_error(
                    "add_walk_from_admin: invalid successor map edge."
                );
            }
        }

        /*
         * Add the administrative edge, oriented away from the original
         * administrative root.
         */
        admin_tree[parent_admin].push_back(
            current_admin
        );

        /*
         * Add the corresponding exact map edge, oriented away from the
         * original map-level root.
         */
        map_tree[parent_map_vertex].push_back(
            child_map_vertex
        );

        /*
         * Every nonroot administrative unit receives exactly one map-level
         * root: the endpoint of its edge toward its administrative parent.
         *
         * The later within-county Wilson call grows toward this vertex.
         */
        map_visited[child_map_vertex] = true;
        // Mark this vertex as the one in this unit 
        // Notice this is a directed edge (parent_map_vertex, child_map_vertex)
        // and (parent_admin_unit, child_admin_unit)
        admin_roots[current_admin] = child_map_vertex;

        /*
         * The current administrative unit has now joined the tree.
         */
        admin_visited[current_admin] = true;

        current_admin = parent_admin;
        ++added;
    }

    return added;
}


// We assume that ignore has been properly set and nothing else
// has been cleared 
template <bool SkipUnsplittableTrees>
SampleSubUSTResult sample_full_ust(
    MapParams const &map_params, 
    Tree &tree, // FlatGraph &tree, 
    double const lower, double const upper, 
    std::vector<bool> &visited, const std::vector<bool> &ignore, 
    Tree &county_tree, 
    WilsonGraphScratch &g_scratch,
    WilsonMultiGraphScratch &mg_scratch,
    RNGState &rng_state
    ) {
    // auto t1_start = std::chrono::steady_clock::now();
    int const V = map_params.V;
    int const n_county = map_params.num_counties;

    // pick root
    int const root = find_unvisited_vertex(V, visited, g_scratch.smallest_v_seen);
    if (root < 0) {
        throw std::runtime_error(
            "sample_full_ust: no unvisited root exists."
        );
    }
    if constexpr(DEBUG_WILSON_VERBOSE){
        std::cerr << "Starting with root of " << root << std::endl;
    }
    visited[root] = true;
    g_scratch.remaining--;

    int root_county = map_params.counties[root] - 1;
    mg_scratch.c_visited[root_county] = true;
    mg_scratch.c_remaining--;
    mg_scratch.admin_roots[root_county] = root;

    

    constexpr int admin_tries_multiple = 500;
    // optional check for overflow
    if constexpr (perf_config::bounds_checking) {
        std::int64_t const max_try_wide =
            static_cast<std::int64_t>(n_county) *
            static_cast<std::int64_t>(admin_tries_multiple);
        if (max_try_wide > std::numeric_limits<int>::max()) {
            std::ostringstream oss;
            oss << "sample_full_ust: max_try exceeds int maximum.\n";
            oss << "remaining=" << g_scratch.remaining << "\n";
            oss << "max_try=" << max_try_wide << "\n";

            throw std::overflow_error(oss.str());
        }
    }

    int const max_admin_tries = admin_tries_multiple * n_county;

    // Connect counties
    while (mg_scratch.c_remaining > 0) {
        
        int unvisited_county = find_unvisited_vertex(
            mg_scratch.num_admin_units, 
            mg_scratch.c_visited, 
            mg_scratch.smallest_county_seen    
        );
        if (unvisited_county < 0) {
            std::ostringstream oss;
            oss << "sample_full_ust: c_remaining="
                << mg_scratch.c_remaining
                << " but no unvisited administrative unit exists.";

            throw std::runtime_error(oss.str());
        }
        // update visited list and constructed tree
        int admins_added = add_walk_from_admin(
            map_params.cg,
            county_tree,
            tree,
            unvisited_county,
            max_admin_tries,
            mg_scratch.next_admin_edge,
            mg_scratch.admin_roots,
            mg_scratch.c_visited,
            visited, ignore, rng_state
        );
        if (admins_added == 0) { // bail
            return SampleSubUSTResult{1, root};
        }
        /*
        * Each newly added county contributes:
        *
        *   - one fewer unvisited county;
        *   - one newly visited map vertex serving as that county's root.
        */
        mg_scratch.c_remaining -= admins_added;
        g_scratch.remaining -= admins_added;
    }

    // optional toggle of skipping unsplittable trees 
    if constexpr (SkipUnsplittableTrees) {
    // figure out which counties will not need to be split
    if (n_county > 1) {
        // don't need to fill pop below since it gets reset
        OLD_TO_UPDATE_get_tree_pops_below(
            county_tree, map_params.counties[root] - 1, 
            mg_scratch.county_stack, mg_scratch.county_pop,
            mg_scratch.cty_pop_below);
        for (int county_i = 0; county_i < n_county; county_i++) {
            if (mg_scratch.admin_ignore[county_i]) continue;
            // check child counties
            int children = county_tree[county_i].size();
            int split_ub = mg_scratch.cty_pop_below[county_i];
            int split_lb = split_ub - mg_scratch.county_pop[county_i];
            if (lower - 1 < mg_scratch.county_pop[county_i])
                split_lb = (int)lower;
            for (int j = 0; j < children; j++) {
                int pop_child = mg_scratch.cty_pop_below[county_tree[county_i][j]];
                if (pop_child >= 0 && pop_child < split_lb) {
                    split_lb = pop_child;
                }
            }

            // split_lb < split_ub so the smallest possible population is
            // min(split_lb, tot_pop - split_ub)
            // its impossible to split if smallest possible size is bigger than largest ub
            // bool miss_first = std::min(split_lb, tot_pop - split_ub) > upper;
            // // biggest possible population is max(split_ub, total_pop - split_lb)
            // // its impossible to split if largest possible size is smaller than smallest lb
            // bool miss_second = std::max(split_ub, tot_pop - split_lb) < lower;

            bool miss_first = split_ub < lower || split_lb > upper;
            bool miss_second = (mg_scratch.total_pop - split_lb) < lower || (mg_scratch.total_pop - split_ub) > upper;

            // impossible for this county to need to be split
            // so we fill in a dummy tree 
            if (mg_scratch.cty_pop_below[county_i] >= 0 && (miss_first && miss_second)) {
                int const vertices_added = add_county_to_tree_dfs(
                    map_params.county_restricted_graph,
                    county_i, // pass in the county CAREFUL COUNTIES ARE 1-indexed
                    mg_scratch.admin_roots[county_i],
                    visited,
                    g_scratch.dummy_county_tree_queue,
                    tree
                );
                // note we filled this one in deterministically assuming not a singleton county
                if (vertices_added > 0) {
                    mg_scratch.deterministic_counties.push_back(
                        county_i
                    );
                }
                g_scratch.remaining -= vertices_added; // function does not count the already visited county root
            }
        }
    }
    }

    if constexpr(DEBUG_WILSON_VERBOSE){
        std::cerr << "Starting to draw within county trees!\n";
    }
    
    // Generate tree within each county
    if (g_scratch.remaining > 0) {
        constexpr int graph_tries_multiple = 500;
        // optional check for overflow
        if constexpr (perf_config::bounds_checking) {
            std::int64_t const max_try_wide =
                static_cast<std::int64_t>(graph_tries_multiple) *
                static_cast<std::int64_t>(g_scratch.remaining) *
                (static_cast<int>(std::log(g_scratch.remaining)) + 1);
            if (max_try_wide > std::numeric_limits<int>::max()) {
                std::ostringstream oss;
                oss << "sample_full_ust: max_try exceeds int maximum.\n";
                oss << "remaining=" << g_scratch.remaining << "\n";
                oss << "max_try=" << max_try_wide << "\n";

                throw std::overflow_error(oss.str());
            }
        }

        int max_try = graph_tries_multiple * g_scratch.remaining * (static_cast<int>(std::log(g_scratch.remaining)) + 1);

        // int unvisited_vertex = root;
        while (g_scratch.remaining > 0) {
            int unvisited_vertex = find_unvisited_vertex(V, visited, g_scratch.smallest_v_seen);
            // Profiled on NY and the comment out method seemed a little slower
            // unvisited_vertex = find_unvisited_walk_vertex(map_params.county_restricted_graph,
            //     V, visited, unvisited_vertex,  g_scratch.smallest_v_seen);

            if (unvisited_vertex < 0) {
                std::ostringstream oss;
                oss << "sample_full_ust: unvisited_vertex="
                    << g_scratch.remaining
                    << " but no unvisited vertex exists.";

                throw std::runtime_error(oss.str());
            }
            // random walk from `unvisited_vertex` until we hit the path
            int vertices_added = add_walk_from(
                map_params.county_restricted_flat_graph,
                tree,
                unvisited_vertex,
                max_try,
                g_scratch.next_vertex,
                visited, ignore, rng_state
            );
            // update visited list and constructed tree
            if (vertices_added == 0) { // bail
                return SampleSubUSTResult{1, root};
            }
            g_scratch.remaining -= vertices_added; // no correction needed because walk until doesn't count the hit
            // vertex already in the tree 
        }
    }

    return SampleSubUSTResult{0, root};
}

}


// sample_non_hierarchical_ust


/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2025/3
 * Purpose: Encapsulation of uniform spanning tree sampler functions
 ********************************************************/


 // The idea is to move all this outside of sample sub ust

// Bigger picture refactor 
// - Pick the root in the USTSampler call, not wilson
// - Since we prep the graph anyways just make the root the smallest non ignore variable seen 
// - Non hierarchical wilson
// - Add draw_fresh_ust and fill_in_ust 


// TODO
// - Make a helper `sample_subgraph_ust` function 
// - Make a helper `sample_submultigraph_ust` function 
// - `sample_ust` should involve calls to these functions 
// - Add something tracking counties skipped and a vertex in the county


// We assume that ignore has been correctly set, nothing else
// INCLUDING THE TREE HAS NOT BEEN CLEARED!


template <typename IsActive>
int USTSampler::prep_fresh_ust_call(IsActive const &is_active){
    // We will now 
    // - properly set the visited vector 

    // - TODO: set up county multigraph stuff 

    // First we clear the multigraph stuff
    mg_scratch.c_remaining = 0;
    mg_scratch.total_pop = 0;
    // Set 
    std::fill(
        mg_scratch.c_visited.begin(),
        mg_scratch.c_visited.end(),
        true
    );
    std::fill(
        mg_scratch.admin_ignore.begin(),
        mg_scratch.admin_ignore.end(),
        true
    );
    mg_scratch.deterministic_counties.clear();
    

    bool first_county_seen = true;

    g_scratch.remaining = 0;
    int const V = map_params.V;
    for (int v = 0; v < V; v++)
    {
        // function tells us if this vertex is in the subgraph
        // we're drawing the tree on
        bool const active = is_active(v);
        ignore[v] = !active;
        if (!active) {
            // if we're ignoring we just mark as visited and continue
            visited[v] = true;
            // there's nothing more to do 
            continue;
        } 
        // Else it means v is part of the subgraph we care about 
        visited[v] = false;
        g_scratch.remaining++;
        auto v_county = map_params.counties[v] - 1;
        if (mg_scratch.c_visited[v_county]) {
            mg_scratch.c_visited[v_county] = false;
            mg_scratch.admin_ignore[v_county] = false;
            ++mg_scratch.c_remaining;
            // clear the infor associated with this county
            mg_scratch.county_pop[v_county] = 0.0;
            // clear the tree
            county_tree[v_county].clear();

            // This will be our initial start vertex since this will 
            // TODO add another check 
            if(first_county_seen){
                g_scratch.smallest_v_seen = v;
                first_county_seen = false;
                mg_scratch.smallest_county_seen = v_county;
            }
            // check if this is a smaller county than we've seen before
            if(v_county < mg_scratch.smallest_county_seen){
                mg_scratch.smallest_county_seen = v_county;
            }

            // TODO add this vertex to county root vector 
            // TODO in the future clear the county multigraph and add this edge 
        }
        mg_scratch.total_pop += map_params.pop[v];
        mg_scratch.county_pop[v_county] += map_params.pop[v];
        // Now we clear the ust at this vertex
        // ust.clear_vertex(i);
        ust[v].clear();
        // Now we clear the subgraph 
        // wilson_submap.clear_vertex(v);
        // // now we add all the neighbors not ignored and in the same county 
        // for (auto const u: map_params.g[v])
        // {
        //     if (!ignore[u] && map_params.counties[u] - 1 == v_county){
        //         wilson_submap.add_directed_edge(v, u);
        //     }
        // }
        
    }

    return g_scratch.remaining;
}


template <bool SkipUnsplittableTrees, typename IsActive>
USTDrawResult USTSampler::draw_fresh_ust(
    double const lower,
    double const upper,
    RNGState &rng_state,
    IsActive const &is_active
) {
    int const num_tree_vertices =
        prep_fresh_ust_call(is_active);

    auto const result =
        sample_full_ust<SkipUnsplittableTrees>(
            map_params,
            ust,
            lower,
            upper,
            visited,
            ignore,
            county_tree,
            g_scratch,
            mg_scratch,
            rng_state
        );

    if constexpr (
        perf_config::object_integrity_checking
    ) {
        if (result.code == 0) {
            check_tree_integrity(
                get_vertex_tree(),
                "Just called sample_full_ust in draw_fresh_ust.\n",
                result.root,
                num_tree_vertices,
                true
            );
        }
    }

    return USTDrawResult{
        result.code == 0,
        num_tree_vertices,
        result.root
    };
}

// USTDrawResult USTSampler::draw_fresh_ust(
//       double const lower, double const upper,
//         RNGState &rng_state) {


//     // prep the inputs and get the number of vertices 
//     int const num_tree_vertices = prep_fresh_ust_call(
//         [&](int const v) noexcept {
//             return plan.region_ids[v] ==
//                    region_to_draw_tree_on;
//         }
//     );
//     // REprintf("%d Vertices!\n", num_tree_vertices);

//     auto const result = sample_full_ust<true>(map_params, ust, 
//         lower, upper, visited, ignore,
//         county_tree, g_scratch, mg_scratch, rng_state
//     );

//     if constexpr(perf_config::object_integrity_checking){
//         if (result.code == 0){
//             check_tree_integrity(
//                 get_vertex_tree(),
//                 "Just called `sample_sub_ust` in attempt_to_draw_tree_on_region\n",
//                 result.root,
//                 num_tree_vertices,
//                 true
//             );
//         }
//     }

//     // result.code == 0 means it was successful
//     return USTDrawResult{result.code == 0, num_tree_vertices, result.root};
// }



USTDrawResult USTSampler::attempt_to_draw_tree_on_region(
    RNGState &rng_state, Plan const &plan, const int region_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;
    double sample_sub_ust_lower, sample_sub_ust_upper;

    if (use_custom_bounds){
        // if using custom upper and lwoer bounds set those
        sample_sub_ust_lower = custom_sample_sub_ust_lower;
        sample_sub_ust_upper = custom_sample_sub_ust_upper;
    }else{
        // else just do upper and lower scaled by min and max possible cut size
        // get upper and lower bounds on region pops
        auto min_max_pair = splitting_schedule.all_regions_min_and_max_possible_cut_sizes
                                [plan.region_sizes[region_to_draw_tree_on]];
        sample_sub_ust_lower = min_max_pair.first * map_params.lower; // lower
        sample_sub_ust_upper = min_max_pair.second * map_params.upper; // upper
    }
    // Now sample a uniform spanning tree drawn on that region
    return draw_fresh_ust<true>(
        sample_sub_ust_lower,
        sample_sub_ust_upper,
        rng_state,
        [&plan, region_to_draw_tree_on](
            int const v
        ) noexcept {
            return plan.region_ids[v] ==
                region_to_draw_tree_on;
        }
    );
}


USTDrawResult USTSampler::attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
                                                       const int region1_to_draw_tree_on,
                                                       const int region2_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;

    double sample_sub_ust_lower, sample_sub_ust_upper;

    if (use_custom_bounds){
        // if using custom upper and lwoer bounds set those
        sample_sub_ust_lower = custom_sample_sub_ust_lower;
        sample_sub_ust_upper = custom_sample_sub_ust_upper;
    }else{
        // else just do upper and lower scaled by min and max possible cut size
        // get upper and lower bounds on region pops
        int merged_region_size =
            plan.region_sizes[region1_to_draw_tree_on] + plan.region_sizes[region2_to_draw_tree_on];
        // get upper and lower bounds on region pops
        auto min_max_pair =
            splitting_schedule.all_regions_min_and_max_possible_cut_sizes[merged_region_size];

        sample_sub_ust_lower = min_max_pair.first * map_params.lower; // lower
        sample_sub_ust_upper = min_max_pair.second * map_params.upper; // upper
    }


    // Now sample a uniform spanning tree drawn on that region
    return draw_fresh_ust<true>(
        sample_sub_ust_lower,
        sample_sub_ust_upper,
        rng_state,
        [
            &plan,
            region1_to_draw_tree_on,
            region2_to_draw_tree_on
        ](int const v) noexcept {
            RegionID const region =
                plan.region_ids[v];

            return
                region == region1_to_draw_tree_on ||
                region == region2_to_draw_tree_on;
        }
    );
}

bool USTSampler::fill_in_skipped_subtrees(
    EdgeBitset &packed_forest_edges,
    RNGState &rng_state, 
    int const max_tries_multiple
){
    if (mg_scratch.deterministic_counties.empty()) {
        return true;
    }

    if constexpr (perf_config::supposedly_safe_input_checks) {
        if (max_tries_multiple <= 0) {
            throw std::runtime_error(
                "fill_in_skipped_subtrees: max_tries must be positive."
            );
        }
    }

    int const V = map_params.V;

    // We loop over each of the deterministically filled counties 
    /*
     * Remove every deterministic within-county tree.
     *
     * Cross-county children are retained. The edge entering the county is
     * stored in the parent county's adjacency list, so it is not touched here.
     */
    for (int const county_i: mg_scratch.deterministic_counties){
        int county_remaining = 0;
        std::size_t lower_index = 0;
        if constexpr (perf_config::bounds_checking) {
            if (
                county_i < 0 ||
                county_i >= map_params.num_counties
            ) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: invalid county."
                );
            }
            int county_root = mg_scratch.admin_roots[county_i];
            if (
                county_root < 0 ||
                county_root >= V
            ) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: invalid county root."
                );
            }
            if (
                static_cast<int>(
                    map_params.counties[county_root]
                ) - 1 != county_i
            ) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: root belongs to wrong county."
                );
            }
            if (ignore[county_root]) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: county root is ignored."
                );
            }
        }
        // get the county root 
        int county_root = mg_scratch.admin_roots[county_i];
        visited[county_root] = true;
        // g_scratch.smallest_v_seen = county_root;
        // We keep the county root as visited 
        // Now we walk through `ust` and erase all within county edges 
        // We also set those edges in the packed forest to be false 
        auto &queue =
            g_scratch.dummy_county_tree_queue;

        queue.clear();
        queue.push(county_root);

        while (!queue.empty()) {
            int const v = queue.pop();
            auto &children = ust[v];

        /*
        * Compact the child vector in place:
        *
        * - within-county children are removed;
        * - cross-county children are retained.
        */
        // NOTE when moving to multilevel trees might need to think if we want to preserve edges
        // in which case 
        //  - define a std::size_t write_pos = 0; 
        //  - children[write_pos++] = child; when the child is cross county edge
        //  - children.resize(write_pos); at the end after going through all children 
        
            for (int const child : children) {
                int const child_county = static_cast<int>(map_params.counties[child]) - 1;

                if (child_county == county_i) {
                    // if (child < g_scratch.smallest_v_seen){
                    //     g_scratch.smallest_v_seen = child;
                    // }
                    /*
                     * Follow the old deterministic tree before deleting its
                     * edge, so we discover the entire county subtree.
                     */
                    queue.push(child);
                    /*
                     * Clear the edge in the packed forest 
                     */
                    packed_forest_edges.clear_edge(
                        v,
                        child,
                        map_params.graph_edge_index
                    );
                    /*
                     * Every nonroot county vertex must be resampled.
                     */
                    visited[child] = false;
                    ++county_remaining;
                }
            }
            // now clear all these edges since we're erasing all within county edges
            // and we don't care about across count edges 
            children.clear();
        }

        if (county_remaining == 0) {
            continue;
        }

        if constexpr (perf_config::bounds_checking) {
            std::int64_t const max_try_wide =
                static_cast<std::int64_t>(max_tries_multiple) *
                static_cast<std::int64_t>(county_remaining) *
                (static_cast<int>(std::log(county_remaining)) + 1);
            if (max_try_wide > std::numeric_limits<int>::max()) {
                std::ostringstream oss;
                oss << "sample_full_ust: max_try exceeds int maximum.\n";
                oss << "remaining=" << county_remaining << "\n";
                oss << "max_try=" << max_try_wide << "\n";

                throw std::overflow_error(oss.str());
            }
        }

        int max_tries = max_tries_multiple * county_remaining * (static_cast<int>(std::log(county_remaining)) + 1);

        while (county_remaining > 0){
            int unvisited_vertex = find_unvisited_county_vertex(
                    map_params.county_vertices[county_i],
                    visited,
                    lower_index
            );
            if (unvisited_vertex < 0) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: "
                    "county_remaining is positive but no "
                    "unvisited county vertex exists."
                );
            }
            // perform another loop erased walk 
            int const vertices_added = add_walk_from(
                map_params.county_restricted_flat_graph,
                ust,
                unvisited_vertex,
                max_tries,
                g_scratch.next_vertex,
                visited, ignore, rng_state
            );

            // if nothing added then we throw an error
            if (vertices_added == 0) {
                throw std::runtime_error(
                    "fill_in_skipped_subtrees: "
                    "Wilson failed while filling county."
                );
            }

            county_remaining -= vertices_added;

        }
        // Now starting from the root we can traverse the tree and 
        // fill in the packed forest 
        // TODO confirm the tree is still a properly directed tree
        queue.push(county_root);
        while (!queue.empty()) {
            int const v = queue.pop();
            for (int const child : ust[v]) {
                int const child_county = map_params.counties[child] - 1;

                if (child_county != county_i) {
                    continue;
                }
                // set the edge
                packed_forest_edges.set_edge(
                    v,
                    child,
                    map_params.graph_edge_index
                );
                queue.push(child);
            }
        }
    }
    mg_scratch.deterministic_counties.clear();
    return true;
}

std::pair<bool, EdgeCut> USTSampler::try_to_find_and_erase_splittable_edge(
    Plan const &plan, int const split_region1, int const split_region2, int const root,
    ScoringFunction const &scoring_function, RNGState &rng_state, TreeSplitter &tree_splitter,
    int const region_populations, int const region_size, bool const save_selection_prob) {
    // We assume a tree has already been successfully drawn so
    // try to find a valid cut
    auto cut_size_bounds =
        splitting_schedule.all_regions_min_and_max_possible_cut_sizes[region_size];
    int min_possible_cut_size = cut_size_bounds.first;
    int max_possible_cut_size = cut_size_bounds.second;

    // REprintf("Remainder Size: %d Cut sizes:", region_size);
    // for(auto const &v: splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size]){
    //     REprintf("%d, ", v);
    // }
    // REprintf("\n");

    std::pair<bool, EdgeCut> edge_search_result = tree_splitter.attempt_to_find_edge_to_cut(
        map_params, scoring_function, rng_state, plan, split_region1, split_region2, ust, root,
        pops_below_vertex, region_populations, region_size,
        min_possible_cut_size, max_possible_cut_size,
        splitting_schedule.all_regions_smaller_cut_sizes_to_try[region_size],
        save_selection_prob);

    bool search_successful = std::get<0>(edge_search_result);
    // return false if unsuccessful
    if (!search_successful)
        return std::make_pair(false, EdgeCut());

    // If successful extract the edge cut info
    EdgeCut cut_edge = std::get<1>(edge_search_result);
    // Now erase the cut edge in the tree
    // ust.erase_directed_edge(cut_edge);
    erase_tree_edge(ust, cut_edge);

    return edge_search_result;
}

std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_split(
    RNGState &rng_state, ScoringFunction const &scoring_function, TreeSplitter &tree_splitter,
    Plan const &plan, int const region_to_split, int const new_region_id,
    bool const save_selection_prob) {
    // Try to draw a tree
    auto const result = attempt_to_draw_tree_on_region(rng_state, plan, region_to_split);
    // return false if unsuccessful
    if (!result.successful)
        return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size = plan.region_sizes[region_to_split];
    int region_to_split_population = plan.region_pops[region_to_split];

    return try_to_find_and_erase_splittable_edge(plan, region_to_split, new_region_id, result.root,
        scoring_function,
                                         rng_state, tree_splitter, region_to_split_population,
                                         region_to_split_size, save_selection_prob);
}

std::pair<bool, EdgeCut> USTSampler::attempt_to_find_valid_tree_mergesplit(
    RNGState &rng_state, ScoringFunction const &scoring_function, TreeSplitter &tree_splitter,
    Plan const &plan, int const merge_region1, int const merge_region2,
    bool const save_selection_prob) {
    // Try to draw a tree
    auto const result = attempt_to_draw_tree_on_merged_region(rng_state, plan, merge_region1, merge_region2);
    // return false if unsuccessful
    if (!result.successful)
        return std::make_pair(false, EdgeCut());

    // Else try to find a valid cut
    int region_to_split_size =
        plan.region_sizes[merge_region1] + plan.region_sizes[merge_region2];
    int region_to_split_population =
        plan.region_pops[merge_region1] + plan.region_pops[merge_region2];

    return try_to_find_and_erase_splittable_edge(plan, merge_region1, merge_region2, result.root, scoring_function,
                                         rng_state, tree_splitter, region_to_split_population,
                                         region_to_split_size, save_selection_prob);
}


void USTSampler::check_tree_integrity(
      Tree const &a_ust,
      std::string_view where,
      int root,
      int expected_tree_vertices,
      bool check_vertex_count
    ){
    // check no garbage vertices in the tree 
    check_directed_tree_edges_are_graph_edges(
        a_ust, map_params.graph_edge_index,
        where
    );
    // now check the tree returned is actually a directed tree
    check_is_directed_tree(
        a_ust,
        where,
        root,
        expected_tree_vertices,
        check_vertex_count,
        visited,
        stack
    );
}

std::pair<bool, int>  USTSampler::draw_tree_on_subgraph(
      RNGState &rng_state, std::vector<bool> const &vertices_to_ignore,
      bool const skip_unsplittable_subtrees, 
      const double lower, const double upper,
      WilsonTimes &wilson_times
){
    auto prep_start_time = std::chrono::steady_clock::now();

    int V = map_params.V;
    // prepare the inputs for a `sample_ust` call
    // The function marks a vertex as ignore if true in `vertices_to_ignore` 
    int const num_vertices = prep_fresh_ust_call(
        [&vertices_to_ignore](
            int const v
        ) noexcept {
            return !vertices_to_ignore[v];
        }
    );

    auto prep_end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::ratio<1>> prep_time_diff = prep_end_time - prep_start_time;
    wilson_times.input_prep_time += prep_time_diff.count();

    auto sample_ust_start_time = std::chrono::steady_clock::now();
    // Now sample a uniform spanning tree drawn on that region

    std::pair<bool, int> result;

    if (skip_unsplittable_subtrees){
        auto const ust_result = sample_full_ust<true>(map_params, ust, 
            lower, upper, visited, ignore,
            county_tree, g_scratch, mg_scratch, rng_state
        );
        result = std::make_pair(ust_result.code == 0, ust_result.root);
    }else{
        auto const ust_result = sample_full_ust<false>(map_params, ust, 
            lower, upper, visited, ignore,
            county_tree, g_scratch, mg_scratch, rng_state
        );
        result = std::make_pair(ust_result.code == 0, ust_result.root);
    }
    auto sample_ust_end_time = std::chrono::steady_clock::now();
    std::chrono::duration<double, std::ratio<1>> ust_time_diff = sample_ust_end_time - sample_ust_start_time;
    wilson_times.sub_ust_call_time += ust_time_diff.count();
    
    bool const valid_tree = result.first;

    if constexpr(perf_config::object_integrity_checking){
        if (valid_tree){
            check_tree_integrity(
                get_vertex_tree(),
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_region\n",
                result.second,
                num_vertices,
                true
            );
        }
    }

    return result;
}

