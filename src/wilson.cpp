#include <cstdio>
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
    return lower;
}


/*
 * Builds a deterministic spanning tree on a county using depth first search
 */
static void add_county_to_tree_dfs(
    Graph const &county_restricted_g,
    int const county_id,
    std::vector<int> const &county_vertices, // vector of vertices in the county 
    std::vector<bool> &visited,
    DummyTreeQueue &queue,
    Tree &ust // FlatGraph &ust
) {
    int const V = static_cast<int>(county_restricted_g.size());
    int const n_vtx = static_cast<int>(county_vertices.size());

    // Optional checks 
    if constexpr(perf_config::supposedly_safe_input_checks){
        if (n_vtx <= 0) {
            return;
        }
    }

    // set the root as the first county vertex 
    int root = -1;

    // Walk through the vertices in the county 
    for (int const v : county_vertices) {
        // optional bounds checking 
        if constexpr (perf_config::bounds_checking){
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "add_county_dfs_tree_edges got invalid county vertex.\n";
                oss << "county_id=" << county_id << "\n";
                oss << "v=" << v << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        // Check if we have already visited this vertex in the Wilson call (sample_sub_ust)
        // If yes then that means this is the vertex to use as the root for this county's portion
        // of the tree. 
        if (visited[v]) {
            // Checks we haven't already found a root 
            if constexpr(perf_config::supposedly_safe_input_checks){
                // This shouldn't be neccesary if the tree is truly a directed tree
                if (root != -1) {
                    std::ostringstream oss;
                    oss << "add_county_dfs_tree_edges found multiple already-visited "
                        << "vertices in county.\n";
                    oss << "county_id=" << county_id << "\n";
                    oss << "first_root=" << root << "\n";
                    oss << "second_root=" << v << "\n";
                    throw std::runtime_error(oss.str());
                }
            }
            // we've found the vertex to use as the root so we can exit if not error checking
            root = v;
            // if not error checking then break out now since the tree is supposed to be
            // a directed tree 
            if constexpr(!perf_config::supposedly_safe_input_checks){
                break;
            }
        }

    }

    if constexpr(perf_config::bounds_checking){
        if (root < 0) {
            std::ostringstream oss;
            oss << "add_county_dfs_tree_edges could not find already-visited county root.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "county_vertices.size()=" << county_vertices.size() << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    // clear the stack and start from the root
    queue.clear();
    queue.push(root);

    int seen_count = 1;

    // perform DFS
    while (!queue.empty()) {
        int const v = queue.pop();

        if constexpr(perf_config::bounds_checking){
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "add_county_dfs_tree_edges popped invalid vertex.\n";
                oss << "county_id=" << county_id << "\n";
                oss << "v=" << v << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }


        for (auto const u : county_restricted_g[v]) {
            if constexpr(perf_config::bounds_checking){
                if (u < 0 || u >= V) {
                    std::ostringstream oss;
                    oss << "add_county_dfs_tree_edges saw invalid graph neighbor.\n";
                    oss << "county_id=" << county_id << "\n";
                    oss << "v=" << v << "\n";
                    oss << "u=" << u << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            // because this graph only contains within county edges we don't need to check
            // county status

            // ignore if we've already visited this vertex since it was added to the stack
            if (visited[u]) {
                continue;
            }

            // This is a real graph edge v--u, oriented away from the county root.
            // ust.add_directed_edge(v, u);
            ust[v].push_back(u);
            // mark as visited and add this to stack
            visited[u] = true;
            queue.push(u);

            if constexpr(perf_config::redundancy_checks) ++seen_count;
        }
    }

    // This just checks we visited the number of vertices we expected to visit 
    if constexpr(perf_config::redundancy_checks){
        if (seen_count != n_vtx) {
            std::ostringstream oss;
            oss << "add_county_dfs_tree_edges could not span county-induced subgraph.\n";
            oss << "county_id=" << county_id << "\n";
            oss << "root=" << root << "\n";
            oss << "seen_count=" << seen_count << "\n";
            oss << "n_vtx=" << n_vtx << "\n";

            oss << "Unvisited county vertices: ";
            for (int const v : county_vertices) {
                if (!visited[v]) {
                    oss << v << " ";
                }
            }
            oss << "\n";

            throw std::runtime_error(oss.str());
        }
    }

}

/*
 * Generate a random vertex (integer) among unvisited vertices
 * `lower` is a lower bound (inclusive) on the index of the first unvisited element
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining, int &lower,
         RNGState &rng_state) {
    int idx = rng_state.r_int(remaining);
    int accuml = 0;
    bool seen_one = false;
    for (int i = lower; i < size - 1; i++) {
        accuml += 1 - visited[i];
        if (!seen_one && !visited[i]) {
            seen_one = true;
            lower = i;
        }
        if (accuml - 1 == idx)
            return i;
    }
    return size - 1;
}

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
void loop_erase(std::vector<int> &path, int proposal);

/*
 * Random walk along `county_restricted_graph` from `root` until something in `visited` is hit
 * county_restricted_graph is a graph with all edges going across administrative boundaries removed
 * so we don't need to worry about sampling an edge across counties 
 */
// TESTED
int walk_until(const FlatGraph &county_restricted_graph, 
    int root, std::vector<int> &path, int MAX,
               const std::vector<bool> &visited, const std::vector<bool> &ignore,
               RNGState &rng_state) {
    path[0] = root;
    // walk until we hit something in `visited`
    int curr = root;
    int added = 1; // cursor
    int i;
    for (i = 0; i < MAX; i++) {
        int proposal = rnbor(county_restricted_graph, curr, rng_state);
        if (ignore[proposal]) {
            continue;
        } else if (!visited[proposal]) {
            for (int j = added - 1; j >= 0; j--) {
                if (path[j] == proposal) { // if yes, restart from there
                    added = j;
                    break;
                }
            }
            path[added++] = proposal;
        } else { // reached something in `visited`
            path[added++] = proposal;
            break;
        }
        curr = proposal;
    }
    if (i == MAX) {
        added = 0;
    }

    return added;
}

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::array<int, 3>> &path, int proposal, int root) {
    int length = path.size();
    if (proposal == root) {
        path.erase(path.begin(), path.begin() + length);
        return;
    }

    int idx;
    for (idx = 0; idx < length - 1; idx++) {
        if (path[idx][0] == proposal)
            break;
    }

    if (idx != length - 1) { // a loop
        path.erase(path.begin() + idx + 1, path.begin() + length);
    }
}


/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(const Multigraph &mg, int root, std::vector<std::array<int, 3>> &path,
                    const std::vector<bool> &visited, const std::vector<bool> &ignore,
                    RNGState &rng_state) {
    path.clear();

    // walk until we hit something in `visited`
    int curr = root;
    // while (true) {
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        int prop_idx = rng_state.r_int((int)mg[curr].size());
        int proposal = mg[curr][prop_idx][0];
        if (ignore[mg[curr][prop_idx][2]] || ignore[mg[curr][prop_idx][1]]) {
            continue;
        } else if (!visited[proposal]) {
            path.push_back(mg[curr][prop_idx]);
            loop_erase_cty(path, proposal, root);
        } else {
            path.push_back(mg[curr][prop_idx]);
            break;
        }
        curr = proposal;
    }
    if (i == max) {
        path.clear();
    }
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


// NOTE: in the future add this template option to make it 
// easy to turn off fake tree generation
/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
// TESTED
template <bool SkipUnsplittableTrees>
SampleSubUSTResult sample_sub_ust(MapParams const &map_params, Tree &tree, // FlatGraph &tree, 
                   double const lower,
                   double const upper, std::vector<bool> &visited,
                   const std::vector<bool> &ignore, Tree &county_tree, TreePopStack &county_stack,
                   DummyTreeQueue &dummy_county_tree_queue,
                   std::vector<unsigned int> &county_pop, std::vector<std::vector<int>> &county_members,
                   std::vector<bool> &c_visited, std::vector<int> &cty_pop_below,
                   std::vector<std::array<int, 3>> &county_path, std::vector<int> &path,
                   RNGState &rng_state) {
    // auto t1_start = std::chrono::steady_clock::now();
    int const n_county = map_params.num_counties;
    int tot_pop = 0;
    // reset the county members inner vectors
    // and zero out the county pops
    // and mark all counties as zero
    for (size_t i = 0; i < map_params.num_counties; i++) {
        county_members[i].clear();
        c_visited[i] = true;
        county_pop[i] = 0.0;
    }


    int remaining = 0;
    int const V = map_params.V;
    for (int i = 0; i < V; i++) {
        if (ignore[i]) {
            visited[i] = true;
        } else {
            visited[i] = false;
            remaining++;
            int county = map_params.counties[i] - 1;
            tot_pop += map_params.pop[i];
            county_pop[county] += map_params.pop[i];
            if (c_visited[county]) {
                c_visited[county] = false;
            }
            county_members[county].push_back(i);
        }
    }

    int c_remaining = 0;
    for (int i = 0; i < n_county; i++) {
        c_remaining += 1 - c_visited[i];
    }

    // pick root
    int lower_i = 0;
    int lower_c = 0;
    // int const root = rvtx(visited, V, remaining, lower_i, rng_state);
    int const root = find_unvisited_vertex(V, visited, lower_i);
    visited[root] = true;
    remaining--;
    c_visited[map_params.counties[root] - 1] = true;
    c_remaining--;

    // Connect counties
    // clear the tree
    clear_tree(county_tree);
    county_path.clear();
    while (c_remaining > 0) {
        int add = rvtx(c_visited, n_county, c_remaining, lower_c, rng_state);
        // random walk from `add` until we hit the path
        walk_until_cty(map_params.cg, add, county_path, c_visited, ignore, rng_state);
        // update visited list and constructed tree
        int added = county_path.size();
        if (added == 0) { // bail
            return SampleSubUSTResult{1, root};
        }
        c_remaining -= added;
        c_visited[add] = true;
        for (int i = 0; i < added; i++) {
            c_visited[county_path[i][0]] = true;
            // reverse path so that arrows point away from root
            // tree.add_directed_edge(county_path[i][2], county_path[i][1]);
            tree[county_path[i][2]].push_back(county_path[i][1]);
            county_tree[county_path[i][0]]
                .push_back(map_params.counties[county_path[i][1]] - 1);

            visited[county_path[i][1]] = true; // root for next district
            remaining--;
        }
    }

    // optional toggle of skipping unsplittable trees 
    if constexpr (SkipUnsplittableTrees) {
    // figure out which counties will not need to be split
    if (n_county > 1) {
        // don't need to fill pop below since it gets reset
        OLD_TO_UPDATE_get_tree_pops_below(county_tree, map_params.counties[root] - 1, county_stack, county_pop,
                            cty_pop_below);
        for (int i = 0; i < n_county; i++) {
            int n_vtx = county_members[i].size();
            if (n_vtx <= 1)
                continue;
            // check child counties
            int children = county_tree[i].size();
            int split_ub = cty_pop_below[i];
            int split_lb = split_ub - county_pop[i];
            if (lower - 1 < county_pop[i])
                split_lb = (int)lower;
            for (int j = 0; j < children; j++) {
                int pop_child = cty_pop_below[county_tree[i][j]];
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
            bool miss_second = (tot_pop - split_lb) < lower || (tot_pop - split_ub) > upper;

            // impossible for this county to need to be split
            // so we fill in a dummy tree 
            if (cty_pop_below[i] >= 0 && (miss_first && miss_second)) {
                // check if we need the dummy tree to actually be a subset of g
                add_county_to_tree_dfs(
                    map_params.county_restricted_graph,
                    map_params.counties[county_members[i][0]], // pass in the county CAREFUL COUNTIES ARE 1-indexed
                    county_members[i],
                    visited,
                    dummy_county_tree_queue,
                    tree
                );
                remaining -= n_vtx - 1; // already visited county root
            }
        }
    }
    }

    // Generate tree within each county
    if (remaining > 0) {
        path.clear();
        path.resize(remaining + 2);
        int max_try = 50 * remaining * (static_cast<int>(std::log(remaining)) + 1);
        while (remaining > 0) {
            // int add = rvtx(visited, V, remaining, lower_i, rng_state);
            int add = find_unvisited_vertex(V, visited, lower_i);
            // random walk from `add` until we hit the path
            int added = walk_until(map_params.county_restricted_flat_graph, add, path, max_try, visited, ignore,
                                   rng_state);
            // update visited list and constructed tree
            if (added == 0) { // bail
                return SampleSubUSTResult{1, root};
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited[path[i]] = true;
                // reverse path so that arrows point away from root
                // tree.add_directed_edge(path[i + 1], path[i]);
                tree[path[i + 1]].push_back(path[i]);
            }
        }
    }

    return SampleSubUSTResult{0, root};
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
        int const parent =
            next_vertex[curr_vertex];

        if constexpr (perf_config::supposedly_safe_input_checks) {
            if (visited[parent]) {
                throw std::runtime_error(
                    "add_walk_from: active path contains an already visited vertex."
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


// We assume that ignore has been properly set and nothing else
// has been cleared 
template <bool SkipUnsplittableTrees>
SampleSubUSTResult sample_ust(
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
    if constexpr(DEBUG_WILSON_VERBOSE){
        std::cerr << "Starting with root of " << root << std::endl;
    }
    visited[root] = true;
    g_scratch.remaining--;

    int lower_c = 0;
    mg_scratch.c_visited[map_params.counties[root] - 1] = true;
    mg_scratch.c_remaining--;

    // Connect counties
    // clear the tree
    clear_tree(county_tree);
    mg_scratch.county_path.clear();
    while (mg_scratch.c_remaining > 0) {
        int add = rvtx(mg_scratch.c_visited, n_county, mg_scratch.c_remaining, lower_c, rng_state);
        // random walk from `add` until we hit the path
        walk_until_cty(map_params.cg, add, mg_scratch.county_path, mg_scratch.c_visited, ignore, rng_state);
        // update visited list and constructed tree
        int added = mg_scratch.county_path.size();
        if (added == 0) { // bail
            return SampleSubUSTResult{1, root};
        }
        mg_scratch.c_remaining -= added;
        mg_scratch.c_visited[add] = true;
        for (int i = 0; i < added; i++) {
            mg_scratch.c_visited[mg_scratch.county_path[i][0]] = true;
            // reverse path so that arrows point away from root
            // tree.add_directed_edge(county_path[i][2], county_path[i][1]);
            tree[mg_scratch.county_path[i][2]].push_back(mg_scratch.county_path[i][1]);
            county_tree[mg_scratch.county_path[i][0]]
                .push_back(map_params.counties[mg_scratch.county_path[i][1]] - 1);

            visited[mg_scratch.county_path[i][1]] = true; // root for next district
            g_scratch.remaining--;
        }
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
        for (int i = 0; i < n_county; i++) {
            int n_vtx = mg_scratch.county_members[i].size();
            if (n_vtx <= 1)
                continue;
            // check child counties
            int children = county_tree[i].size();
            int split_ub = mg_scratch.cty_pop_below[i];
            int split_lb = split_ub - mg_scratch.county_pop[i];
            if (lower - 1 < mg_scratch.county_pop[i])
                split_lb = (int)lower;
            for (int j = 0; j < children; j++) {
                int pop_child = mg_scratch.cty_pop_below[county_tree[i][j]];
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
            if (mg_scratch.cty_pop_below[i] >= 0 && (miss_first && miss_second)) {
                // check if we need the dummy tree to actually be a subset of g
                add_county_to_tree_dfs(
                    map_params.county_restricted_graph,
                    map_params.counties[mg_scratch.county_members[i][0]], // pass in the county CAREFUL COUNTIES ARE 1-indexed
                    mg_scratch.county_members[i],
                    visited,
                    g_scratch.dummy_county_tree_queue,
                    tree
                );
                g_scratch.remaining -= n_vtx - 1; // already visited county root
            }
        }
    }
    }

    if constexpr(DEBUG_WILSON_VERBOSE){
        std::cerr << "Starting to draw within county trees!\n";
    }
    
    // Generate tree within each county
    if (g_scratch.remaining > 0) {
        g_scratch.path.clear();
        g_scratch.path.resize(g_scratch.remaining + 2);
        int max_try = 50 * g_scratch.remaining * (static_cast<int>(std::log(g_scratch.remaining)) + 1);
        while (g_scratch.remaining > 0) {
            int unvisited_vertex = find_unvisited_vertex(V, visited, g_scratch.smallest_v_seen);
            // int unvisited_vertex = rvtx(visited, V, g_scratch.remaining, g_scratch.smallest_v_seen, rng_state);
            // random walk from `unvisited_vertex` until we hit the path
            int vertices_added = add_walk_from(
                map_params.county_restricted_flat_graph,
                tree,
                unvisited_vertex,
                max_try,
                g_scratch.path,
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

int USTSampler::prep_fresh_ust_call(){
    // We assume that ignore has already been properly set 
    // AND NOTHING ELSE. No trees have been cleared or graphs touched
    // We will now 
    // - properly set the visited vector 
    // - Set the wilson_submap so its the subgraph restricted to vertices and 
    //   edges solely contained in the `ignore[v] == false`
    // - TODO: set up county multigraph stuff 

    mg_scratch.c_remaining = 0;
    mg_scratch.total_pop = 0;
    // reset the county members inner vectors
    // and zero out the county pops
    // and mark all counties as zero
    for (size_t i = 0; i < map_params.num_counties; i++) {
        mg_scratch.county_members[i].clear();
        mg_scratch.c_visited[i] = true;
        mg_scratch.county_pop[i] = 0.0;
    }

    bool first_county_seen = false;

    g_scratch.remaining = 0;
    int const V = map_params.V;
    for (int v = 0; v < V; v++)
    {
        if (ignore[v]) {
            // if we're ignoring we just mark as visited and continue
            visited[v] = true;
            // there's nothing more to do 
            continue;
        } 
        // Else it means v is part of the subgraph we care about 
        visited[v] = false;
        g_scratch.remaining++;
        auto v_county = map_params.counties[v] - 1;
        mg_scratch.total_pop += map_params.pop[v];
        mg_scratch.county_pop[v_county] += map_params.pop[v];
        if (mg_scratch.c_visited[v_county]) {
            mg_scratch.c_visited[v_county] = false;
            ++mg_scratch.c_remaining;
            // This will be our initial start vertex since this will 
            // TODO add another check 
            if(!first_county_seen){
                g_scratch.smallest_v_seen = v;
                first_county_seen = true;
            }

            // TODO add this vertex to county root vector 
            // TODO in the future clear the county multigraph and add this edge 
        }
        mg_scratch.county_members[v_county].push_back(v);
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

USTDrawResult USTSampler::draw_fresh_ust(
      double const lower, double const upper,
        RNGState &rng_state) {


    // prep the inputs and get the number of vertices 
    int const num_tree_vertices = prep_fresh_ust_call();
    // REprintf("%d Vertices!\n", num_tree_vertices);

    auto const result = sample_ust<true>(map_params, ust, 
        lower, upper, visited, ignore,
        county_tree, g_scratch, mg_scratch, rng_state
    );

    if constexpr(perf_config::object_integrity_checking){
        if (result.code == 0){
            check_tree_integrity(
                get_vertex_tree(),
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_region\n",
                result.root,
                num_tree_vertices,
                true
            );
        }
    }

    // result.code == 0 means it was successful
    return USTDrawResult{result.code == 0, num_tree_vertices, result.root};
}
 

std::pair<bool, int> USTSampler::draw_ust(
      double const lower, double const upper,
        RNGState &rng_state) {
    // We assume that ignore has already been properly set 
    // We also assume the tree has been properly cleared
    // Now get a uniform spanning tree drawn on the subgraph denoted by the ignore
    // vertices 
    auto const result = sample_sub_ust<true>(map_params, ust, lower, upper, visited, ignore,
                                county_tree, mg_scratch.county_stack, g_scratch.dummy_county_tree_queue,
                                mg_scratch.county_pop, mg_scratch.county_members,
                                mg_scratch.c_visited, mg_scratch.cty_pop_below, 
                                mg_scratch.county_path, g_scratch.path, rng_state
                                );
    // result.code == 0 means it was successful
    return std::make_pair(result.code == 0, result.root);
}


USTDrawResult USTSampler::OLD_attempt_to_draw_tree_on_region(
    RNGState &rng_state, Plan const &plan, const int region_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;

    int num_region_vertices = 0;
    // Mark it as ignore if its not in the region to split and clear those vertices in 
    // the region
    for (int i = 0; i < V; ++i) {
        // check if in the region
        bool const in_region = plan.region_ids[i] == region_to_draw_tree_on;
        ignore[i] = !in_region;
        // clear if in the region
        if (in_region) {
            // ust.clear_vertex(i);
            ust[i].clear();
            ++num_region_vertices;
        }
    }

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
    auto const result = draw_ust(
        sample_sub_ust_lower, // lower, 
        sample_sub_ust_upper, // upper
        rng_state
        );
    bool const valid_tree = result.first;

    if constexpr(perf_config::object_integrity_checking){
        if (valid_tree){
            check_tree_integrity(
                get_vertex_tree(),
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_region\n",
                result.second,
                num_region_vertices,
                true
            );
        }
    }

    return USTDrawResult{result.first, num_region_vertices, result.second};
}


USTDrawResult USTSampler::attempt_to_draw_tree_on_region(
    RNGState &rng_state, Plan const &plan, const int region_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;

    // Mark it as ignore if its not in the region to split and clear those vertices in 
    // the region
    for (int i = 0; i < V; ++i) {
        // check if in the region
        // TODO: in the future make this all a single pass through
        ignore[i] = plan.region_ids[i] != region_to_draw_tree_on;
    }

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
    auto const result = draw_fresh_ust(
        sample_sub_ust_lower, // lower, 
        sample_sub_ust_upper, // upper
        rng_state
        );

    return result;
}


USTDrawResult USTSampler::OLD_attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
                                                       const int region1_to_draw_tree_on,
                                                       const int region2_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;

    // optional for checking 
    int num_merged_region_vertices = 0;

    // mark the ignore values and clear vertices in either region
    for (int i = 0; i < V; ++i) {
        // check if in the region
        bool const in_region = plan.region_ids[i] == region1_to_draw_tree_on ||
                    plan.region_ids[i] == region2_to_draw_tree_on;

        ignore[i] = !in_region;
        // clear if in the region
        if (in_region) {
            // ust.clear_vertex(i);
            ust[i].clear();
            ++num_merged_region_vertices;
        }
    }

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
    auto const result = draw_ust(
        sample_sub_ust_lower, // lower, 
        sample_sub_ust_upper, // upper
        rng_state
        );
    bool const valid_tree = result.first;

    if constexpr(perf_config::object_integrity_checking){
        if (valid_tree){
            check_tree_integrity(
                get_vertex_tree(),
                "Just called `sample_sub_ust` in attempt_to_draw_tree_on_merged_region\n",
                result.second,
                num_merged_region_vertices,
                true
            );
        }
    }

    return USTDrawResult{valid_tree, num_merged_region_vertices, result.second};
}


USTDrawResult USTSampler::attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
                                                       const int region1_to_draw_tree_on,
                                                       const int region2_to_draw_tree_on,
    bool const use_custom_bounds, double const custom_sample_sub_ust_lower,
    double const custom_sample_sub_ust_upper) {
    int V = map_params.V;


    // mark the ignore values
    for (int i = 0; i < V; ++i) {
        // check if in the region
        ignore[i] = plan.region_ids[i] != region1_to_draw_tree_on &&
                    plan.region_ids[i] != region2_to_draw_tree_on;
    }

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
    auto const result = draw_fresh_ust(
        sample_sub_ust_lower, // lower, 
        sample_sub_ust_upper, // upper
        rng_state
        );

    return result;
}

void USTSampler::fill_in_skipped_subtrees(RNGState &rng_state, int max_tries){
    // TODO
    // Need to mark all the vertices to fill in as not visited, everything else is fine
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
        pops_below_vertex, ignore, region_populations, region_size,
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
    int num_vertices = 0;

    // Mark it as ignore if its not in the region to split and clear those vertices in 
    // the region
    for (int i = 0; i < V; ++i) {
        ignore[i] = vertices_to_ignore[i];
        // check if we're ignoring it
        if (!ignore[i]) {
            // ust.clear_vertex(i);
            ust[i].clear();
            num_vertices++;
        }
    }

    // prepare the inputs for a `sample_ust` call
    prep_fresh_ust_call();

    auto prep_end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::ratio<1>> prep_time_diff = prep_end_time - prep_start_time;
    wilson_times.input_prep_time += prep_time_diff.count();

    std::pair<bool, int> result;

    auto sample_ust_start_time = std::chrono::steady_clock::now();
    // Now sample a uniform spanning tree drawn on that region
    if (skip_unsplittable_subtrees){
        auto const ust_result = sample_ust<true>(map_params, ust, 
            lower, upper, visited, ignore,
            county_tree, g_scratch, mg_scratch, rng_state
        );
        result = std::make_pair(ust_result.code == 0, ust_result.root);
    }else{
        auto const ust_result = sample_ust<false>(map_params, ust, 
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


std::pair<bool, int>  USTSampler::OLD_draw_tree_on_subgraph(
      RNGState &rng_state, std::vector<bool> const &vertices_to_ignore,
      bool const skip_unsplittable_subtrees, 
      const double lower, const double upper,
      WilsonTimes &wilson_times
){
    auto prep_start_time = std::chrono::steady_clock::now();

    int V = map_params.V;
    int num_vertices = 0;

    // Mark it as ignore if its not in the region to split and clear those vertices in 
    // the region
    for (int i = 0; i < V; ++i) {
        ignore[i] = vertices_to_ignore[i];
        // check if we're ignoring it
        if (!ignore[i]) {
            // ust.clear_vertex(i);
            ust[i].clear();
            num_vertices++;
        }
    }

    auto prep_end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::ratio<1>> prep_time_diff = prep_end_time - prep_start_time;
    wilson_times.input_prep_time += prep_time_diff.count();

    std::pair<bool, int> result;

    auto sample_ust_start_time = std::chrono::steady_clock::now();
    // result = draw_fresh_ust(lower, upper, rng_state);
    // auto sample_ust_end_time = std::chrono::steady_clock::now();
    // std::chrono::duration<double, std::ratio<1>> ust_time_diff = sample_ust_end_time - sample_ust_start_time;
    // wilson_times.sub_ust_call_time += ust_time_diff.count();
    // return result;

    // Now sample a uniform spanning tree drawn on that region
    if (skip_unsplittable_subtrees){
        auto const ust_result = sample_sub_ust<true>(map_params, ust, lower, upper, visited, ignore,
                                county_tree, mg_scratch.county_stack, g_scratch.dummy_county_tree_queue,
                                mg_scratch.county_pop, mg_scratch.county_members,
                                mg_scratch.c_visited, mg_scratch.cty_pop_below, 
                                mg_scratch.county_path, g_scratch.path, rng_state
                                    );
        result = std::make_pair(ust_result.code == 0, ust_result.root);
    }else{
        auto const ust_result = sample_sub_ust<false>(map_params, ust, lower, upper, visited, ignore,
                                county_tree, mg_scratch.county_stack, g_scratch.dummy_county_tree_queue,
                                mg_scratch.county_pop, mg_scratch.county_members,
                                mg_scratch.c_visited, mg_scratch.cty_pop_below, 
                                mg_scratch.county_path, g_scratch.path, rng_state
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