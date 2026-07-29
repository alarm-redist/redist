#include <RcppArmadillo.h>
#include "wilson.h"

#include "splitting_schedule_types.h"
#include "random.h"
#include "base_plan_type.h"
#include "tree_splitting.h"

namespace{


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
    int const root = rvtx(visited, V, remaining, lower_i, rng_state);
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
            int add = rvtx(visited, V, remaining, lower_i, rng_state);
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
 * Advance the epoch stamp used for identifying vertices touched during the
 * current loop-erased random walk.
 *
 * The idea:
 *
 *   walk_epochs[v] == current_walk_epoch
 *
 * means that vertex v has been touched during the current walk. This lets us
 * avoid clearing a length-V array at the beginning of every walk.
 *
 * However, uint32_t can eventually overflow. If we simply wrapped from
 * UINT32_MAX back to 1 without clearing `walk_epochs`, old vertices with stamp
 * 1 could be incorrectly interpreted as belonging to the new walk.
 *
 * So on overflow we clear all stamps back to zero and restart at epoch 1.
 * In practice, overflow is extremely unlikely, but this makes the logic
 * formally safe.
 */
inline void advance_walk_epoch(
    std::vector<std::uint32_t> &walk_epochs,
    std::uint32_t &current_walk_epoch
) {
    if (current_walk_epoch == std::numeric_limits<std::uint32_t>::max()) {
        std::fill(walk_epochs.begin(), walk_epochs.end(), 0);
        current_walk_epoch = 1;
    } else {
        ++current_walk_epoch;
    }
}

/*
 * Return true if vertex `v` is currently on the active loop-erased walk path.
 *
 * This function is the heart of the no-explicit-erase logic.
 *
 * We maintain a reusable path buffer:
 *
 *   path[0], path[1], ..., path[active_len - 1]
 *
 * These entries are the current loop-erased path.
 *
 * Entries path[active_len], path[active_len + 1], ... may contain stale
 * vertices from earlier loops. We do not physically erase them.
 *
 * For each vertex v, path_pos[v] stores the last slot where v was placed.
 * But because we do not clear stale positions, path_pos[v] alone is not enough.
 *
 * A vertex v is truly on the current path iff all three conditions hold:
 *
 *   1. walk_epochs[v] == current_walk_epoch
 *      v has been touched during this walk.
 *
 *   2. path_pos[v] < active_len
 *      v's recorded slot is inside the currently active prefix.
 *
 *   3. path[path_pos[v]] == v
 *      the slot still actually contains v. This rules out stale vertices whose
 *      old slot has since been overwritten by another vertex.
 *
 * The third condition is important.
 *
 * Example:
 *
 *   current path: 1 -> 2 -> 3 -> 4
 *   loop to 2
 *
 * Then active_len becomes 2, so the active path is:
 *
 *   1 -> 2
 *
 * The entries 3 and 4 are stale. Later the walk may add vertex 5 at slot 2,
 * overwriting path[2]. Without checking path[path_pos[3]] == 3, vertex 3
 * might falsely look active once active_len grows beyond 2 again.
 */
inline bool is_on_current_walk_path(
    int const v,
    std::vector<std::uint32_t> const &walk_epochs,
    std::uint32_t const current_walk_epoch,
    std::vector<int> const &path,
    std::vector<int> const &path_pos,
    int const active_len
) {
    if (walk_epochs[v] != current_walk_epoch) {
        return false;
    }

    int const pos = path_pos[v];

    return pos >= 0 &&
           pos < active_len &&
           path[pos] == v;
}

/*
 * Perform one Wilson loop-erased random walk and add the resulting path to
 * an undirected tree.
 *
 * High-level purpose:
 *
 *   Starting from `start_vertex`, which is not yet in the tree, perform a
 *   random walk until the walk hits a vertex that is already in the tree.
 *   Erase loops chronologically as they appear. Once the walk hits the
 *   existing tree, add the loop-erased path to `tree`.
 *
 * Important assumptions:
 *
 *   1. `tree` is an undirected adjacency-list tree with size graph.size().
 *
 *   2. `visited[v] != 0` means vertex v is already in the current tree.
 *
 *   3. `start_vertex` must not already be visited.
 *
 *   4. `graph` must be the graph on which the random walk is supposed to run.
 *
 *      If `graph` is the induced subgraph for the current district/region,
 *      then no extra active/ignore check is needed.
 *
 *      If `graph` is only county-restricted, but still includes vertices from
 *      other regions, then this function is not enough as written. In that
 *      case you must add an active-region check after sampling `curr_nbor`.
 *
 * Scratch vectors:
 *
 *   walk_epochs:
 *      length V. Stamp array used to identify vertices touched in the current
 *      walk without clearing all V entries each time.
 *
 *   current_walk_epoch:
 *      shared epoch counter. This function advances it once at the beginning
 *      of each walk.
 *
 *   path:
 *      reusable buffer storing the current loop-erased path in forward order:
 *
 *          path[0]                  = start_vertex
 *          path[active_len - 1]     = current endpoint of walk
 *
 *      Only the prefix path[0:active_len) is active. Entries beyond active_len
 *      are stale and ignored.
 *
 *   path_pos:
 *      length V. path_pos[v] is the last slot where v was written into `path`.
 *      This is only meaningful when combined with walk_epochs and the slot
 *      verification path[path_pos[v]] == v.
 *
 *   visited:
 *      length V. Byte-valued flag array. The function marks all newly-added
 *      path vertices as visited after successfully hitting the tree.
 *
 * Return value:
 *
 *   true:
 *      The walk hit the existing tree and the loop-erased path was added.
 *
 *   false:
 *      The walk did not hit the existing tree within max_tries proposals.
 *      In this case, `tree` and `visited` are not modified.
 */
bool add_walk_from(
    FlatGraph const &graph,
    Tree &tree,
    int const start_vertex,
    std::vector<std::uint32_t> &walk_epochs,
    std::uint32_t &current_walk_epoch,
    int const max_tries,
    std::vector<int> &path,
    std::vector<int> &path_pos,
    std::vector<bool> &visited,
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

        if (static_cast<int>(walk_epochs.size()) != V) {
            throw std::runtime_error("add_walk_from: walk_epochs has wrong size.");
        }

        if (static_cast<int>(path_pos.size()) != V) {
            throw std::runtime_error("add_walk_from: path_pos has wrong size.");
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

    /*
     * Start a new walk epoch.
     *
     * This lets us distinguish vertices touched during this walk from vertices
     * touched during previous calls to add_walk_from, without clearing the
     * entire `walk_epochs` vector.
     */
    advance_walk_epoch(walk_epochs, current_walk_epoch);

    /*
     * Initialize the current loop-erased path.
     *
     * We reuse `path` as a buffer. If it already has storage, overwrite slot 0.
     * Otherwise push the start vertex.
     *
     * active_len is the logical length of the current loop-erased path.
     * Only path[0], ..., path[active_len - 1] are part of the current path.
     */
    int active_len = 1;

    if (path.empty()) {
        path.push_back(start_vertex);
    } else {
        path[0] = start_vertex;
    }

    /*
     * Record that start_vertex is on the current walk path at position 0.
     */
    path_pos[start_vertex] = 0;
    walk_epochs[start_vertex] = current_walk_epoch;

    /*
     * curr_vertex is always the last vertex in the active path:
     *
     *   curr_vertex == path[active_len - 1]
     *
     * We maintain this invariant throughout the random walk.
     */
    int curr_vertex = start_vertex;

    /*
     * hit_vertex will become the first already-visited tree vertex hit by the
     * random walk.
     *
     * We do not include hit_vertex inside `path`, because it is already in the
     * tree. At the end, the last path vertex is connected to hit_vertex.
     */
    int hit_vertex = -1;

    /*
     * Main random-walk loop.
     *
     * Each iteration proposes one random neighbor of curr_vertex.
     *
     * There are three possible cases:
     *
     *   Case 1: proposed neighbor is already in the tree.
     *           The walk is complete.
     *
     *   Case 2: proposed neighbor is already on the current walk path.
     *           A loop has formed. Shorten active_len to erase the loop.
     *
     *   Case 3: proposed neighbor is new to the current active path.
     *           Append it to the path.
     */
    for (int tries = 0; tries < max_tries; ++tries) {
        /*
         * Draw a random neighbor from the graph.
         *
         * This assumes rnbor(graph, curr_vertex, rng_state) returns a uniformly
         * random neighbor of curr_vertex in the graph.
         *
         * If graph is not already induced to the active region, insert an
         * active/ignore rejection check immediately after this line.
         */
        int const curr_nbor = rnbor(graph, curr_vertex, rng_state);

        /*
         * Case 1: the random walk has hit the existing tree.
         *
         * We stop immediately. The current active path connects to hit_vertex.
         */
        if (visited[curr_nbor]) {
            hit_vertex = curr_nbor;
            break;
        }

        /*
         * Case 2: curr_nbor is already on the current loop-erased path.
         *
         * This means the random walk has formed a loop:
         *
         *   curr_nbor -> ... -> curr_vertex -> curr_nbor
         *
         * Chronological loop erasure removes everything after the first
         * occurrence of curr_nbor. Because curr_nbor is already in `path`, and
         * path_pos[curr_nbor] gives its location, this is just:
         *
         *   active_len = path_pos[curr_nbor] + 1
         *
         * We do not clear stale path entries. We simply shrink the active prefix.
         *
         * Example:
         *
         *   path:       1 2 3 4
         *   active_len: 4
         *   proposal:   2
         *
         * New active path:
         *
         *   path:       1 2 [3 4 are stale]
         *   active_len: 2
         */
        if (is_on_current_walk_path(
                curr_nbor,
                walk_epochs,
                current_walk_epoch,
                path,
                path_pos,
                active_len
            )) {
            active_len = path_pos[curr_nbor] + 1;
            curr_vertex = curr_nbor;
            continue;
        }

        /*
         * Case 3: curr_nbor is not currently on the active path.
         *
         * It may be either:
         *
         *   - a vertex never touched in this walk, or
         *   - a stale vertex from an erased loop, or
         *   - a vertex from a previous walk epoch.
         *
         * In all cases, we append it to the current active path.
         *
         * Because `path` is alawys guaranteed to be of length V
         * we just overrite path[active_len]
         *
        path[active_len] = curr_nbor;
        

        /*
         * Record curr_nbor's current position and stamp.
         *
         * The slot verification path[path_pos[v]] == v is what makes stale
         * entries safe, even if old path_pos values remain.
         */
        path_pos[curr_nbor] = active_len;
        walk_epochs[curr_nbor] = current_walk_epoch;

        /*
         * The active path is now one vertex longer, and curr_nbor is the new
         * endpoint.
         */
        ++active_len;
        curr_vertex = curr_nbor;
    }

    /*
     * If we did not hit the existing tree within max_tries, report failure.
     *
     * Notice that up to this point, we have only modified scratch arrays:
     *
     *   path
     *   path_pos
     *   walk_epochs
     *
     * We have not modified:
     *
     *   tree
     *   visited
     *
     * So returning false does not leave a partially-added tree path.
     */
    if (hit_vertex == -1) {
        return false;
    }

    /*
     * Add the loop-erased path to the tree.
     *
     * The active path is stored in forward random-walk order:
     *
     *   path[0]              = start_vertex
     *   path[active_len - 1] = last unvisited vertex before hitting the tree
     *
     * and we know:
     *
     *   path[active_len - 1] is adjacent to hit_vertex
     *
     * Since hit_vertex is already in the tree, the easiest way to attach the
     * path is to walk backward:
     *
     *   hit_vertex
     *      ^
     *      |
     *   path[active_len - 1]
     *      ^
     *      |
     *   path[active_len - 2]
     *      ^
     *      |
     *   ...
     *      ^
     *      |
     *   path[0]
     *
     * For each child-parent pair, add both directions because `tree` is
     * undirected.
     */
    int parent = hit_vertex;

    for (int i = active_len - 1; i >= 0; --i) {
        int const child = path[i];

        /*
         * Defensive checks.
         *
         * The child should not already be in the tree. If it were, then the walk
         * should have stopped when it first proposed that vertex.
         */
        if constexpr (perf_config::supposedly_safe_input_checks) {
            if (visited[child]) {
                throw std::runtime_error(
                    "add_walk_from: active path contains an already visited vertex."
                );
            }
        }

        /*
         * Add the undirected tree edge child -- parent.
         */
        tree[child].push_back(parent);
        tree[parent].push_back(child);

        /*
         * Mark child as now belonging to the tree.
         *
         * This is essential. Future Wilson walks must treat this vertex as a
         * valid stopping point.
         */
        visited[child] = true;

        /*
         * Move backward along the path. On the next iteration, this child will
         * be the parent of the previous path vertex.
         */
        parent = child;
    }

    return true;
}


// We assume that ignore has been properly set and nothing else
// has been cleared 
// template <bool SkipUnsplittableTrees>
// SampleSubUSTResult sample_ust(
//     MapParams const &map_params, 
//     FlatGraph &wilson_submap, 
//     Tree &tree, // FlatGraph &tree, 
//     double const lower, double const upper, 
//     std::vector<bool> &visited, const std::vector<bool> &ignore, 
//     Tree &county_tree, 
//     WilsonGraphScratch &g_scratch,
//     WilsonMultiGraphScratch &mg_scratch,
//     RNGState &rng_state
//     ) {
//     // auto t1_start = std::chrono::steady_clock::now();
//     int const n_county = map_params.num_counties;
//     int tot_pop = 0;
//     // reset the county members inner vectors
//     // and zero out the county pops
//     // and mark all counties as zero
//     for (size_t i = 0; i < map_params.num_counties; i++) {
//         county_members[i].clear();
//         c_visited[i] = true;
//         county_pop[i] = 0.0;
//     }

//     int c_remaining = 0;
//     int remaining = 0;
//     int const V = map_params.V;
//     for (int v = 0; v < V; v++)
//     {
//         if (ignore[v]) {
//             // if we're ignoring we just continue as 
//             // there's nothing more to do 
//             continue;
//         } 
//         // Else it means v is part of the subgraph we care about 
//         visited[v] = false;
//         remaining++;
//         auto county = map_params.counties[v] - 1;
//         tot_pop += map_params.pop[v];
//         county_pop[county] += map_params.pop[v];
//         if (c_visited[county]) {
//             c_visited[county] = false;
//             c_remaining++;
//         }
//         county_members[county].push_back(v);
//         // Now we clear the ust at this vertex
//         // ust.clear_vertex(i);
//         tree[v].clear();
//         // Now we clear the subgraph 
//         wilson_submap.clear_vertex(v);
//         // now we add all the neighbors not ignored 
//         for (auto const v_nbor: map_params.g[v])
//         {
//             if (!ignore[v_nbor]){
//                 wilson_submap.add_directed_edge(v, v_nbor);
//             }
//         }
        
//     }



//     // pick root
//     int lower_i = 0;
//     int lower_c = 0;
//     int const root = rvtx(visited, V, remaining, lower_i, rng_state);
//     visited[root] = true;
//     remaining--;
//     c_visited[map_params.counties[root] - 1] = true;
//     c_remaining--;

//     // Connect counties
//     // clear the tree
//     clear_tree(county_tree);
//     county_path.clear();
//     while (c_remaining > 0) {
//         int add = rvtx(c_visited, n_county, c_remaining, lower_c, rng_state);
//         // random walk from `add` until we hit the path
//         walk_until_cty(map_params.cg, add, county_path, c_visited, ignore, rng_state);
//         // update visited list and constructed tree
//         int added = county_path.size();
//         if (added == 0) { // bail
//             return SampleSubUSTResult{1, root};
//         }
//         c_remaining -= added;
//         c_visited[add] = true;
//         for (int i = 0; i < added; i++) {
//             c_visited[county_path[i][0]] = true;
//             // reverse path so that arrows point away from root
//             // tree.add_directed_edge(county_path[i][2], county_path[i][1]);
//             tree[county_path[i][2]].push_back(county_path[i][1]);
//             county_tree[county_path[i][0]]
//                 .push_back(map_params.counties[county_path[i][1]] - 1);

//             visited[county_path[i][1]] = true; // root for next district
//             remaining--;
//         }
//     }

//     // optional toggle of skipping unsplittable trees 
//     // if constexpr (SkipUnsplittableTrees) {
//     // figure out which counties will not need to be split
//     if (n_county > 1) {
//         // don't need to fill pop below since it gets reset
//         OLD_TO_UPDATE_get_tree_pops_below(county_tree, map_params.counties[root] - 1, county_stack, county_pop,
//                             cty_pop_below);
//         for (int i = 0; i < n_county; i++) {
//             int n_vtx = county_members[i].size();
//             if (n_vtx <= 1)
//                 continue;
//             // check child counties
//             int children = county_tree[i].size();
//             int split_ub = cty_pop_below[i];
//             int split_lb = split_ub - county_pop[i];
//             if (lower - 1 < county_pop[i])
//                 split_lb = (int)lower;
//             for (int j = 0; j < children; j++) {
//                 int pop_child = cty_pop_below[county_tree[i][j]];
//                 if (pop_child >= 0 && pop_child < split_lb) {
//                     split_lb = pop_child;
//                 }
//             }

//             // split_lb < split_ub so the smallest possible population is
//             // min(split_lb, tot_pop - split_ub)
//             // its impossible to split if smallest possible size is bigger than largest ub
//             // bool miss_first = std::min(split_lb, tot_pop - split_ub) > upper;
//             // // biggest possible population is max(split_ub, total_pop - split_lb)
//             // // its impossible to split if largest possible size is smaller than smallest lb
//             // bool miss_second = std::max(split_ub, tot_pop - split_lb) < lower;

//             bool miss_first = split_ub < lower || split_lb > upper;
//             bool miss_second = (tot_pop - split_lb) < lower || (tot_pop - split_ub) > upper;

//             // impossible for this county to need to be split
//             // so we fill in a dummy tree 
//             if (cty_pop_below[i] >= 0 && (miss_first && miss_second)) {
//                 // check if we need the dummy tree to actually be a subset of g
//                 add_county_to_tree_dfs(
//                     map_params.county_restricted_graph,
//                     map_params.counties[county_members[i][0]], // pass in the county CAREFUL COUNTIES ARE 1-indexed
//                     county_members[i],
//                     visited,
//                     dummy_county_tree_queue,
//                     tree
//                 );
//                 remaining -= n_vtx - 1; // already visited county root
//             }
//         }
//     }
//     // }

//     // Generate tree within each county
//     if (remaining > 0) {
//         path.clear();
//         path.resize(remaining + 2);
//         int max_try = 50 * remaining * (static_cast<int>(std::log(remaining)) + 1);
//         while (remaining > 0) {
//             int add = rvtx(visited, V, remaining, lower_i, rng_state);
//             // random walk from `add` until we hit the path
//             int added = walk_until(map_params.county_restricted_flat_graph, add, path, max_try, visited, ignore,
//                                    rng_state);
//             // update visited list and constructed tree
//             if (added == 0) { // bail
//                 return SampleSubUSTResult{1, root};
//             }
//             remaining -= added - 1; // minus 1 because ending vertex already in tree
//             for (int i = 0; i < added - 1; i++) {
//                 visited[path[i]] = true;
//                 // reverse path so that arrows point away from root
//                 // tree.add_directed_edge(path[i + 1], path[i]);
//                 tree[path[i + 1]].push_back(path[i]);
//             }
//         }
//     }

//     return SampleSubUSTResult{0, root};
// }

}





/********************************************************
 * Author: Philip O'Sullivan'
 * Institution: Harvard University
 * Date Created: 2025/3
 * Purpose: Encapsulation of uniform spanning tree sampler functions
 ********************************************************/


 // The idea is to move all this outside of sample sub ust






// We assume that ignore has been correctly set, nothing else
// INCLUDING THE TREE HAS NOT BEEN CLEARED!
// std::pair<bool, int> USTSampler::NEW_draw_ust(
//       double const lower, double const upper,
//         RNGState &rng_state) {
//     // We assume that ignore has already been properly set 
//     // AND NOTHING ELSE. No trees have been cleared or graphs touched

//     // We will now 
//     // - properly set the visited vector 
//     // - Set the wilson_submap so its the subgraph restricted to vertices and 
//     //   edges solely contained in the `ignore[v] == false`
//     // - TODO: set up county multigraph stuff 

//     int const n_county = map_params.num_counties;
//     int tot_pop = 0;
//     // reset the county members inner vectors
//     // and zero out the county pops
//     // and mark all counties as zero
//     for (size_t i = 0; i < map_params.num_counties; i++) {
//         county_members[i].clear();
//         c_visited[i] = true;
//         county_pop[i] = 0.0;
//     }

//     int remaining = 0;
//     int const V = map_params.V;
//     for (int v = 0; v < V; v++)
//     {
//         if (ignore[v]) {
//             // if we're ignoring we just continue as 
//             // there's nothing more to do 
//             continue;
//         } 
//         // Else it means v is part of the subgraph we care about 
//         visited[v] = false;
//         remaining++;
//         auto county = map_params.counties[v] - 1;
//         tot_pop += map_params.pop[v];
//         county_pop[county] += map_params.pop[v];
//         if (c_visited[county]) {
//             c_visited[county] = false;
//         }
//         county_members[county].push_back(v);
//         // Now we clear the ust at this vertex
//         // ust.clear_vertex(i);
//         ust[v].clear();
//         // Now we clear the subgraph 
//         wilson_submap.clear_vertex(v);
//         // now we add all the neighbors not ignored 
//         for (auto const v_nbor: map_params.g[v])
//         {
//             if (!ignore[v_nbor]){
//                 wilson_submap.add_directed_edge(v, v_nbor);
//             }
//         }
        
//     }

    

//     // We also assume the tree has been properly cleared
//     // Now get a uniform spanning tree drawn on the subgraph denoted by the ignore
//     // vertices 
//     auto const result = sample_sub_ust<true>(map_params, ust, lower, upper, visited, ignore,
//                                 county_tree, county_stack, dummy_county_tree_queue,
//                                 county_pop, county_members,
//                                 c_visited, cty_pop_below, county_path, path, rng_state
//                                 );
//     // result.code == 0 means it was successful
//     return std::make_pair(result.code == 0, result.root);
// }
 

std::pair<bool, int> USTSampler::draw_ust(
      double const lower, double const upper,
        RNGState &rng_state) {
    // We assume that ignore has already been properly set 
    // We also assume the tree has been properly cleared
    // Now get a uniform spanning tree drawn on the subgraph denoted by the ignore
    // vertices 
    auto const result = sample_sub_ust<false>(map_params, ust, lower, upper, visited, ignore,
                                county_tree, mg_scratch.county_stack, g_scratch.dummy_county_tree_queue,
                                mg_scratch.county_pop, mg_scratch.county_members,
                                mg_scratch.c_visited, mg_scratch.cty_pop_below, 
                                mg_scratch.county_path, g_scratch.path, rng_state
                                );
    // result.code == 0 means it was successful
    return std::make_pair(result.code == 0, result.root);
}


USTDrawResult USTSampler::attempt_to_draw_tree_on_region(
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


USTDrawResult USTSampler::NEW_attempt_to_draw_tree_on_region(
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


USTDrawResult USTSampler::attempt_to_draw_tree_on_merged_region(RNGState &rng_state, Plan const &plan,
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

    auto prep_end_time = std::chrono::steady_clock::now();

    std::chrono::duration<double, std::ratio<1>> prep_time_diff = prep_end_time - prep_start_time;
    wilson_times.input_prep_time += prep_time_diff.count();

    std::pair<bool, int> result;

    auto sample_ust_start_time = std::chrono::steady_clock::now();
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