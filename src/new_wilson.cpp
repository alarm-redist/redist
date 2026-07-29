#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "advanced_types.h"
#include "random.h"

namespace{

// idea to prepare is 
// walk through the graph and only add 

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

}