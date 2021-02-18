#include "wilson.h"
#include "tree_op.h"

/*
 * Sample a uniform spanning tree using Wilson's algorithm
 */
// TESTED
Tree sample_ust(List l, int &root, const uvec &counties) {
    Graph g = list_to_graph(l);
    int V = g.size();
    Tree tree = init_tree(V);
    Multigraph mg = county_graph(g, counties);

    // initialize 'hit' list with sampled vertex
    std::vector<bool> ignore (V, false);
    return sample_sub_ust(g, tree, V, root, ignore, counties, mg);
}

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
// TESTED
Tree sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                    const std::vector<bool> &ignore,
                    const uvec &counties, Multigraph &mg) {
    int n_county = mg.size();
    std::vector<bool> visited(V, false);
    std::vector<bool> c_visited(n_county, true);
    int remaining = 0;
    for (int i = 0; i < V; i++) {
        if (ignore.at(i)) {
            visited[i] = true;
        } else {
            remaining++;
            c_visited.at(counties(i) - 1) = false;
        }
    }

    int c_remaining = 0;
    for (int i = 0; i < n_county; i++) {
        c_remaining += 1 - c_visited.at(i);
    }

    // pick root
    root = rvtx(visited, V, remaining);
    visited[root] = true;
    remaining--;
    c_visited.at(counties[root] - 1) = true;
    c_remaining--;

    // Connect counties
    while (c_remaining > 0) {
        int add = rvtx(c_visited, n_county, c_remaining);
        // random walk from `add` until we hit the path
        std::vector<std::vector<int>> path = walk_until_cty(mg, add,
                                                            c_visited, ignore);
        // update visited list and constructed tree
        int added = path.size();
        if (added == 0) { // bail
            Tree null_tree;
            return null_tree;
        }
        c_remaining -= added;
        c_visited.at(add) = true;
        for (int i = 0; i < added; i++) {
            c_visited.at(path[i][0]) = true;
            tree.at(path[i][2]).push_back(path[i][1]);
            visited.at(path[i][1]) = true; // root for next district
            remaining--;
        }
    }

    // Generate tree within each county
    do {
        int add = rvtx(visited, V, remaining);
        // random walk from `add` until we hit the path
        std::vector<int> path = walk_until(g, add, visited, ignore, counties);
        // update visited list and constructed tree
        int added = path.size();
        if (added == 0) { // bail
            Tree null_tree;
            return null_tree;
        }
        remaining -= added - 1; // minus 1 because ending vertex already in tree
        for (int i = 0; i < added - 1; i++) {
            visited.at(path[i]) = true;
            tree.at(path[i+1]).push_back(path[i]);
        }
    } while (remaining > 0);

    return tree;
}



/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
std::vector<int> walk_until(const Graph &g, int root,
                            const std::vector<bool> &visited,
                            const std::vector<bool> &ignore,
                            const uvec &counties) {
    std::vector<int> path = {root};
    // walk until we hit something in `visited`
    int curr = root;
    int county = counties(root);
    //while (true) {
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        int proposal = rnbor(g, curr);
        if (ignore.at(proposal) || counties(proposal) != county) {
            continue;
        } else if (!visited.at(proposal)) {
            loop_erase(path, proposal);
            path.push_back(proposal);
        } else {
            path.push_back(proposal);
            break;
        }
        curr = proposal;
    }
    if (i == max) {
        path.erase(path.begin(), path.end());
    }

    return path;
}


/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase(std::vector<int> &path, int proposal) {
    int idx;
    int length = path.size();
    for (idx = 0; idx < length; idx++) {
        if (path[idx] == proposal) break;
    }

    if (idx != length) { // a loop
        path.erase(path.begin() + idx, path.begin() + length);
    }
}


/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
std::vector<std::vector<int>> walk_until_cty(Multigraph &mg, int root,
                                             const std::vector<bool> &visited,
                                             const std::vector<bool> &ignore) {
    std::vector<std::vector<int>> path;

    // walk until we hit something in `visited`
    int curr = root;
    //while (true) {
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        int prop_idx = rint((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        if (ignore[mg[curr][prop_idx][2]] || ignore[mg[curr][prop_idx][1]]) {
            continue;
        } else if (!visited.at(proposal)) {
            path.push_back(mg[curr][prop_idx]);
            loop_erase_cty(path, proposal, root);
        } else {
            path.push_back(mg[curr][prop_idx]);
            break;
        }
        curr = proposal;
    }
    if (i == max) {
        path.erase(path.begin(), path.end());
    }

    return path;
}


/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::vector<int>> &path, int proposal, int root) {
    int length = path.size();
    if (proposal == root) {
        path.erase(path.begin(), path.begin() + length);
        return;
    }

    int idx;
    for (idx = 0; idx < length - 1; idx++) {
        if (path[idx][0] == proposal) break;
    }

    if (idx != length - 1) { // a loop
        path.erase(path.begin() + idx + 1, path.begin() + length);
    }
}
