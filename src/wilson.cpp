#include "wilson.h"

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root,
               std::vector<int> &path, int MAX,
               const std::vector<bool> &visited,
               const std::vector<bool> &ignore,
               const uvec &counties);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
void loop_erase(std::vector<int> &path, int proposal);

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
std::vector<std::vector<int>> walk_until_cty(Multigraph &mg, int root,
                                             const std::vector<bool> &visited,
                                             const std::vector<bool> &ignore);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::vector<int>> &path, int proposal, int root);


// [[Rcpp::export]]
Tree sample_ust(List l, const arma::uvec &pop, double lower, double upper,
                const arma::uvec &counties) {
    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    Tree tree = init_tree(V);
    int root;
    const std::vector<bool> ignore(V, false);
    return sample_sub_ust(g, tree, V, root, ignore, pop, lower, upper, counties, cg);
}

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
// TESTED
Tree sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                    const std::vector<bool> &ignore, const uvec &pop,
                    double lower, double upper,
                    const uvec &counties, Multigraph &mg) {
    int n_county = mg.size();
    std::vector<bool> visited(V, false);
    std::vector<bool> c_visited(n_county, true);
    uvec county_pop(n_county, fill::zeros);
    int tot_pop = 0;
    std::vector<std::vector<int>> county_members(n_county);
    int remaining = 0;
    for (int i = 0; i < V; i++) {
        if (ignore.at(i)) {
            visited[i] = true;
        } else {
            remaining++;
            int county = counties(i) - 1;
            tot_pop += pop(i);
            county_pop(county) += pop[i];
            if (c_visited[county]) {
                c_visited[county] = false;
                county_members[county] = std::vector<int>();
                county_members[county].reserve(16);
            }
            county_members[county].push_back(i);
        }
    }

    int c_remaining = 0;
    for (int i = 0; i < n_county; i++) {
        c_remaining += 1 - c_visited.at(i);
    }

    // pick root
    int lower_i = 0;
    int lower_c = 0;
    root = rvtx(visited, V, remaining, lower_i);
    visited[root] = true;
    remaining--;
    c_visited.at(counties[root] - 1) = true;
    c_remaining--;

    // Connect counties
    Tree cty_tree = init_tree(n_county);
    while (c_remaining > 0) {
        int add = rvtx(c_visited, n_county, c_remaining, lower_c);
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
            // reverse path so that arrows point away from root
            tree.at(path[i][2]).push_back(path[i][1]);
            cty_tree.at(path[i][0]).push_back(counties(path[i][1])-1);

            visited.at(path[i][1]) = true; // root for next district
            remaining--;
        }
    }

    // figure out which counties will not need to be split
    if (n_county > 1) {
    std::vector<int> cty_pop_below(n_county, -1);
    std::vector<int> cty_parent(n_county);
    tree_pop(cty_tree, counties[root] - 1, county_pop, cty_pop_below, cty_parent);
    for (int i = 0; i < n_county; i++) {
        int n_vtx = county_members[i].size();
        if (n_vtx <= 1) continue;
        // check child counties
        int children = cty_tree[i].size();
        int split_ub = cty_pop_below[i];
        int split_lb = split_ub - county_pop[i];
        if (lower-1 <  county_pop[i]) split_lb = (int) lower;
        for (int j = 0; j < children; j++) {
            int pop_child = cty_pop_below[cty_tree[i][j]];
            if (pop_child >= 0 && pop_child < split_lb) {
                split_lb = pop_child;
            }
        }
        // whether the range of split populations misses the 3 possible target intervals
        bool miss_first = split_ub < lower || split_lb > upper;
        bool miss_second = (tot_pop - split_lb) < lower || (tot_pop - split_ub) > upper;

        // impossible for this county to need to be split
        if (cty_pop_below[i] >= 0 && (miss_first && miss_second)) {
            // fill in with a dummy tree
            remaining -= n_vtx - 1; // already visited county root
            int cty_root = -1;
            for (int j = 0; j < n_vtx; j++) {
                int vtx_idx = county_members[i][j];
                if (visited.at(vtx_idx)) { // county root
                    cty_root = j;
                }
                if (j > 0 && j != cty_root + 1) {
                    tree.at(vtx_idx).push_back(county_members[i][j-1]);
                }
                visited.at(vtx_idx) = true;
            }

            if (cty_root < n_vtx - 1) {
                tree.at(county_members[i][cty_root]).push_back(county_members[i][n_vtx-1]);
            }
        }
    }
    }

    // Generate tree within each county
    if (remaining > 0) {
        std::vector<int> path(remaining + 2);
        int max_try = 50 * remaining * ((int) std::log(remaining));
        while (remaining > 0) {
            int add = rvtx(visited, V, remaining, lower_i);
            // random walk from `add` until we hit the path
            int added = walk_until(g, add, path, max_try, visited, ignore, counties);
            // update visited list and constructed tree
            if (added == 0) { // bail
                Tree null_tree;
                return null_tree;
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited.at(path[i]) = true;
                // reverse path so that arrows point away from root
                tree.at(path[i+1]).push_back(path[i]);
            }
        }
    }

    return tree;
}



/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root,
               std::vector<int> &path, int MAX,
               const std::vector<bool> &visited,
               const std::vector<bool> &ignore,
               const uvec &counties) {
    path[0] = root;
    // walk until we hit something in `visited`
    int curr = root;
    int county = counties[root];
    int added = 1; // cursor
    int i;
    for (i = 0; i < MAX; i++) {
        int proposal = rnbor(g, curr);
        if (ignore[proposal] || counties[proposal] != county) {
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
        int prop_idx = r_int((int) mg.at(curr).size());
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
