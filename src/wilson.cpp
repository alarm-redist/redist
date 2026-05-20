#include "wilson.h"
#include <array>

using LevelEdge = std::array<int, 3>;
using LevelPath = std::vector<LevelEdge>;
using ActiveMultigraph = std::vector<std::vector<LevelEdge>>;

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
void walk_until_cty(Multigraph &mg, int root,
                    std::vector<std::vector<int>> &path,
                    const std::vector<bool> &visited,
                    const std::vector<bool> &ignore);

void walk_until_level(Multigraph &mg, int root,
                      std::vector<std::vector<int>> &path,
                      const std::vector<bool> &visited);

void walk_until_active_level(ActiveMultigraph &mg, int root,
                             LevelPath &path,
                             const std::vector<bool> &visited);

ActiveMultigraph active_level_graph(const Graph &g, const std::vector<int> &level,
                                    const std::vector<bool> &ignore,
                                    const std::vector<int> &active_vertices,
                                    int n_group,
                                    const std::vector<int> &parent);

void loop_erase_level(LevelPath &path, int proposal, int root);

int rvtx_subset(const std::vector<int> &vertices,
                const std::vector<bool> &visited,
                int remaining, int &lower);

void split_active_components(const Graph &g,
                             const std::vector<uvec> &levels,
                             const std::vector<bool> &ignore,
                             const std::vector<int> &active_vertices,
                             std::vector<std::vector<int>> &active_levels,
                             std::vector<int> &n_groups);

int walk_until_label(const Graph &g, int root,
                     std::vector<int> &path, int MAX,
                     const std::vector<bool> &visited,
                     const std::vector<bool> &ignore,
                     const std::vector<int> &labels);

void prune_finest_groups(Tree &tree,
                         const std::vector<int> &active_vertices,
                         const std::vector<int> &finest_level,
                         int n_group,
                         int root,
                         std::vector<bool> &visited,
                         int &remaining,
                         const uvec &pop,
                         double lower, double upper);

/*
 * Erase loops in `path` that would be created by adding `proposal` to path
 */
// TESTED
void loop_erase_cty(std::vector<std::vector<int>> &path, int proposal, int root);


// [[Rcpp::export]]
Tree sample_ust(List l, const arma::uvec &pop, double lower, double upper,
                const arma::uvec &counties, const std::vector<bool> ignore) {
    Graph g = list_to_graph(l);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    Tree tree = init_tree(V);
    int root;
    std::vector<bool> visited(V);
    sample_sub_ust(g, tree, V, root, visited, ignore, pop, lower, upper, counties, cg);
    return tree;
}

/*
 * Sample a uniform spanning subtree of unvisited nodes using Wilson's algorithm
 */
// TESTED
int sample_sub_ust(const Graph &g, Tree &tree, int V, int &root,
                   std::vector<bool> &visited,
                   const std::vector<bool> &ignore, const uvec &pop,
                   double lower, double upper,
                   const uvec &counties, Multigraph &mg) {
    int n_county = mg.size();
    std::vector<bool> c_visited(n_county, true);
    uvec county_pop(n_county, fill::zeros);
    int tot_pop = 0;
    std::vector<std::vector<int>> county_members(n_county);
    int remaining = 0;
    for (int i = 0; i < V; i++) {
        if (ignore.at(i)) {
            visited[i] = true;
        } else {
            visited[i] = false;
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
    std::vector<std::vector<int>> path;
    while (c_remaining > 0) {
        int add = rvtx(c_visited, n_county, c_remaining, lower_c);
        // random walk from `add` until we hit the path
        walk_until_cty(mg, add, path, c_visited, ignore);
        // update visited list and constructed tree
        int added = path.size();
        if (added == 0) { // bail
            return 1;
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
        int max_try = std::max(50, 50 * remaining * ((int) std::log(remaining)));
        while (remaining > 0) {
            int add = rvtx(visited, V, remaining, lower_i);
            // random walk from `add` until we hit the path
            int added = walk_until(g, add, path, max_try, visited, ignore, counties);
            // update visited list and constructed tree
            if (added == 0) { // bail
                return 1;
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited.at(path[i]) = true;
                // reverse path so that arrows point away from root
                tree.at(path[i+1]).push_back(path[i]);
            }
        }
    }

    return 0;
}

/*
 * Sample a spanning subtree using nested levels ordered finest -> coarsest.
 * The coarsest graph is connected first, then progressively finer graphs are
 * connected within the parent group selected at the next coarser level. Once
 * each finest group has a root in the tree, Wilson walks fill in each group.
 */
int sample_sub_ust_hier(const Graph &g, Tree &tree, int V, int &root,
                        std::vector<bool> &visited,
                        const std::vector<bool> &ignore,
                        const std::vector<int> &active_vertices,
                        const uvec &pop,
                        double lower, double upper,
                        const std::vector<uvec> &levels) {
    int n_levels = levels.size();
    if (n_levels == 0) return 1;

    std::vector<std::vector<int>> active_levels;
    std::vector<int> n_groups;
    split_active_components(g, levels, ignore, active_vertices,
                            active_levels, n_groups);

    std::fill(visited.begin(), visited.end(), true);
    for (int v : active_vertices) visited[v] = false;
    int remaining = active_vertices.size();
    if (remaining == 0) return 1;

    int lower_i = 0;
    root = rvtx_subset(active_vertices, visited, remaining, lower_i);
    visited[root] = true;
    remaining--;

    LevelPath path;
    for (int lev = n_levels - 1; lev >= 0; lev--) {
        int n_group = n_groups[lev];
        std::vector<bool> group_active(n_group, false);
        std::vector<bool> group_has_root(n_group, false);
        std::vector<int> group_parent(n_group, -1);

        for (int i : active_vertices) {
            int grp = active_levels[lev][i] - 1;
            group_active[grp] = true;
            if (visited[i]) {
                group_has_root[grp] = true;
            }
            if (lev < n_levels - 1 && group_parent[grp] < 0) {
                group_parent[grp] = active_levels[lev + 1][i] - 1;
            }
        }

        ActiveMultigraph active_mg = active_level_graph(g, active_levels[lev], ignore,
                                                        active_vertices, n_group,
                                                        group_parent);

        if (lev == n_levels - 1) {
            std::vector<bool> group_visited(n_group, true);
            int group_remaining = 0;
            for (int grp = 0; grp < n_group; grp++) {
                group_visited[grp] = !group_active[grp] || group_has_root[grp];
                if (!group_visited[grp]) group_remaining++;
            }

            int lower_g = 0;
            while (group_remaining > 0) {
                int add = rvtx(group_visited, n_group, group_remaining, lower_g);

                walk_until_active_level(active_mg, add, path, group_visited);
                int added = path.size();
                if (added == 0) return 1;

                group_remaining -= added;
                group_visited[add] = true;
                for (int i = 0; i < added; i++) {
                    group_visited[path[i][0]] = true;
                    tree[path[i][2]].push_back(path[i][1]);
                    visited[path[i][1]] = true;
                    remaining--;
                }
            }
        } else {
            std::vector<bool> parent_seen(n_groups[lev + 1], false);
            for (int grp = 0; grp < n_group; grp++) {
                if (!group_active[grp]) continue;
                int parent_id = group_parent[grp];
                if (parent_id < 0 || parent_seen[parent_id]) continue;
                parent_seen[parent_id] = true;

                std::vector<bool> group_visited(n_group, true);
                int group_remaining = 0;
                bool has_root = false;
                for (int grp2 = 0; grp2 < n_group; grp2++) {
                    if (!group_active[grp2] || group_parent[grp2] != parent_id) {
                        group_visited[grp2] = true;
                        continue;
                    }
                    group_visited[grp2] = group_has_root[grp2];
                    if (group_has_root[grp2]) {
                        has_root = true;
                    } else {
                        group_remaining++;
                    }
                }
                if (!has_root) return 1;

                int lower_g = 0;
                while (group_remaining > 0) {
                    int add = rvtx(group_visited, n_group, group_remaining, lower_g);

                    walk_until_active_level(active_mg, add, path, group_visited);
                    int added = path.size();
                    if (added == 0) return 1;

                    group_remaining -= added;
                    group_visited[add] = true;
                    for (int i = 0; i < added; i++) {
                        group_visited[path[i][0]] = true;
                        tree[path[i][2]].push_back(path[i][1]);
                        visited[path[i][1]] = true;
                        remaining--;
                    }
                }
            }
        }
    }

    prune_finest_groups(tree, active_vertices, active_levels[0], n_groups[0],
                        root, visited, remaining, pop, lower, upper);

    if (remaining > 0) {
        std::vector<int> walk_path(remaining + 2);
        int max_try = std::max(50, 50 * remaining * ((int) std::log(remaining)));
        while (remaining > 0) {
            int add = rvtx_subset(active_vertices, visited, remaining, lower_i);
            int added = walk_until_label(g, add, walk_path, max_try, visited,
                                         ignore, active_levels[0]);
            if (added == 0) return 1;
            remaining -= added - 1;
            for (int i = 0; i < added - 1; i++) {
                visited[walk_path[i]] = true;
                tree[walk_path[i+1]].push_back(walk_path[i]);
            }
        }
    }

    return 0;
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
void walk_until_cty(Multigraph &mg, int root,
                    std::vector<std::vector<int>> &path,
                    const std::vector<bool> &visited,
                    const std::vector<bool> &ignore) {
    path.clear();

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
        path.clear();
    }
}

void walk_until_level(Multigraph &mg, int root,
                      std::vector<std::vector<int>> &path,
                      const std::vector<bool> &visited) {
    path.clear();

    int curr = root;
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        if (mg.at(curr).empty()) {
            path.clear();
            return;
        }
        int prop_idx = r_int((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        if (!visited.at(proposal)) {
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

void walk_until_active_level(ActiveMultigraph &mg, int root,
                             LevelPath &path,
                             const std::vector<bool> &visited) {
    path.clear();

    int curr = root;
    int i;
    int max = visited.size() * 500;
    for (i = 0; i < max; i++) {
        if (mg.at(curr).empty()) {
            path.clear();
            return;
        }
        int prop_idx = r_int((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        if (!visited.at(proposal)) {
            path.push_back(mg[curr][prop_idx]);
            loop_erase_level(path, proposal, root);
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

ActiveMultigraph active_level_graph(const Graph &g, const std::vector<int> &level,
                                    const std::vector<bool> &ignore,
                                    const std::vector<int> &active_vertices,
                                    int n_group,
                                    const std::vector<int> &parent) {
    ActiveMultigraph mg(n_group);
    for (int i : active_vertices) {
        int grp = level[i] - 1;
        for (int nbor : g[i]) {
            if (ignore[nbor]) continue;
            int nbor_grp = level[nbor] - 1;
            if (grp == nbor_grp) continue;
            if (parent[grp] >= 0 && parent[nbor_grp] != parent[grp]) continue;
            LevelEdge el = {nbor_grp, i, nbor};
            mg[grp].push_back(el);
        }
    }
    return mg;
}

void loop_erase_level(LevelPath &path, int proposal, int root) {
    int length = path.size();
    if (proposal == root) {
        path.erase(path.begin(), path.begin() + length);
        return;
    }

    int idx;
    for (idx = 0; idx < length - 1; idx++) {
        if (path[idx][0] == proposal) break;
    }

    if (idx != length - 1) {
        path.erase(path.begin() + idx + 1, path.begin() + length);
    }
}

int rvtx_subset(const std::vector<int> &vertices,
                const std::vector<bool> &visited,
                int remaining, int &lower) {
    int idx = r_int(remaining);
    int accuml = 0;
    bool seen_one = false;
    int last = vertices.back();
    for (int i = lower; i < (int) vertices.size(); i++) {
        int v = vertices[i];
        accuml += 1 - visited[v];
        if (!seen_one && !visited[v]) {
            seen_one = true;
            lower = i;
        }
        if (accuml - 1 == idx) return v;
    }
    return last;
}

void split_active_components(const Graph &g,
                             const std::vector<uvec> &levels,
                             const std::vector<bool> &ignore,
                             const std::vector<int> &active_vertices,
                             std::vector<std::vector<int>> &active_levels,
                             std::vector<int> &n_groups) {
    int n_levels = levels.size();
    int V = ignore.size();
    active_levels.assign(n_levels, std::vector<int>(V, 0));
    n_groups.assign(n_levels, 0);
    std::vector<int> stack;
    stack.reserve(active_vertices.size());

    for (int lev = 0; lev < n_levels; lev++) {
        std::vector<int> &label = active_levels[lev];
        if (lev > 0) {
            int max_raw = 0;
            for (int v : active_vertices) {
                int raw_label = (int) levels[lev](v);
                if (raw_label > max_raw) max_raw = raw_label;
            }
            std::vector<int> relabel(max_raw + 1, 0);
            for (int v : active_vertices) {
                int raw_label = (int) levels[lev](v);
                if (relabel[raw_label] == 0) {
                    relabel[raw_label] = ++n_groups[lev];
                }
                label[v] = relabel[raw_label];
            }
            continue;
        }

        for (int start : active_vertices) {
            if (label[start] != 0) continue;
            int base = levels[lev](start);
            int comp = ++n_groups[lev];
            label[start] = comp;
            stack.clear();
            stack.push_back(start);
            while (!stack.empty()) {
                int v = stack.back();
                stack.pop_back();
                for (int nbor : g[v]) {
                    if (ignore[nbor] || label[nbor] != 0) continue;
                    if ((int) levels[lev](nbor) != base) continue;
                    label[nbor] = comp;
                    stack.push_back(nbor);
                }
            }
        }
    }
}

int walk_until_label(const Graph &g, int root,
                     std::vector<int> &path, int MAX,
                     const std::vector<bool> &visited,
                     const std::vector<bool> &ignore,
                     const std::vector<int> &labels) {
    path[0] = root;
    int curr = root;
    int label = labels[root];
    int added = 1;
    int i;
    for (i = 0; i < MAX; i++) {
        int proposal = rnbor(g, curr);
        if (ignore[proposal] || labels[proposal] != label) {
            continue;
        } else if (!visited[proposal]) {
            for (int j = added - 1; j >= 0; j--) {
                if (path[j] == proposal) {
                    added = j;
                    break;
                }
            }
            path[added++] = proposal;
        } else {
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

void prune_finest_groups(Tree &tree,
                         const std::vector<int> &active_vertices,
                         const std::vector<int> &finest_level,
                         int n_group,
                         int root,
                         std::vector<bool> &visited,
                         int &remaining,
                         const uvec &pop,
                         double lower, double upper) {
    if (n_group <= 1 || remaining <= 0) return;

    Tree group_tree = init_tree(n_group);
    uvec group_pop(n_group, fill::zeros);
    std::vector<std::vector<int>> members(n_group);
    int total_pop = 0;

    for (int v : active_vertices) {
        int grp = finest_level[v] - 1;
        group_pop(grp) += pop[v];
        total_pop += pop[v];
        members[grp].push_back(v);
        for (int child : tree[v]) {
            int child_grp = finest_level[child] - 1;
            if (child_grp != grp) {
                group_tree[grp].push_back(child_grp);
            }
        }
    }

    std::vector<int> group_pop_below(n_group, -1);
    std::vector<int> group_parent(n_group);
    int group_root = finest_level[root] - 1;
    tree_pop(group_tree, group_root, group_pop, group_pop_below, group_parent);

    for (int grp = 0; grp < n_group; grp++) {
        int n_vtx = members[grp].size();
        if (n_vtx <= 1) continue;

        int split_ub = group_pop_below[grp];
        if (split_ub < 0) continue;
        int split_lb = split_ub - group_pop[grp];
        if (lower - 1 < group_pop[grp]) split_lb = (int) lower;

        for (int child_grp : group_tree[grp]) {
            int pop_child = group_pop_below[child_grp];
            if (pop_child >= 0 && pop_child < split_lb) {
                split_lb = pop_child;
            }
        }

        bool miss_first = split_ub < lower || split_lb > upper;
        bool miss_second = (total_pop - split_lb) < lower ||
                           (total_pop - split_ub) > upper;
        if (!(miss_first && miss_second)) continue;

        int group_root_idx = -1;
        for (int j = 0; j < n_vtx; j++) {
            if (visited[members[grp][j]]) {
                group_root_idx = j;
                break;
            }
        }
        if (group_root_idx < 0) continue;

        for (int j = 0; j < n_vtx; j++) {
            int vtx_idx = members[grp][j];
            if (!visited[vtx_idx]) {
                remaining--;
            }
            if (j > 0 && j != group_root_idx + 1) {
                tree[vtx_idx].push_back(members[grp][j - 1]);
            }
            visited[vtx_idx] = true;
        }

        if (group_root_idx < n_vtx - 1) {
            tree[members[grp][group_root_idx]].push_back(members[grp][n_vtx - 1]);
        }
    }
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
