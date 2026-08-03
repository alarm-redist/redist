#include "wilson.h"

int rvtx_subset(const std::vector<int> &vertices,
                const std::vector<bool> &visited,
                int &lower);

class Frontier {
  public:
    explicit Frontier(int size) : queued_(size, false) {}

    void add(int vertex, const std::vector<bool> &visited) {
        if (!visited[vertex] && !queued_[vertex]) {
            vertices_.push_back(vertex);
            queued_[vertex] = true;
        }
    }

    void clear() {
        vertices_.clear();
        std::fill(queued_.begin(), queued_.end(), false);
        next_vertex_ = 0;
    }

    int next(const std::vector<bool> &visited, int size, int &lower) {
        while (next_vertex_ < static_cast<int>(vertices_.size()) &&
               visited[vertices_[next_vertex_]]) {
            next_vertex_++;
        }
        if (next_vertex_ < static_cast<int>(vertices_.size())) {
            return vertices_[next_vertex_++];
        }
        return rvtx(visited, size, lower);
    }

    int next_subset(const std::vector<int> &vertices,
                    const std::vector<bool> &visited, int &lower) {
        while (next_vertex_ < static_cast<int>(vertices_.size()) &&
               visited[vertices_[next_vertex_]]) {
            next_vertex_++;
        }
        if (next_vertex_ < static_cast<int>(vertices_.size())) {
            return vertices_[next_vertex_++];
        }
        return rvtx_subset(vertices, visited, lower);
    }

  private:
    std::vector<int> vertices_;
    std::vector<bool> queued_;
    int next_vertex_ = 0;
};

static void add_county_frontier(const Multigraph &mg, int county,
                                const std::vector<bool> &visited,
                                const std::vector<bool> &ignore,
                                Frontier &frontier) {
    for (const std::vector<int> &edge : mg[county]) {
        if (!ignore[edge[1]] && !ignore[edge[2]]) {
            frontier.add(edge[0], visited);
        }
    }
}

static void add_active_level_frontier(const ActiveMultigraph &mg, int group,
                                      const std::vector<bool> &visited,
                                      Frontier &frontier) {
    for (const LevelEdge &edge : mg[group]) {
        frontier.add(edge[0], visited);
    }
}

static void add_precinct_frontier(const Graph &g, int vertex,
                                  const std::vector<bool> &visited,
                                  const std::vector<bool> &ignore,
                                  const uvec &counties, Frontier &frontier) {
    if (ignore[vertex]) {
        return;
    }
    for (int neighbor : g[vertex]) {
        if (!ignore[neighbor] && counties[neighbor] == counties[vertex]) {
            frontier.add(neighbor, visited);
        }
    }
}

static void add_flat_frontier(const std::vector<int> &flat_adj,
                              const std::vector<int> &flat_off,
                              int vertex,
                              const std::vector<bool> &visited,
                              const std::vector<bool> &ignore,
                              Frontier &frontier) {
    if (ignore[vertex]) return;
    for (int i = flat_off[vertex]; i < flat_off[vertex + 1]; i++) {
        int neighbor = flat_adj[i];
        if (!ignore[neighbor]) frontier.add(neighbor, visited);
    }
}

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root,
               std::vector<int> &path, int MAX,
               std::vector<int> &next_vertex,
               const std::vector<bool> &visited,
               const std::vector<bool> &ignore,
               const uvec &counties);

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(Multigraph &mg, int root,
                    std::vector<std::vector<int>> &path,
                    std::vector<int> &next_group,
                    std::vector<int> &next_edge,
                    const std::vector<bool> &visited,
                    const std::vector<bool> &ignore);

void walk_until_level(Multigraph &mg, int root,
                      std::vector<std::vector<int>> &path,
                      std::vector<int> &next_group,
                      std::vector<int> &next_edge,
                      const std::vector<bool> &visited);

void walk_until_active_level(ActiveMultigraph &mg, int root,
                             LevelPath &path,
                             std::vector<int> &next_group,
                             std::vector<int> &next_edge,
                             const std::vector<bool> &visited);

void active_level_graph(
                        const Graph &g,
                        const std::vector<int> &level,
                        const std::vector<bool> &ignore,
                        const std::vector<int> &active_vertices,
                        int n_group,
                        const std::vector<int> &parent,
                        ActiveMultigraph &mg);

int rvtx_subset(const std::vector<int> &vertices,
                const std::vector<bool> &visited,
                int &lower);

int walk_until_flat(const std::vector<int> &flat_adj,
                    const std::vector<int> &flat_off,
                    int root,
                    std::vector<int> &path,
                    int MAX,
                    std::vector<std::int8_t> &status,
                    std::vector<int> &path_pos);

void prune_finest_groups(Tree &tree,
                         const std::vector<int> &active_vertices,
                         const std::vector<int> &finest_level,
                         int n_group,
                         int root,
                         std::vector<bool> &visited,
                         int &remaining,
                         const uvec &pop,
                         double lower, double upper,
                         HierarchicalSamplerWorkspace &workspace);

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
    root = rvtx(visited, V, lower_i);
    visited[root] = true;
    remaining--;
    c_visited.at(counties[root] - 1) = true;
    c_remaining--;

    // Connect counties
    Tree cty_tree = init_tree(n_county);
    std::vector<std::vector<int>> path;
    std::vector<int> next_group(n_county, -1);
    std::vector<int> next_edge(n_county, -1);
    Frontier county_frontier(n_county);
    add_county_frontier(mg, counties[root] - 1, c_visited, ignore, county_frontier);
    while (c_remaining > 0) {
        int add = county_frontier.next(c_visited, n_county, lower_c);
        // random walk from `add` until we hit the path
        walk_until_cty(mg, add, path, next_group, next_edge, c_visited, ignore);
        // update visited list and constructed tree
        int added = path.size();
        if (added == 0) { // bail
            return 1;
        }
        c_remaining -= added;
        c_visited.at(add) = true;
        add_county_frontier(mg, add, c_visited, ignore, county_frontier);
        for (int i = 0; i < added; i++) {
            c_visited.at(path[i][0]) = true;
            add_county_frontier(mg, path[i][0], c_visited, ignore, county_frontier);
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
        std::vector<int> next_vertex(V, -1);
        Frontier precinct_frontier(V);
        for (int i = 0; i < V; i++) {
            if (visited[i]) {
                add_precinct_frontier(g, i, visited, ignore, counties, precinct_frontier);
            }
        }
        int max_try = std::max(50, 50 * remaining * ((int) std::log(remaining)));
        while (remaining > 0) {
            int add = precinct_frontier.next(visited, V, lower_i);
            // random walk from `add` until we hit the path
            int added = walk_until(g, add, path, max_try, next_vertex,
                                   visited, ignore, counties);
            // update visited list and constructed tree
            if (added == 0) { // bail
                return 1;
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited.at(path[i]) = true;
                add_precinct_frontier(g, path[i], visited, ignore, counties,
                                      precinct_frontier);
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
                        const std::vector<std::vector<int>> &levels,
                        const std::vector<int> &n_groups,
                        const std::vector<std::vector<int>> &parents,
                        HierarchicalSamplerWorkspace &workspace,
                        const std::vector<int> &finest_adj,
                        const std::vector<int> &finest_off) {
    int n_levels = levels.size();
    if (n_levels == 0 || n_groups.size() != n_levels ||
        parents.size() != n_levels) {
        return 1;
    }

    std::fill(visited.begin(), visited.end(), true);
    for (int v : active_vertices) visited[v] = false;
    int remaining = active_vertices.size();
    if (remaining == 0) return 1;

    int lower_i = 0;
    root = rvtx_subset(active_vertices, visited, lower_i);
    visited[root] = true;
    remaining--;

    LevelPath &path = workspace.level_path;
    path.clear();
    for (int lev = n_levels - 1; lev >= 0; lev--) {
        int n_group = n_groups[lev];
        workspace.group_active.assign(n_group, false);
        workspace.group_has_root.assign(n_group, false);
        const std::vector<int> &group_parent = parents[lev];

        for (int i : active_vertices) {
            int grp = levels[lev][i] - 1;
            workspace.group_active[grp] = true;
            if (visited[i]) {
                workspace.group_has_root[grp] = true;
            }
        }

        active_level_graph(g, levels[lev], ignore, active_vertices, n_group,
                           group_parent, workspace.active_graph);
        workspace.next_group.assign(n_group, -1);
        workspace.next_edge.assign(n_group, -1);

        if (lev == n_levels - 1) {
            workspace.group_visited.assign(n_group, true);
            std::vector<bool> &group_visited = workspace.group_visited;
            int group_remaining = 0;
            for (int grp = 0; grp < n_group; grp++) {
                group_visited[grp] = !workspace.group_active[grp] ||
                                     workspace.group_has_root[grp];
                if (!group_visited[grp]) group_remaining++;
            }

            Frontier group_frontier(n_group);
            for (int grp = 0; grp < n_group; grp++) {
                if (group_visited[grp]) {
                    add_active_level_frontier(workspace.active_graph, grp,
                                              group_visited,
                                              group_frontier);
                }
            }
            int lower_g = 0;
            while (group_remaining > 0) {
                int add = group_frontier.next(group_visited, n_group, lower_g);

                walk_until_active_level(workspace.active_graph, add, path,
                                        workspace.next_group, workspace.next_edge,
                                        group_visited);
                int added = path.size();
                if (added == 0) return 1;

                group_remaining -= added;
                group_visited[add] = true;
                add_active_level_frontier(workspace.active_graph, add,
                                          group_visited,
                                           group_frontier);
                for (int i = 0; i < added; i++) {
                    group_visited[path[i][0]] = true;
                    add_active_level_frontier(workspace.active_graph, path[i][0],
                                              group_visited,
                                              group_frontier);
                    tree[path[i][2]].push_back(path[i][1]);
                    visited[path[i][1]] = true;
                    remaining--;
                }
            }
        } else {
            workspace.parent_seen.assign(n_groups[lev + 1], false);
            std::vector<bool> &parent_seen = workspace.parent_seen;
            Frontier group_frontier(n_group);
            for (int grp = 0; grp < n_group; grp++) {
                if (!workspace.group_active[grp]) continue;
                int parent_id = group_parent[grp];
                if (parent_id < 0 || parent_seen[parent_id]) continue;
                parent_seen[parent_id] = true;

                workspace.group_visited.assign(n_group, true);
                std::vector<bool> &group_visited = workspace.group_visited;
                int group_remaining = 0;
                bool has_root = false;
                for (int grp2 = 0; grp2 < n_group; grp2++) {
                    if (!workspace.group_active[grp2] ||
                        group_parent[grp2] != parent_id) {
                        group_visited[grp2] = true;
                        continue;
                    }
                    group_visited[grp2] = workspace.group_has_root[grp2];
                    if (workspace.group_has_root[grp2]) {
                        has_root = true;
                    } else {
                        group_remaining++;
                    }
                }
                if (!has_root) return 1;

                group_frontier.clear();
                for (int grp2 = 0; grp2 < n_group; grp2++) {
                    if (workspace.group_has_root[grp2]) {
                        add_active_level_frontier(workspace.active_graph, grp2,
                                                  group_visited, group_frontier);
                    }
                }
                int lower_g = 0;
                while (group_remaining > 0) {
                    int add = group_frontier.next(group_visited, n_group, lower_g);

                    walk_until_active_level(workspace.active_graph, add, path,
                                            workspace.next_group,
                                            workspace.next_edge, group_visited);
                    int added = path.size();
                    if (added == 0) return 1;

                    group_remaining -= added;
                    group_visited[add] = true;
                    add_active_level_frontier(workspace.active_graph, add,
                                              group_visited, group_frontier);
                    for (int i = 0; i < added; i++) {
                        group_visited[path[i][0]] = true;
                        add_active_level_frontier(workspace.active_graph,
                                                  path[i][0], group_visited,
                                                  group_frontier);
                        tree[path[i][2]].push_back(path[i][1]);
                        visited[path[i][1]] = true;
                        remaining--;
                    }
                }
            }
        }
    }

    prune_finest_groups(tree, active_vertices, levels[0], n_groups[0],
                        root, visited, remaining, pop, lower, upper,
                        workspace);

    if (remaining > 0) {
        workspace.walk_path.resize(remaining + 2);
        workspace.path_pos.assign(V, -1);
        workspace.status.assign(V, 2);
        for (int vertex : active_vertices) {
            workspace.status[vertex] = visited[vertex] ? 1 : 0;
        }
        Frontier finest_frontier(V);
        for (int vertex : active_vertices) {
            if (visited[vertex]) {
                add_flat_frontier(finest_adj, finest_off, vertex, visited,
                                  ignore, finest_frontier);
            }
        }
        int max_try = std::max(50, 50 * remaining * ((int) std::log(remaining)));
        while (remaining > 0) {
            int add = finest_frontier.next_subset(active_vertices, visited, lower_i);
            int added = walk_until_flat(finest_adj, finest_off, add,
                                        workspace.walk_path, max_try,
                                        workspace.status, workspace.path_pos);
            if (added == 0) return 1;
            remaining -= added - 1;
            for (int i = 0; i < added - 1; i++) {
                visited[workspace.walk_path[i]] = true;
                workspace.status[workspace.walk_path[i]] = 1;
                add_flat_frontier(finest_adj, finest_off,
                                  workspace.walk_path[i], visited, ignore,
                                  finest_frontier);
                tree[workspace.walk_path[i+1]].push_back(workspace.walk_path[i]);
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
               std::vector<int> &next_vertex,
               const std::vector<bool> &visited,
               const std::vector<bool> &ignore,
               const uvec &counties) {
    int curr = root;
    int county = counties[root];
    bool hit_tree = false;
    for (int i = 0; i < MAX; i++) {
        int proposal = rnbor(g, curr);
        if (ignore[proposal] || counties[proposal] != county) {
            continue;
        }

        next_vertex[curr] = proposal;
        if (visited[proposal]) {
            hit_tree = true;
            break;
        }
        curr = proposal;
    }

    if (!hit_tree) {
        return 0;
    }

    int added = 0;
    curr = root;
    while (!visited[curr]) {
        path[added++] = curr;
        curr = next_vertex[curr];
    }
    path[added++] = curr;
    return added;
}



/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(Multigraph &mg, int root,
                    std::vector<std::vector<int>> &path,
                    std::vector<int> &next_group,
                    std::vector<int> &next_edge,
                    const std::vector<bool> &visited,
                    const std::vector<bool> &ignore) {
    path.clear();

    int curr = root;
    bool hit_tree = false;
    int max = visited.size() * 500;
    for (int i = 0; i < max; i++) {
        if (mg.at(curr).empty()) {
            break;
        }
        int prop_idx = r_int((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        if (ignore[mg[curr][prop_idx][2]] || ignore[mg[curr][prop_idx][1]]) {
            continue;
        }

        next_group[curr] = proposal;
        next_edge[curr] = prop_idx;
        if (visited.at(proposal)) {
            hit_tree = true;
            break;
        }
        curr = proposal;
    }

    if (!hit_tree) {
        return;
    }

    curr = root;
    while (!visited[curr]) {
        path.push_back(mg[curr][next_edge[curr]]);
        curr = next_group[curr];
    }
}

void walk_until_level(Multigraph &mg, int root,
                      std::vector<std::vector<int>> &path,
                      std::vector<int> &next_group,
                      std::vector<int> &next_edge,
                      const std::vector<bool> &visited) {
    path.clear();

    int curr = root;
    bool hit_tree = false;
    int max = visited.size() * 500;
    for (int i = 0; i < max; i++) {
        if (mg.at(curr).empty()) {
            return;
        }
        int prop_idx = r_int((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        next_group[curr] = proposal;
        next_edge[curr] = prop_idx;
        if (visited.at(proposal)) {
            hit_tree = true;
            break;
        }
        curr = proposal;
    }

    if (!hit_tree) {
        return;
    }

    curr = root;
    while (!visited[curr]) {
        path.push_back(mg[curr][next_edge[curr]]);
        curr = next_group[curr];
    }
}

void walk_until_active_level(ActiveMultigraph &mg, int root,
                             LevelPath &path,
                             std::vector<int> &next_group,
                             std::vector<int> &next_edge,
                             const std::vector<bool> &visited) {
    path.clear();

    int curr = root;
    bool hit_tree = false;
    int max = visited.size() * 500;
    for (int i = 0; i < max; i++) {
        if (mg.at(curr).empty()) {
            return;
        }
        int prop_idx = r_int((int) mg.at(curr).size());
        int proposal = mg[curr][prop_idx][0];
        next_group[curr] = proposal;
        next_edge[curr] = prop_idx;
        if (visited.at(proposal)) {
            hit_tree = true;
            break;
        }
        curr = proposal;
    }

    if (!hit_tree) {
        return;
    }

    curr = root;
    while (!visited[curr]) {
        path.push_back(mg[curr][next_edge[curr]]);
        curr = next_group[curr];
    }
}

void active_level_graph(const Graph &g, const std::vector<int> &level,
                        const std::vector<bool> &ignore,
                        const std::vector<int> &active_vertices,
                        int n_group,
                        const std::vector<int> &parent,
                        ActiveMultigraph &mg) {
    mg.resize(n_group);
    for (std::vector<LevelEdge> &group_edges : mg) group_edges.clear();
    for (int i : active_vertices) {
        int grp = level[i] - 1;
        for (int nbor : g[i]) {
            if (ignore[nbor]) continue;
            int nbor_grp = level[nbor] - 1;
            if (grp == nbor_grp) continue;
            if (parent[grp] >= 0 && parent[nbor_grp] != parent[grp]) continue;
            LevelEdge edge = {nbor_grp, i, nbor};
            mg[grp].push_back(edge);
        }
    }
}

int rvtx_subset(const std::vector<int> &vertices,
                const std::vector<bool> &visited,
                int &lower) {
    for (int i = lower; i < (int) vertices.size(); i++) {
        int v = vertices[i];
        if (!visited[v]) {
            lower = i;
            return v;
        }
    }
    return -1;
}

int walk_until_flat(const std::vector<int> &flat_adj,
                    const std::vector<int> &flat_off,
                    int root,
                    std::vector<int> &path,
                    int MAX,
                    std::vector<std::int8_t> &status,
                    std::vector<int> &path_pos) {
    path[0] = root;
    path_pos[root] = 0;
    int curr = root;
    int added = 1;
    int i;
    for (i = 0; i < MAX; i++) {
        int off = flat_off[curr];
        int degree = flat_off[curr + 1] - off;
        if (degree == 0) return 0;
        int proposal = flat_adj[off + r_int(degree)];
        std::int8_t proposal_status = status[proposal];
        if (proposal_status == 2) {
            continue;
        } else if (proposal_status == 0) {
            int position = path_pos[proposal];
            if (position >= 0 && position < added && path[position] == proposal) {
                added = position;
            }
            path[added] = proposal;
            path_pos[proposal] = added++;
        } else {
            path[added++] = proposal;
            break;
        }
        curr = proposal;
    }
    if (i == MAX) return 0;
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
                         double lower, double upper,
                         HierarchicalSamplerWorkspace &workspace) {
    if (n_group <= 1 || remaining <= 0) return;

    Tree &group_tree = workspace.prune_group_tree;
    if (static_cast<int>(group_tree.size()) != n_group) {
        group_tree = init_tree(n_group);
    } else {
        for (std::vector<int> &children : group_tree) children.clear();
    }
    uvec &group_pop = workspace.prune_group_pop;
    group_pop.set_size(n_group);
    group_pop.zeros();
    std::vector<std::vector<int>> &members = workspace.prune_members;
    members.resize(n_group);
    for (std::vector<int> &group_members : members) group_members.clear();
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

    std::vector<int> &group_pop_below = workspace.prune_group_pop_below;
    group_pop_below.assign(n_group, -1);
    std::vector<int> &group_parent = workspace.prune_group_parent;
    group_parent.resize(n_group);
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
