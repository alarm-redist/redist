#include "wilson.h"

namespace {

class Frontier {
  public:
    explicit Frontier(int size) : queued_(size, false) {}

    void add(int vertex, const std::vector<bool> &visited) {
        if (!visited[vertex] && !queued_[vertex]) {
            vertices_.push_back(vertex);
            queued_[vertex] = true;
        }
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

  private:
    std::vector<int> vertices_;
    std::vector<bool> queued_;
    int next_vertex_ = 0;
};

void add_precinct_frontier(const Graph &g, int vertex, const std::vector<bool> &visited,
                           const std::vector<bool> &ignore, const uvec &counties,
                           Frontier &frontier) {
    if (ignore[vertex]) {
        return;
    }

    for (int neighbor : g[vertex]) {
        if (!ignore[neighbor] && counties[neighbor] == counties[vertex]) {
            frontier.add(neighbor, visited);
        }
    }
}

void add_county_frontier(const Multigraph &mg, int county,
                         const std::vector<bool> &visited,
                         const std::vector<bool> &ignore, Frontier &frontier) {
    for (const std::vector<int> &edge : mg[county]) {
        if (!ignore[edge[1]] && !ignore[edge[2]]) {
            frontier.add(edge[0], visited);
        }
    }
}

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
int walk_until(const Graph &g, int root, std::vector<int> &path, int MAX,
               const std::vector<bool> &visited, const std::vector<bool> &ignore,
               const uvec &counties, std::vector<int> &next_vertex);

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
void walk_until_cty(Multigraph &mg, int root, std::vector<std::vector<int>> &path,
                    const std::vector<bool> &visited, const std::vector<bool> &ignore,
                    std::vector<int> &next_vertex, std::vector<int> &next_edge);

} // namespace

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
int sample_sub_ust(const Graph &g, Tree &tree, int V, int &root, std::vector<bool> &visited,
                   const std::vector<bool> &ignore, const uvec &pop, double lower, double upper,
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
    std::vector<int> cty_next_vertex(n_county, -1);
    std::vector<int> cty_next_edge(n_county, -1);
    Frontier cty_frontier(n_county);
    add_county_frontier(mg, counties[root] - 1, c_visited, ignore, cty_frontier);
    while (c_remaining > 0) {
        int add = cty_frontier.next(c_visited, n_county, lower_c);
        // random walk from `add` until we hit the path
        walk_until_cty(mg, add, path, c_visited, ignore, cty_next_vertex, cty_next_edge);
        // update visited list and constructed tree
        int added = path.size();
        if (added == 0) { // bail
            return 1;
        }
        c_remaining -= added;
        c_visited.at(add) = true;
        add_county_frontier(mg, add, c_visited, ignore, cty_frontier);
        for (int i = 0; i < added; i++) {
            c_visited.at(path[i][0]) = true;
            add_county_frontier(mg, path[i][0], c_visited, ignore, cty_frontier);
            // reverse path so that arrows point away from root
            tree.at(path[i][2]).push_back(path[i][1]);
            cty_tree.at(path[i][0]).push_back(counties(path[i][1]) - 1);

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
            if (n_vtx <= 1)
                continue;
            // check child counties
            int children = cty_tree[i].size();
            int split_ub = cty_pop_below[i];
            int split_lb = split_ub - county_pop[i];
            if (lower - 1 < county_pop[i])
                split_lb = (int)lower;
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
                        tree.at(vtx_idx).push_back(county_members[i][j - 1]);
                    }
                    visited.at(vtx_idx) = true;
                }

                if (cty_root < n_vtx - 1) {
                    tree.at(county_members[i][cty_root])
                        .push_back(county_members[i][n_vtx - 1]);
                }
            }
        }
    }

    // Generate tree within each county
    if (remaining > 0) {
        std::vector<int> path(remaining + 2);
        int max_try = std::max(50, 50 * remaining * ((int)std::log(remaining)));
        std::vector<int> next_vertex(V, -1);
        Frontier precinct_frontier(V);
        for (int i = 0; i < V; i++) {
            if (visited[i]) {
                add_precinct_frontier(g, i, visited, ignore, counties, precinct_frontier);
            }
        }
        while (remaining > 0) {
            int add = precinct_frontier.next(visited, V, lower_i);
            // random walk from `add` until we hit the path
            int added =
                walk_until(g, add, path, max_try, visited, ignore, counties, next_vertex);
            // update visited list and constructed tree
            if (added == 0) { // bail
                return 1;
            }
            remaining -= added - 1; // minus 1 because ending vertex already in tree
            for (int i = 0; i < added - 1; i++) {
                visited.at(path[i]) = true;
                add_precinct_frontier(g, path[i], visited, ignore, counties, precinct_frontier);
                // reverse path so that arrows point away from root
                tree.at(path[i + 1]).push_back(path[i]);
            }
        }
    }

    return 0;
}

/*
 * Random walk along `g` from `root` until something in `visited` is hit
 */
// TESTED
namespace {

int walk_until(const Graph &g, int root, std::vector<int> &path, int MAX,
               const std::vector<bool> &visited, const std::vector<bool> &ignore,
               const uvec &counties, std::vector<int> &next_vertex) {
    int curr = root;
    bool hit_tree = false;
    int county = static_cast<int>(counties[root]);
    for (int tries = 0; tries < MAX; tries++) {
        int proposal = rnbor(g, curr);
        if (ignore[proposal] || static_cast<int>(counties[proposal]) != county) {
            continue;
        }

        // Recording the latest exit from each vertex erases loops implicitly.
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
void walk_until_cty(Multigraph &mg, int root, std::vector<std::vector<int>> &path,
                    const std::vector<bool> &visited, const std::vector<bool> &ignore,
                    std::vector<int> &next_vertex, std::vector<int> &next_edge) {
    path.clear();

    int curr = root;
    bool hit_tree = false;
    int max_tries = visited.size() * 500;
    for (int tries = 0; tries < max_tries; tries++) {
        if (mg[curr].empty()) {
            break;
        }

        int edge = r_int((int)mg[curr].size());
        int proposal = mg[curr][edge][0];
        if (ignore[mg[curr][edge][2]] || ignore[mg[curr][edge][1]]) {
            continue;
        }

        // Recording the latest exit from each county erases loops implicitly.
        next_vertex[curr] = proposal;
        next_edge[curr] = edge;
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
        curr = next_vertex[curr];
    }
}

} // namespace
