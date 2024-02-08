#include "tree_op.h"

/*
 * Generate a random vertex (integer) among unvisited vertices
 * `lower` is a lower bound (inclusive) on the index of the first unvisited element
 */
// TESTED
int rvtx(const std::vector<bool> &visited, int size, int remaining, int &lower) {
    int idx = r_int(remaining);
    int accuml = 0;
    bool seen_one = false;
    for (int i = lower; i < size - 1; i++) {
        accuml += 1 - visited[i];
        if (!seen_one && !visited[i]) {
            seen_one = true;
            lower = i;
        }
        if (accuml - 1 == idx) return i;
    }
    return size - 1;
}

/*
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
int rnbor(const Graph &g, int vtx) {
    int n_nbors = g[vtx].size();
    return g[vtx][r_int(n_nbors)];
}

/*
 * Make a county graph from a precinct graph and list of counties
 * County graph is list of list of 3: <cty of nbor, index of vtx, index of nbor>
 */
// TESTED
Multigraph county_graph(const Graph &g, const uvec &counties) {
    int n_county = max(counties);
    Multigraph cg = init_multigraph(n_county);

    int V = g.size();
    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        int length = nbors.size();
        int county = counties.at(i) - 1;
        for (int j = 0; j < length; j++) {
            int nbor_cty = counties.at(nbors[j]) - 1;
            if (county == nbor_cty) continue;
            std::vector<int> el = {nbor_cty, i, nbors[j]};
            cg.at(county).push_back(el);
        }
    }

    return cg;
}


/*
 * Make the district adjacency graph for `plan` from the overall precinct graph `g`
 */
// TESTED
Graph district_graph(const Graph &g, const uvec &plan, int nd, bool zero) {
    int V = g.size();
    std::vector<std::vector<bool>> gr_bool;
    for (int i = 0; i < nd; i++) {
        std::vector<bool> tmp(nd, false);
        gr_bool.push_back(tmp);
    }

    for (int i = 0; i < V; i++) {
        std::vector<int> nbors = g[i];
        int dist_i = plan[i] - 1 + zero;
        for (int nbor : nbors) {
            int dist_j = plan[nbor] - 1 + zero;
            if (dist_j != dist_i) {
                gr_bool[dist_i][dist_j] = true;
            }
        }
    }

    Graph out;
    for (int i = 0; i < nd; i++) {
        std::vector<int> tmp;
        for (int j = 0; j < nd; j++) {
            if (gr_bool[i][j]) {
                tmp.push_back(j);
            }
        }
        out.push_back(tmp);
    }

    return out;
}


/*
 * Initialize empty multigraph structure on graph with `V` vertices
 */
// TESTED
Multigraph init_multigraph(int V) {
    Multigraph g;
    for (int i = 0; i < V; i++) {
        std::vector<std::vector<int>> el;
        g.push_back(el);
    }
    return g;
}

/*
 * Initialize empty tree structure on graph with `V` vertices
 */
// TESTED
Tree init_tree(int V) {
    Tree tree;
    for (int i = 0; i < V; i++) {
        tree.push_back(std::vector<int>());
    }
    return tree;
}

/*
 * Initialize empty tree structure on graph with `V` vertices
 */
// TESTED
void clear_tree(Tree &tree) {
    for (auto &nodes : tree) {
        nodes.clear();
    }
}

/*
 * Convert R adjacency list to Graph object (vector of vectors of ints).
 */
Graph list_to_graph(const List &l) {
    int V = l.size();
    Graph g;
    for (int i = 0; i < V; i++) {
        g.push_back(as<std::vector<int>>((IntegerVector) l[i]));
    }
    return g;
}


/*
 * Count population below each node in tree
 */
// TESTED
// [[Rcpp::export]]
int tree_pop(Tree &ust, int vtx, const arma::uvec &pop,
             std::vector<int> &pop_below, std::vector<int> &parent) {
    int pop_at = pop(vtx);
    const std::vector<int> *nbors = &ust[vtx];
    int length = nbors->size();
    for (int j = 0; j < length; j++) {
        int nbor = (*nbors)[j];
        parent.at(nbor) = vtx;
        pop_at += tree_pop(ust, nbor, pop, pop_below, parent);
    }

    pop_below.at(vtx) = pop_at;
    return pop_at;
}

/*
 * Assign `district` to all descendants of `root` in `ust`
 */
// TESTED
void assign_district(const Tree &ust, subview_col<uword> &districts,
                     int root, int district) {
    districts(root) = district;
    int n_desc = ust.at(root).size();
    for (int i = 0; i < n_desc; i++) {
        assign_district(ust, districts, ust.at(root).at(i), district);
    }
}

/*
 * Find the root of a subtree.
 */
// TESTED
int find_subroot(const Tree &ust, const std::vector<bool> &ignore) {
    int V = ust.size();
    std::vector<bool> visited(V, false);
    for (int j = 0; j < V; j++) {
        const std::vector<int>* nbors = &ust[j];
        for (int k = 0; k < nbors->size(); k++) {
            visited[(*nbors)[k]] = true;
        }
    }
    int root;
    for (root = 0; root < V; root++) {
        if (!visited[root] && !ignore.at(root)) break;
    }
    return root;
}
