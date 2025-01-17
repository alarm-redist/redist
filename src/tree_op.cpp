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
 * Count population below each node in tree and get each node's parent
 */
// TESTED
// [[Rcpp::export]]
int tree_pop(Tree &ust, int vtx, const arma::uvec &pop,
             std::vector<int> &pop_below, std::vector<int> &parent) {
    // REprintf("At vertex %d with pop %d\n", vtx, pop(vtx));
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




void assign_region_id_from_tree(const Tree &ust, 
                    arma::subview_col<arma::uword> &region_ids,
                    int root,
                    const int new_region_id) {
    region_ids(root) = new_region_id;
    int n_desc = ust.at(root).size();
    for (int i = 0; i < n_desc; i++) {
        assign_region_id_from_tree(ust, region_ids, ust.at(root).at(i), new_region_id);
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

/*  
 *  Erases an edge from a tree
 * 
 *  Erases the directed edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`)
 *  from the tree `ust`. The directed edge here means we have `child_vertex` being one of 
 *  the values in `ust[parent_vertex]`.
 * 
 * 
 *  @param ust A directed spanning tree passed by reference
 *  @param cut_edge An `EdgeCut` object representing the edge cut
 * 
 *  @details Modifications
 *     - The edge (`cut_edge.cut_vertex_parent`, `cut_edge.cut_vertex`) 
 *     is removed from `ust`
 * 
 * 
 */ 
void erase_tree_edge(Tree &ust, EdgeCut cut_edge){
    // Get all of the descendents of `cut_vertex_parent` 
    std::vector<int> *siblings = &ust[cut_edge.cut_vertex_parent];
    int length = siblings->size();
    int j;
    // find index of `cut_vertex` among `cut_vertex_parent`'s children
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_edge.cut_vertex) break;
    }

    // Now remove the edge corresponding to `cut_vertex` from the tree
    siblings->erase(siblings->begin()+j);
}