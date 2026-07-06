#include "tree_op.h"

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
 * Generate a random neighbor to a vertex, except for the `last` vertex.
 */
// TESTED
int rnbor(const Graph &g, int vtx, RNGState &rng_state) {
    int n_nbors = g[vtx].size();
    return g[vtx][rng_state.r_int(n_nbors)];
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

void print_tree(Tree const &ust) {
    Rprintf("Printing Tree:\n");
    for (int i = 0; i < ust.size(); i++) {
        Rprintf("%d: (", i);
        for (auto &node : ust.at(i)) {
            Rprintf("%d ", node);
        }
        Rprintf(")\n");
    }
}


/*
 * Count population below each node in tree and get each node's parent
 */
// TESTED
// [[Rcpp::export]]
int tree_pop(Tree &ust, int vtx, const arma::uvec &pop, std::vector<int> &pop_below,
             std::vector<int> &parent) {
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

void get_tree_pops_below(const Tree &ust, const int root, TreePopStack &stack,
                         const arma::uvec &pop, std::vector<int> &pop_below) {
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

/*
 * Just Count population below each node in tree
 */
// TESTED
int get_tree_pops_below(const Tree &ust, const int vtx, const arma::uvec &pop,
                        std::vector<int> &pop_below) {
    int pop_at = pop[vtx];
    for (auto const nbor : ust[vtx]) {
        pop_at += get_tree_pops_below(ust, nbor, pop, pop_below);
    }

    pop_below[vtx] = pop_at;
    return pop_at;
}

void assign_region_id_from_tree(Tree const &ust, PlanVector &region_ids, int const root,
                                const int new_region_id,
                                CircularQueue<std::pair<int, int>> &vertex_queue) {
    // clear the queue
    vertex_queue.clear();

    // update root and add its children to queue
    region_ids[root] = new_region_id;
    for (auto const &child_vertex : ust[root]) {
        vertex_queue.push({child_vertex, 0});
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, dont_care] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = new_region_id;
        // add children
        for (auto const &child_vertex : ust[vertex]) {
            vertex_queue.push({child_vertex, 0});
        }
    }

    return;
}

// updates both the vertex labels and the forest adjacency from a directed tree
void assign_region_id_and_forest_from_tree(const Tree &ust, PlanVector &region_ids,
                                           Tree &forest_graph, EdgeBitset &forest_edges,
                                           int root,
                                           const int new_region_id,
                                           const GraphEdgeIndex &graph_edge_index,
                                           CircularQueue<std::pair<int, int>> &vertex_queue) {

    // clear the queue of vertex, parent
    vertex_queue.clear();

    // update root region id
    region_ids[root] = new_region_id;
    // and its forest vertices
    int n_desc = ust[root].size();
    // clear this vertices neighbors in the graph and reserve size for children
    forest_graph[root].clear();
    forest_graph[root].reserve(n_desc);

    // add roots children to queue
    for (auto const &child_vertex : ust[root]) {
        vertex_queue.push({child_vertex, root});
        forest_graph[root].push_back(child_vertex);
        // Now add this edge to the packed forest 
        forest_edges.set_edge(root, child_vertex, graph_edge_index);
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, parent_vertex] = vertex_queue.pop();

        // update region ids
        region_ids[vertex] = new_region_id;
        // clear this vertices neighbors in the graph and reserve size for children and parent
        int n_desc = ust[vertex].size();
        forest_graph[vertex].clear();
        forest_graph[vertex].reserve(n_desc + 1);
        // add the edge from vertex to parent
        forest_graph[vertex].push_back(parent_vertex);
        // set the packed edge 
        forest_edges.set_edge(vertex, parent_vertex, graph_edge_index);

        for (auto const &child_vertex : ust[vertex]) {
            // add children to queue
            vertex_queue.push({child_vertex, vertex});
            // add this edge from vertex to its children
            forest_graph[vertex].push_back(child_vertex);
            // set the packed edge
            forest_edges.set_edge(vertex, child_vertex, graph_edge_index);
        }
    }

    return;
}


// updates both the vertex labels and the forest adjacency
// NOTE: This assumes that the region has already been reset in the packed forest
void assign_region_id_and_forest_from_tree_NEW(const Tree &ust, PlanVector &region_ids,
                                           EdgeBitset &forest_edges, int root,
                                           const int new_region_id,
                                           const GraphEdgeIndex &graph_edge_index,
                                           CircularQueue<std::pair<int, int>> &vertex_queue) {

    // clear the queue of vertex, parent
    vertex_queue.clear();

    // update root region id
    region_ids[root] = new_region_id;


    // add roots children to queue
    for (auto const &child_vertex : ust[root]) {
        vertex_queue.push({child_vertex, root});
        // Now add this edge to the packed forest 
        forest_edges.set_edge(root, child_vertex, graph_edge_index);
    }


    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, parent_vertex] = vertex_queue.pop();

        // update region ids
        region_ids[vertex] = new_region_id;
        // add this edge to the packed forest
        forest_edges.set_edge(vertex, parent_vertex, graph_edge_index);

        // TODO check if set edge is neccesary here. I don't think so but leaving to be safe 
        for (auto const &child_vertex : ust[vertex]) {
            // add children to queue
            vertex_queue.push({child_vertex, vertex});
            // add this edge from vertex to its children
            forest_edges.set_edge(vertex, child_vertex, graph_edge_index);
        }
    }

    return;
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
void erase_tree_edge(Tree &ust, EdgeCut cut_edge) {
    // Get all of the descendents of `cut_vertex_parent`
    std::vector<int> *siblings = &ust[cut_edge.cut_vertex_parent];
    int length = siblings->size();
    int j;
    // find index of `cut_vertex` among `cut_vertex_parent`'s children
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_edge.cut_vertex)
            break;
    }

    // Now remove the edge corresponding to `cut_vertex` from the tree
    siblings->erase(siblings->begin() + j);
}




void EdgeBitset::clear_region_tree(
        PlanVector const &region_ids,
        RegionID const region_id,
        GraphEdgeIndex const &edge_index
){
        int const V = static_cast<int>(region_ids.size());

        // iterate through the graph 
        for (int v = 0; v < V; ++v) {
            // skip if not in that region
            if (region_ids[v] != region_id) {
                continue;
            }

            // Now check each of the edges 
            for (auto const &incident_edge : edge_index.incident_edges[v]) {
                int const u = static_cast<int>(incident_edge.neighbor);

                // Avoid clearing the same undirected edge twice.
                if (v < u && region_ids[u] == region_id) {
                    // clears the edge 
                    clear_edge_id(incident_edge.edge_id);
                }

                if(region_ids[u] != region_id && test_edge(v, u, edge_index)){
                    REprintf("Somehow pair (%d, %d) has an edge despite %d not being in region %d (its in region %d)!\n",
                    u, v, u, static_cast<int>(region_id) ,static_cast<int>(region_ids[u]));
                    print_full_tree(edge_index);
                }
            }
        }
}


void EdgeBitset::clear_region_trees(
        PlanVector const &region_ids,
        RegionID const region_id1, RegionID const region_id2,
        GraphEdgeIndex const &edge_index
){
        int const V = static_cast<int>(region_ids.size());

        // iterate through the graph 
        for (int v = 0; v < V; ++v) {
            // skip if not in that region
            if (region_ids[v] != region_id1 && region_ids[v] != region_id2) {
                continue;
            }

            // Now check each of the edges 
            for (auto const &incident_edge : edge_index.incident_edges[v]) {
                int const u = static_cast<int>(incident_edge.neighbor);

                // Avoid clearing the same undirected edge twice.
                if (v < u && (region_ids[u] == region_id1 || region_ids[u] == region_id2)) {
                    // clears the edge 
                    clear_edge_id(incident_edge.edge_id);
                }

                if( (region_ids[u] == region_id1 || region_ids[u] == region_id2) && test_edge(v, u, edge_index)){
                    REprintf("Somehow pair (%d, %d) has an edge despite %d not being in region %d  or %d (its in region %d)!\n",
                    u, v, u, static_cast<int>(region_id1), static_cast<int>(region_id2),
                    static_cast<int>(region_ids[u]));
                    print_full_tree(edge_index);
                }
            }
        }
}


void EdgeBitset::fill_vector_tree(
        GraphEdgeIndex const &edge_index,
        Tree &ust
) const{
    if constexpr (perf_config::unnecessary_input_checks){
        if (static_cast<int>(ust.size()) != edge_index.V) {
            std::ostringstream oss;
            Rcpp::Rcerr << "EdgeBitset::fill_vector_tree received wrong-sized Tree. "
                << "ust.size()=" << ust.size()
                << ", edge_index.V=" << edge_index.V;
            
            throw Rcpp::exception("fill_vector_tree");
            throw std::runtime_error(oss.str());
        }
    }


    // iterate through the graph
    for (int v = 0; v < edge_index.V; ++v) {
        // clear this vertices edges in the tree 
        ust[v].clear();
        // Check each of v's edges 
        for (auto const &incident_edge : edge_index.incident_edges[v]) {
            // if its adjacent in the packed forest then add it. 
            if (test_edge_id(incident_edge.edge_id)) {
                ust[v].push_back(static_cast<int>(incident_edge.neighbor));
            }
        }
    }
}

