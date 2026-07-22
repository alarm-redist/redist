#include "tree_op.h"

#include "random.h"

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
int tree_pop(Tree &ust, int vtx, const std::vector<unsigned int> &pop, std::vector<int> &pop_below,
             std::vector<int> &parent) {
    int pop_at = pop[vtx];
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

void get_tree_pops_below(const Tree &ust, // const FlatGraph &ust, 
    const int root, TreePopStack &stack,
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
            // for (int child : ust.neighbors(vtx)) {
            for (int child : ust[vtx]) {
                total_pop += pop_below[child];
            }
            pop_below[vtx] = total_pop;
        } else {
            // else push again and push the children
            stack.push({vtx, 0, true});
            // for (const auto &child : ust.neighbors(vtx)) {
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


// updates both the vertex labels and the forest adjacency from a directed tree
void assign_region_id_and_forest_from_tree(const Tree &ust,  // const FlatGraph &ust, 
    PlanVector &region_ids,
                                           EdgeBitset &forest_edges,
                                           int root,
                                           const int new_region_id,
                                           const GraphEdgeIndex &graph_edge_index,
                                           CircularQueue<std::pair<int, int>> &vertex_queue) {
    int const V = graph_edge_index.V;
    // optional bounds checking 
    if constexpr (perf_config::bounds_checking) {
        if (root < 0 || root >= V) {
            std::ostringstream oss;
            oss << "assign_region_id_and_forest_from_tree got invalid root.\n";
            oss << "root=" << root << "\n";
            oss << "V=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (static_cast<int>(ust.size()) != V) {
            std::ostringstream oss;
            oss << "assign_region_id_and_forest_from_tree got tree with wrong size.\n";
            oss << "tree.size()=" << ust.size() << "\n";
            oss << "V=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (new_region_id < std::numeric_limits<int8_t>::min() ||
            new_region_id > std::numeric_limits<int8_t>::max()) {
            std::ostringstream oss;
            oss << "assign_region_id_and_forest_from_tree got invalid region id.\n";
            oss << "new_region_id=" << new_region_id << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    // start from the root
    vertex_queue.clear();

    // Handle root separately. The root has no parent edge to preserve or add.
    region_ids[root] = new_region_id;

    // Clear all stale packed forest edges incident to the root.
    for (auto const incident_edge : graph_edge_index.incident_edges[root]) {
        int const u = static_cast<int>(incident_edge.neighbor);

        if constexpr (perf_config::bounds_checking) {
            if (u < 0 || u >= V) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree saw invalid root neighbor.\n";
                oss << "root=" << root << "\n";
                oss << "u=" << u << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        forest_edges.clear_edge_id(incident_edge.edge_id);
    }


    // Push root's children. Every queued vertex now has a real parent.
    // for (int const child : ust.neighbors(root)) {
    for (int const child : ust[root]) {
        if constexpr (perf_config::bounds_checking) {
            if (child < 0 || child >= V) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree saw invalid root child.\n";
                oss << "root=" << root << "\n";
                oss << "child=" << child << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }

            if (child == root) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree saw root self-loop child.\n";
                oss << "root=" << root << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        vertex_queue.push({child, root});
    }

    while (!vertex_queue.empty()) {
        auto const [v, parent] = vertex_queue.pop();

        if constexpr (perf_config::bounds_checking) {
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree popped invalid vertex.\n";
                oss << "v=" << v << "\n";
                oss << "parent=" << parent << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }

            if (parent < -1 || parent >= V) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree got invalid parent.\n";
                oss << "v=" << v << "\n";
                oss << "parent=" << parent << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        region_ids[v] = new_region_id;


        // optional for sanity checking. Should be optimized out by compiler when not checking
        bool found_parent_edge = false;

        // Clear stale packed forest edges incident to v.
        //
        // Do NOT clear the edge to parent, because that is the new tree edge
        // connecting v to the already-processed part of this tree.
        for (auto const incident_edge : graph_edge_index.incident_edges[v]) {
            int const u = static_cast<int>(incident_edge.neighbor);

            if constexpr (perf_config::bounds_checking) {
                if (u < 0 || u >= V) {
                    std::ostringstream oss;
                    oss << "assign_region_id_and_forest_from_tree saw invalid neighbor.\n";
                    oss << "v=" << v << "\n";
                    oss << "u=" << u << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            if (u == parent) {
                // Add the new tree edge from parent to v.
                if constexpr (perf_config::supposedly_safe_input_checks){
                    // this branch adds an extra check that we find a parent edge period. 
                    found_parent_edge = true;
                }else{ 
                    // else we're not checking if a parent edge was found so just set now 
                    forest_edges.set_edge_id(incident_edge.edge_id);
                }
                // continue so we don't remove this edge
                continue;
            }
            // clear the edge
            forest_edges.clear_edge_id(incident_edge.edge_id);
        }
        
        if constexpr (perf_config::supposedly_safe_input_checks) {
            if (!found_parent_edge) {
                std::ostringstream oss;
                oss << "assign_region_id_and_forest_from_tree could not find parent edge.\n";
                oss << "v=" << v << "\n";
                oss << "parent=" << parent << "\n";
                throw std::runtime_error(oss.str());
            }
        }

        // for (int const child : ust.neighbors(v)) {
        for (int const child : ust[v]) {
            if constexpr (perf_config::bounds_checking) {
                if (child < 0 || child >= V) {
                    std::ostringstream oss;
                    oss << "assign_region_id_and_forest_from_tree saw invalid child.\n";
                    oss << "v=" << v << "\n";
                    oss << "child=" << child << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }

                if (child == v) {
                    std::ostringstream oss;
                    oss << "assign_region_id_and_forest_from_tree saw self-loop child.\n";
                    oss << "v=" << v << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            vertex_queue.push({child, v});
        }
    }
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
    if constexpr (perf_config::supposedly_safe_input_checks){
        int const V = static_cast<int>(ust.size());
        if (cut_edge.cut_vertex_parent < 0 || cut_edge.cut_vertex_parent >= V) {
            std::ostringstream oss;
            oss << "erase_tree_edge invalid cut_vertex_parent.\n";
            oss << "cut_vertex_parent=" << cut_edge.cut_vertex_parent << "\n";
            oss << "cut_vertex=" << cut_edge.cut_vertex << "\n";
            oss << "tree_root=" << cut_edge.tree_root << "\n";
            oss << "ust.size()=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (cut_edge.cut_vertex < 0 || cut_edge.cut_vertex >= V) {
            std::ostringstream oss;
            oss << "erase_tree_edge invalid cut_vertex.\n";
            oss << "cut_vertex_parent=" << cut_edge.cut_vertex_parent << "\n";
            oss << "cut_vertex=" << cut_edge.cut_vertex << "\n";
            oss << "tree_root=" << cut_edge.tree_root << "\n";
            oss << "ust.size()=" << V << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    // Get all of the descendents of `cut_vertex_parent`
    std::vector<int> &siblings = ust[cut_edge.cut_vertex_parent];

    // find index of `cut_vertex` among `cut_vertex_parent`'s children
    auto it = std::find(
        siblings.begin(),
        siblings.end(),
        cut_edge.cut_vertex
    );

    if constexpr (perf_config::supposedly_safe_input_checks){
        if (it == siblings.end()) {
            throw std::runtime_error(
                "Actual cut edge not found in valid_edges."
            );
        }
    }

    // Now remove the edge corresponding to `cut_vertex` from the tree
    siblings.erase(it);
}


// checks that every edge in the tree is in the underlying 
// edge list 
void check_directed_tree_edges_are_graph_edges(
    Tree const &tree,
    GraphEdgeIndex const &edge_index,
    std::string_view where
){
    int const V = static_cast<int>(tree.size());

    if (V != edge_index.V) {
        std::ostringstream oss;
        oss << where << ": tree.size() != edge_index.V. "
            << "tree.size()=" << tree.size()
            << ", edge_index.V=" << edge_index.V;
        throw std::runtime_error(oss.str());
    }

    for (int v = 0; v < V; ++v) {
        for (int const u : tree[v]) {
            if (u < 0 || u >= V) {
                std::ostringstream oss;
                oss << where << ": invalid tree neighbor. "
                    << "v=" << v
                    << ", u=" << u
                    << ", V=" << V;
                throw std::runtime_error(oss.str());
            }

            try {
                edge_index.get_edge_id(v, u);
            } catch (std::exception const &e) {
                std::ostringstream oss;
                oss << where << ": tree contains non-graph edge.\n";
                oss << "edge=(" << v << ", " << u << ")\n";
                oss << "underlying error: " << e.what() << "\n";
                throw std::runtime_error(oss.str());
            }
        }
    }
}

void check_is_directed_tree(
    Tree const &ust,
    std::string_view const where,
    int const root,
    int const expected_tree_vertices,
    bool const check_vertex_count,
    std::vector<bool> &visited,
    TreePopStack &stack
) {
    int const V = static_cast<int>(ust.size());

    if (V <= 0) {
        std::ostringstream oss;
        oss << where << ": check_is_directed_tree received empty ust.\n";
        oss << "ust.size()=" << ust.size() << "\n";
        oss << "root=" << root << "\n";
        if(check_vertex_count){
            oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
        }
        throw std::runtime_error(oss.str());
    }

    if (root < 0 || root >= V) {
        std::ostringstream oss;
        oss << where << ": check_is_directed_tree received invalid root.\n";
        oss << "root=" << root << "\n";
        oss << "ust.size()=" << ust.size() << "\n";
        if(check_vertex_count){
            oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
        }
        throw std::runtime_error(oss.str());
    }

    if (check_vertex_count && (expected_tree_vertices <= 0 || expected_tree_vertices > V)) {
        std::ostringstream oss;
        oss << where << ": check_is_directed_tree received invalid expected_tree_vertices.\n";
        oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
        oss << "ust.size()=" << ust.size() << "\n";
        oss << "root=" << root << "\n";
        throw std::runtime_error(oss.str());
    }

    if (static_cast<int>(visited.size()) != V) {
        visited.assign(static_cast<std::size_t>(V), false);
    } else {
        std::fill(visited.begin(), visited.end(), false);
    }

    stack.clear();

    visited[root] = true;
    stack.push({root, -1, false});

    int visited_count = 0;
    int traversed_edge_count = 0;

    while (!stack.empty()) {
        auto const [v, parent, is_revisiting] = stack.pop();

        if (is_revisiting) {
            std::ostringstream oss;
            oss << where << ": check_is_directed_tree unexpectedly popped "
                << "a revisiting stack item.\n";
            oss << "v=" << v << "\n";
            oss << "parent=" << parent << "\n";
            oss << "root=" << root << "\n";
            if(check_vertex_count){
                oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
            }
            throw std::runtime_error(oss.str());
        }

        if (v < 0 || v >= V) {
            std::ostringstream oss;
            oss << where << ": internal traversal reached invalid vertex.\n";
            oss << "v=" << v << "\n";
            oss << "parent=" << parent << "\n";
            oss << "ust.size()=" << ust.size() << "\n";
            oss << "root=" << root << "\n";
            if(check_vertex_count){
                oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
            }
            throw std::runtime_error(oss.str());
        }

        ++visited_count;

        if (check_vertex_count && visited_count > expected_tree_vertices) {
            std::ostringstream oss;
            oss << where << ": directed tree traversal visited too many vertices.\n";
            oss << "visited_count=" << visited_count << "\n";
            oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
            oss << "root=" << root << "\n";
            oss << "current_vertex=" << v << "\n";
            oss << "parent=" << parent << "\n";
            throw std::runtime_error(oss.str());
        }

        for (int const child : ust[v]) {
            if (child < 0 || child >= V) {
                std::ostringstream oss;
                oss << where << ": directed tree contains invalid child vertex.\n";
                oss << "parent=" << v << "\n";
                oss << "child=" << child << "\n";
                oss << "ust.size()=" << ust.size() << "\n";
                oss << "root=" << root << "\n";
                if(check_vertex_count){
                    oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
                }
                throw std::runtime_error(oss.str());
            }

            if (child == v) {
                std::ostringstream oss;
                oss << where << ": directed tree contains self-loop.\n";
                oss << "vertex=" << v << "\n";
                oss << "root=" << root << "\n";
                if(check_vertex_count){
                    oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
                }
                throw std::runtime_error(oss.str());
            }

            ++traversed_edge_count;

            if (visited[child]) {
                std::ostringstream oss;
                oss << where << ": directed tree visits a vertex more than once.\n";
                oss << "duplicate_child=" << child << "\n";
                oss << "encountered_from_parent=" << v << "\n";
                oss << "current_stack_parent=" << parent << "\n";
                oss << "root=" << root << "\n";
                if(check_vertex_count){
                    oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
                }
                oss << "visited_count_so_far=" << visited_count << "\n";

                oss << "Children of duplicate_child: ";
                for (int const grandchild : ust[child]) {
                    oss << grandchild << " ";
                }
                oss << "\n";

                throw std::runtime_error(oss.str());
            }

            visited[child] = true;
            stack.push({child, v, false});
        }
    }

    if (check_vertex_count && visited_count != expected_tree_vertices) {
        std::ostringstream oss;
        oss << where << ": directed tree did not visit expected number of vertices.\n";
        oss << "visited_count=" << visited_count << "\n";
        oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
        oss << "root=" << root << "\n";
        oss << "ust.size()=" << ust.size() << "\n";

        oss << "Visited vertices: ";
        for (int v = 0; v < V; ++v) {
            if (visited[v]) {
                oss << v << " ";
            }
        }
        oss << "\n";

        oss << "Unvisited vertices: ";
        for (int v = 0; v < V; ++v) {
            if (!visited[v]) {
                oss << v << " ";
            }
        }
        oss << "\n";

        throw std::runtime_error(oss.str());
    }

    if (check_vertex_count && traversed_edge_count != expected_tree_vertices - 1) {
        std::ostringstream oss;
        oss << where << ": directed tree has wrong number of traversed edges.\n";
        oss << "traversed_edge_count=" << traversed_edge_count << "\n";
        oss << "expected_edges=" << (expected_tree_vertices - 1) << "\n";
        oss << "visited_count=" << visited_count << "\n";
        oss << "root=" << root << "\n";
        oss << "expected_tree_vertices=" << expected_tree_vertices << "\n";
        throw std::runtime_error(oss.str());
    }
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

                if (region_ids[u] != region_id && test_edge(v, u, edge_index)) {
                    std::ostringstream oss;

                    oss << "EdgeBitset::clear_region_tree found cross-region packed edge.\n";
                    oss << "pair=(" << v << ", " << u << ")\n";
                    oss << "target region=" << static_cast<int>(region_id) << "\n";
                    oss << "region_ids[" << v << "]=" << static_cast<int>(region_ids[v]) << "\n";
                    oss << "region_ids[" << u << "]=" << static_cast<int>(region_ids[u]) << "\n";
                    oss << "Packed forest materialized tree:\n";
                    oss << tree_to_string(get_graph_tree(edge_index));

                    throw std::runtime_error(oss.str());
                }

                // Avoid clearing the same undirected edge twice.
                if (v < u && region_ids[u] == region_id) {
                    // clears the edge 
                    clear_edge_id(incident_edge.edge_id);
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

                if (region_ids[u] != region_id1 &&
                    region_ids[u] != region_id2 &&
                    test_edge(v, u, edge_index)) {

                    std::ostringstream oss;

                    oss << "EdgeBitset::clear_region_trees found cross-region packed edge.\n";
                    oss << "pair=(" << v << ", " << u << ")\n";
                    oss << "target regions=("
                        << static_cast<int>(region_id1) << ", "
                        << static_cast<int>(region_id2) << ")\n";
                    oss << "region_ids[" << v << "]=" << static_cast<int>(region_ids[v]) << "\n";
                    oss << "region_ids[" << u << "]=" << static_cast<int>(region_ids[u]) << "\n";
                    oss << "Packed forest materialized tree:\n";
                    oss << tree_to_string(get_graph_tree(edge_index));

                    throw std::runtime_error(oss.str());
                }

                // Avoid clearing the same undirected edge twice.
                if (v < u && (region_ids[u] == region_id1 || region_ids[u] == region_id2)) {
                    // clears the edge 
                    clear_edge_id(incident_edge.edge_id);
                }


            }
        }
}


void EdgeBitset::fill_flat_tree(
        GraphEdgeIndex const &edge_index,
        FlatGraph &ust
) const {
    if constexpr(perf_config::supposedly_safe_input_checks){
        if (static_cast<int>(ust.size()) != edge_index.V) {
            std::ostringstream oss;

            oss << "EdgeBitset::fill_vector_tree received wrong-sized Tree. "
                << "ust.size()=" << ust.size()
                << ", edge_index.V=" << edge_index.V;

            throw std::runtime_error(oss.str());
        }
    }

    for (int v = 0; v < edge_index.V; ++v) {
        ust.clear_vertex(v);

        for (auto const &incident_edge : edge_index.incident_edges[v]) {
            int const u = static_cast<int>(incident_edge.neighbor);

            if constexpr(perf_config::bounds_checking){
                if (u < 0 || u >= edge_index.V) {
                    std::ostringstream oss;
                    oss << "EdgeBitset::fill_vector_tree saw invalid incident neighbor. "
                        << "v=" << v
                        << ", neighbor=" << u
                        << ", V=" << edge_index.V;

                    throw std::runtime_error(oss.str());
                }
            }

            if (test_edge_id(incident_edge.edge_id)) {
                ust.add_directed_edge(v, u);
            }
        }
    }
}

void EdgeBitset::fill_vector_tree_regions(
    GraphEdgeIndex const &edge_index,
    PlanVector const &region_ids,
    int const region1_id, int const region2_id,
    Tree &ust
) const {
    if constexpr(perf_config::supposedly_safe_input_checks){
        if (static_cast<int>(ust.size()) != edge_index.V) {
            std::ostringstream oss;

            oss << "EdgeBitset::fill_vector_tree received wrong-sized Tree. "
                << "ust.size()=" << ust.size()
                << ", edge_index.V=" << edge_index.V;

            throw std::runtime_error(oss.str());
        }
    }

    for (int v = 0; v < edge_index.V; ++v) {
        auto v_region = region_ids[v];
        // ignore if not the region we want 
        if (v_region != region1_id && v_region != region2_id) continue;
        // clear the vertex
        ust[v].clear();

        for (auto const &incident_edge : edge_index.incident_edges[v]) {
            int const u = static_cast<int>(incident_edge.neighbor);

            if constexpr(perf_config::bounds_checking){
                if (u < 0 || u >= edge_index.V) {
                    std::ostringstream oss;
                    oss << "EdgeBitset::fill_vector_tree saw invalid incident neighbor. "
                        << "v=" << v
                        << ", neighbor=" << u
                        << ", V=" << edge_index.V;

                    throw std::runtime_error(oss.str());
                }
            }

            if (test_edge_id(incident_edge.edge_id)) {
                ust[v].push_back(u);
            }
        }
    }
}

// starting from the root it clears the edges of 
// every vertex in the component 
void EdgeBitset::fill_vector_tree_component_from_root(
    GraphEdgeIndex const &edge_index,
    int const root,
    FlatGraph &ust,
    TreePopStack &stack
) const {
    int const V = edge_index.V;

    if constexpr (perf_config::bounds_checking) {
        if (root < 0 || root >= V) {
            std::ostringstream oss;
            oss << "EdgeBitset::fill_vector_tree_component_from_packed got invalid root.\n";
            oss << "root=" << root << "\n";
            oss << "V=" << V << "\n";
            throw std::runtime_error(oss.str());
        }

        if (static_cast<int>(ust.size()) != V) {
            std::ostringstream oss;
            oss << "EdgeBitset::fill_vector_tree_component_from_packed got wrong-sized Tree.\n";
            oss << "ust.size()=" << ust.size() << "\n";
            oss << "V=" << V << "\n";
            throw std::runtime_error(oss.str());
        }
    }

    stack.clear();
    stack.push({root, -1, false});

    while (!stack.empty()) {
        auto const [v, parent, unused] = stack.pop();

        if constexpr (perf_config::bounds_checking) {
            if (v < 0 || v >= V) {
                std::ostringstream oss;
                oss << "EdgeBitset::fill_vector_tree_component_from_packed popped invalid vertex.\n";
                oss << "v=" << v << "\n";
                oss << "parent=" << parent << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }

            if (parent < -1 || parent >= V) {
                std::ostringstream oss;
                oss << "EdgeBitset::fill_vector_tree_component_from_packed got invalid parent.\n";
                oss << "v=" << v << "\n";
                oss << "parent=" << parent << "\n";
                oss << "V=" << V << "\n";
                throw std::runtime_error(oss.str());
            }
        }
        // clear this vertex 
        ust.clear_vertex(v);
        
        // Now visit every neighbor to add to the tree
        for (auto const &incident_edge : edge_index.incident_edges[v]) {
            int const u = static_cast<int>(incident_edge.neighbor);

            if constexpr (perf_config::bounds_checking) {
                if (u < 0 || u >= V) {
                    std::ostringstream oss;
                    oss << "EdgeBitset::fill_vector_tree_component_from_packed saw invalid neighbor.\n";
                    oss << "v=" << v << "\n";
                    oss << "u=" << u << "\n";
                    oss << "V=" << V << "\n";
                    throw std::runtime_error(oss.str());
                }
            }

            if (!test_edge_id(incident_edge.edge_id)) {
                continue;
            }

            // Fill the vector adjacency for v.
            // This includes the parent edge, so ust is undirected.
            ust.add_directed_edge(v, u);

            // But do not walk back to the parent.
            if (u == parent) {
                continue;
            }

            stack.push({u, v, false});
        }
    }
}

std::string tree_to_string(Tree const &ust) {
    std::ostringstream oss;

    oss << "Printing Tree:\n";

    for (std::size_t i = 0; i < ust.size(); ++i) {
        oss << i << ": (";

        for (int const node : ust[i]) {
            oss << node << " ";
        }

        oss << ")\n";
    }

    return oss.str();
}