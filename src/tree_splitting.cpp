/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Functions for finding an edge in a tree to
remove (splitting the tree)
********************************************************/

#include "tree_splitting.h"
#include "scoring.h"
#include "random.h"
#include "base_plan_type.h"

constexpr bool FINDING_EDGE_CUTS_VERBOSE = false;
constexpr bool FINDING_JOINED_EDGE_CUTS_VERBOSE = false;

namespace{


/*
 * Calls `emit_valid_cut` once for every valid EdgeCut associated with
 * the tree edge (cut_vertex_parent, cut_vertex).
 *
 * This function determines validity only. It does not decide whether the
 * valid cuts should be stored, counted, or weighted. That is all downstream
 * of this so it can be used flexibly.
 *
 * EmitValidCut should be callable as:
 *
 *     emit_valid_split(
                int const cut_region1_size,
                double const region1_pop,
                int const cut_region2_size,
                double const region2_pop
            )
 *
 */
template <typename EmitValidSplit>
inline void for_each_valid_population_split_from_edge(
    int const total_region_size,
    double const below_pop,
    double const above_pop,
    double const lower,
    double const upper,
    std::vector<int> const &smaller_cut_sizes_to_try,
    EmitValidSplit &&emit_valid_split
) {
    for (int const cut_region1_size : smaller_cut_sizes_to_try) {
        // get the size of the other region
        int const cut_region2_size = total_region_size - cut_region1_size;
        // now compute the lower and upper pop bounds for the two new potential regions 
        double const cut_region1_lb = lower * cut_region1_size;
        double const cut_region1_ub = upper * cut_region1_size;
        double const cut_region2_lb = lower * cut_region2_size;
        double const cut_region2_ub = upper * cut_region2_size;

        if constexpr (FINDING_EDGE_CUTS_VERBOSE) {
            REprintf(
                "\tFor (%d, %d): compare %f vs %f vs %f "
                "and %f vs %f vs %f\n",
                cut_region1_size,
                cut_region2_size,
                cut_region1_lb,
                below_pop,
                cut_region1_ub,
                cut_region2_lb,
                above_pop,
                cut_region2_ub
            );
        }

        /*
         * Assign cut_region1_size to the below component and
         * cut_region2_size to the above component.
         */
        bool const cut_below_ok =
            cut_region1_lb <= below_pop && below_pop <= cut_region1_ub &&
            cut_region2_lb <= above_pop && above_pop <= cut_region2_ub;

        // if this cut size is ok then we call `emit_valid_cut` on the associated edgecut 
        if (cut_below_ok) {
            emit_valid_split(
                cut_region1_size,
                below_pop,
                cut_region2_size,
                above_pop
            );
        }

        /*
         * Swapping equal sizes would produce the same EdgeCut, so do
         * not emit it twice.
         */
        if (cut_region1_size == cut_region2_size) {
            continue;
        }

        if constexpr (
            FINDING_EDGE_CUTS_VERBOSE
        ) {
            REprintf(
                "\tFor (%d, %d): compare %f vs %f vs %f "
                "and %f vs %f vs %f\n",
                cut_region2_size,
                cut_region1_size,
                cut_region2_lb,
                below_pop,
                cut_region2_ub,
                cut_region1_lb,
                above_pop,
                cut_region1_ub
            );
        }

        /*
         * Swapped orientation: 
         * Assign cut_region2_size to the below component and
         * cut_region1_size to the above component.
         */
        bool const cut_above_ok =
            cut_region2_lb <= below_pop && below_pop <= cut_region2_ub &&
            cut_region1_lb <= above_pop && above_pop <= cut_region1_ub;

        if (cut_above_ok) {
            emit_valid_split(
                cut_region2_size,
                below_pop,
                cut_region1_size,
                above_pop
            );
        }
    }
}

/*
 * Appends all valid edge cuts for a particular edge and cut region sizes range
 *
 *
 *
 * Given an edge in a spanning tree and a range of region sizes to consider for
 * the two cut regions this appends to the input edge cut vector all the the valid edge cuts
 * (ie an edge and region sizes for the two cuts) that can be made. For strict
 * population bounds this should be only one but if bounds are loose it can
 * be multiple. An empty vector means there are no valid edge cuts.
 *
 *
 *
 * @param existing_cuts A vector of edge cuts 
 * @param root The root vertex of the spanning tree
 * @param cut_vertex The vertex where we are cutting below it
 * @param cut_vertex_parent The parent of `cut_vertex` (if we think of the tree
 * as directed) so we are cutting `(cut_vertex_parent, cut_vertex)`.
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param below_pop The population of the region induced by cutting below
 * `cut_vertex` (ie the region where `cut_vertex` is the root)
 * @param above_pop The population of the region induced by cutting above
 * `cut_vertex` (ie the region where `cut_vertex_parent` is the root)
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 * @param cut_size_loop_start The starting value of the for loop for the range
 * of potential cut region sizes to loop over.
 * @param cut_size_loop_end The final value of the for loop for the range
 * of potential cut region sizes to loop over.
 *
 * @details Appends any valid edge cuts to existing_cuts
 *
 * @return A vector of EdgeCut objects
 *
 */
inline void
get_all_valid_edge_cuts_from_edge(std::vector<EdgeCut> &existing_cuts, int const root, int const cut_vertex,
                                  int const cut_vertex_parent, int const total_region_size,
                                  double const below_pop, double const above_pop,
                                  double const lower, double const upper,
                                  std::vector<int> const &smaller_cut_sizes_to_try) {
    for_each_valid_population_split_from_edge(
        total_region_size,
        below_pop,
        above_pop,
        lower,
        upper,
        smaller_cut_sizes_to_try,
        // the anonymous lambda function just builds an edge cut from the pop info and pushes it back 
        [&](
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
        ) {
            existing_cuts.emplace_back(
                root,
                cut_vertex,
                cut_vertex_parent,
                cut_below_region_size,
                cut_below_pop,
                cut_above_region_size,
                cut_above_pop
            );
        }
    );
    return;
}




// takes an uncut tree and an edge in the tree
// and updates the region ids as if the cut edge was really removed
// NOTE: When Wilson trees switch to undirected this can probably 
// be deleted and replaced with assign_region_ids_from_joined_undirected_tree
void assign_region_ids_from_uncut_tree(Tree const &ust, // FlatGraph const &ust, 
    PlanVector &region_ids, int const root,
                                       int const root_region_id, int const cut_vertex_root,
                                       int const cut_vertex_parent, int const cut_region_id,
                                       CircularQueue<std::pair<int, int>> &vertex_queue) {
    // clear the queue
    vertex_queue.clear();

    // update root and add its children to queue
    region_ids[root] = root_region_id;
    // for (auto const &child_vertex : ust.neighbors(root)) {
    for (auto const &child_vertex : ust[root]) {
        // check if this edge is one to remove in which case we ignore
        // NOTE/TODO When trees become undirected need to check below as well
        // child_vertex == cut_vertex_parent && root == cut_vertex_root
        if (root == cut_vertex_parent && child_vertex == cut_vertex_root)
            continue;

        vertex_queue.push({child_vertex, 0});
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, dont_care] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = root_region_id;
        // add children
        // for (auto const &child_vertex : ust.neighbors(vertex)) {
        for (auto const &child_vertex : ust[vertex]) {
            // check if this edge is one to remove in which case we ignore
            // NOTE/TODO When trees become undirected need to check below as well
            // child_vertex == cut_vertex_parent && vertex == cut_vertex_root
            if (vertex == cut_vertex_parent && child_vertex == cut_vertex_root)
                continue;

            vertex_queue.push({child_vertex, 0});
        }
    }

    // now we update starting at the cut root vertex
    // update root and add its children to queue
    region_ids[cut_vertex_root] = cut_region_id;
    // for (auto const &child_vertex : ust.neighbors(cut_vertex_root)) {
    for (auto const &child_vertex : ust[cut_vertex_root]) {
        // since tree is directed we don't need to check anything
        vertex_queue.push({child_vertex, 0});
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, dont_care] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = cut_region_id;
        // add children
        // for (auto const &child_vertex : ust.neighbors(vertex)) {
        for (auto const &child_vertex : ust[vertex]) {
            // since tree is directed we don't need to check anything
            // TODO/NOTE When trees become undirected we do need to check this 
            vertex_queue.push({child_vertex, 0});
        }
    }

    return;
}

// Given a vector of edge cuts it 
// - removes any edge cuts that don't satisfy hard constraints 
// - Gives the remaining edge cuts an unnormalized weight equal to e^(-score)
arma::vec compute_constraint_edge_cut_weights(
    std::vector<EdgeCut> &valid_edges, ScoringFunction const &scoring_function, Tree const &ust, // FlatGraph const &ust,
    int const num_regions, PlanVector &region_ids, RegionSizes &region_sizes,
    IntPlanAttribute &region_pops, int const split_region_id1, int const split_region_id2,
    CircularQueue<std::pair<int, int>> &vertex_queue) {
    // Allocate enough space for every edge. We shrink this after removing
    // edges that violate a hard constraint.
    arma::vec unnormalized_wgts(valid_edges.size());

    // read_idx identifies the edge currently being evaluated.
    // write_idx identifies where the next valid edge and weight should go.
    std::size_t write_idx = 0;

    for (std::size_t read_idx = 0;
         read_idx < valid_edges.size();
         ++read_idx) {

        EdgeCut const &edge_cut = valid_edges[read_idx];

        // Temporarily update the split-region sizes for this candidate cut.
        region_sizes[split_region_id1] =
            edge_cut.cut_above_region_size;
        region_sizes[split_region_id2] =
            edge_cut.cut_below_region_size;

        // Temporarily update the split-region populations.
        region_pops[split_region_id1] =
            edge_cut.cut_above_pop;
        region_pops[split_region_id2] =
            edge_cut.cut_below_pop;

        // Assign region IDs corresponding to this candidate cut.
        assign_region_ids_from_uncut_tree(
            ust,
            region_ids,
            edge_cut.tree_root,
            split_region_id1,
            edge_cut.cut_vertex,
            edge_cut.cut_vertex_parent,
            split_region_id2,
            vertex_queue
        );

        // first  = whether all applicable hard constraints are satisfied
        // second = total soft-constraint score
        std::pair<bool, double> const score_result =
            scoring_function.compute_full_split_plan_score(
                num_regions,
                region_ids,
                region_sizes,
                region_pops,
                split_region_id1,
                split_region_id2
            );

        // Do not retain this cut if it violates a hard constraint.
        if (!score_result.first) {
            continue;
        }

        // Compact the valid edges toward the beginning of the vector.
        //
        // If no invalid edges have appeared yet, read_idx == write_idx and
        // no move is needed. Once an invalid edge is skipped, later valid
        // edges are moved left to fill the gap.
        if (write_idx != read_idx) {
            valid_edges[write_idx] =
                std::move(valid_edges[read_idx]);
        }

        // The weight index must correspond to the compacted edge index,
        // not the original read index.
        unnormalized_wgts[write_idx] =
            std::exp(-score_result.second);

        ++write_idx;
    }


    // Remove the unused tail containing rejected or moved-from edges.
    valid_edges.resize(write_idx);

    // Keep exactly one weight for every retained edge.
    unnormalized_wgts.resize(write_idx);


    return unnormalized_wgts;
}

// TODO: See if this performs better just appending to the vector instead of retur by value
std::vector<EdgeCut> get_all_valid_edges_in_directed_tree(
    Tree const &a_ust, // FlatGraph const &a_ust, 
    const int root, const std::vector<unsigned int> &pop, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, const int total_region_pop,
    const int total_region_size, const double lower, const double upper, const double target) {

    std::vector<EdgeCut> valid_edges;
    // this is the largest size a region can be
    // If the population above is bigger than this you can terminate the serach
    // since pop above only gets larger as you continue down the tree
    double biggest_upper_bound = upper * max_potential_cut_size;

    // this is the smallest size a region can be
    // If the pop below is below this then you can terminate the search since
    // pop below only gets smaller as you continue along the tree
    double smallest_lower_bound = lower * min_potential_cut_size;

    // Stack for DFS
    // Elements are: vertex, parent, is_revisiting
    // std::stack<std::tuple<int, int, bool>> stack;
    // Clear the stack now
    stack.clear();
    // set the root to false 
    no_valid_edges_vertices[root] = false;

    // Start by adding all the roots children to the stack
    // for (auto const &root_children : a_ust.neighbors(root)) {
    for (auto const &root_children : a_ust[root]) {
        stack.push({root_children, root, false});
    }

    // Loop until the stack is empty
    while (!stack.empty()) {
        // get the top of the stack
        auto popped = stack.pop();

        int const vtx = std::get<0>(popped);
        int const parent = std::get<1>(popped);
        bool const is_revisiting = std::get<2>(popped);

        if (!is_revisiting) { // This is the first time visiting the node
            // reset this vertex to false 
            no_valid_edges_vertices[vtx] = false;

            // Push the vertex back onto the stack as "revisiting"
            stack.push({vtx, parent, true});

            // Push unvisited child vertices onto the stack to get pop below
            // for (const auto &child_vtx : a_ust.neighbors(vtx)) {
            for (const auto &child_vtx : a_ust[vtx]) {
                // else add to the stack
                stack.push({child_vtx, vtx, false});
            }
        } else if (no_valid_edges_vertices[vtx]) {
            // if parent isn't valid then neither is its parent so mark that
            no_valid_edges_vertices[parent] = true;
        } else if (!no_valid_edges_vertices[parent]) {
            // if revisiting it true that means we already visited all the nodes children
            // so we can get pop_below
            // if no valid edges is true we no there's no point in searching up this path
            // anymore

            // All children of this vertex are processed; calculate its population below
            int pop_below_vtx = pop[vtx]; // Start with the vertex's own population
            // Add population below from each child
            // for (const auto &child : a_ust.neighbors(vtx)) {
            for (const auto &child : a_ust[vtx]) {
                pop_below_vtx += pops_below_vertex[child]; // Add population from child vertices
            }
            pops_below_vertex[vtx] = pop_below_vtx;

            // Check if any cut can be made
            // If pop below is too small we need to keep going up
            if (pop_below_vtx < smallest_lower_bound ||
                total_region_pop - pop_below_vtx > biggest_upper_bound) {
                continue;
            } else if (pop_below_vtx > biggest_upper_bound ||
                       total_region_pop - pop_below_vtx < smallest_lower_bound) {
                no_valid_edges_vertices[parent] = true;
                continue;
                // Recall pop below is only increasing for the parent so we can skip this entire
                // lineage if we want
            }

            // See if any valid edge cuts can be made with this edge
            // if there are any the function will automatically add them 
            get_all_valid_edge_cuts_from_edge(
                valid_edges,
                root, vtx, parent, total_region_size, pops_below_vertex[vtx],
                total_region_pop - pops_below_vertex[vtx], lower, upper,
                smaller_cut_sizes_to_try);
        }
    }

    // Return the total population at the root
    return valid_edges;
}

/*
 * Appends all valid edge cuts in the tree to existing_cuts
 *
 *
 * Returns a vector of all the valid edge cuts (ie an edge and regions for
 * the two cuts) where at least one of the regions is between
 * `min_potential_cut_size` and `max_potential_cut_size` inclusive. Returns
 * the edges as EdgeCut objects. An empty vector means there are no
 * valid edges.
 *
 *
 *
 * @param ust A directed edge spanning tree.
 * @param root The root vertex of the spanning tree.
 * @param cut_below_pops The population corresponding to cutting below each vertex.
 * So `cut_below_pops[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param min_potential_cut_size The smallest potential region size at least one of
 * the regions cut must be
 * @param max_potential_cut_size The largest potential region size at least one of
 * the regions cut must be. Setting this to 1 will result in only 1 district splits.
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details Adds valid cuts to existing_cuts
 *
 * @return A vector of EdgeCut objects
 *
 */
void get_all_valid_edges_in_undirected_tree(
    std::vector<EdgeCut> &existing_cuts,
    GraphEdgeIndex const &edge_index,
    EdgeBitset const &forest_edges, 
    const int root, const std::vector<unsigned int> &pop, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, const int total_region_pop,
    const int total_region_size, const double lower, const double upper, const double target) {

    // this is the largest size a region can be
    // If the population above is bigger than this you can terminate the serach
    // since pop above only gets larger as you continue down the tree
    double biggest_upper_bound = upper * max_potential_cut_size;

    // this is the smallest size a region can be
    // If the pop below is below this then you can terminate the search since
    // pop below only gets smaller as you continue along the tree
    double smallest_lower_bound = lower * min_potential_cut_size;

    // Stack for DFS
    // Elements are: vertex, parent, is_revisiting
    stack.clear();
    no_valid_edges_vertices[root] = false;


    // Start by adding all the roots children to the stack
    // This essentially iterates over each vertex which is adjacent
    // to the root in the forest edge and applies the anonymous function
    // so root_child is the neighbor of root
    forest_edges.for_each_neighbor(root, edge_index, [&](int const root_child) {
        stack.push({root_child, root, false});
    });

    // Loop until the stack is empty
    while (!stack.empty()) {
        // get the top of the stack
        auto popped = stack.pop();

        int const vtx = std::get<0>(popped);
        int const parent = std::get<1>(popped);
        bool const is_revisiting = std::get<2>(popped);

        if (!is_revisiting) { // This is the first time visiting the node
            no_valid_edges_vertices[vtx] = false;

            // Push the vertex back onto the stack as "revisiting"
            stack.push({vtx, parent, true});


            // Push unvisited child vertices onto the stack to get pop below
            // Again using the anonymous function
            forest_edges.for_each_neighbor(vtx, edge_index, [&](int const child_vtx) {
                // skip if its the parent
                if (child_vtx == parent) {
                    return;
                }
                stack.push({child_vtx, vtx, false});
            });
        } else if (no_valid_edges_vertices[vtx]) {
            // if parent isn't valid then neither is its parent so mark that
            no_valid_edges_vertices[parent] = true;
        } else if (!no_valid_edges_vertices[parent]) {
            // if revisiting it true that means we already visited all the nodes children
            // so we can get pop_below
            // if no valid edges is true we no there's no point in searching up this path
            // anymore

            // All children of this vertex are processed; calculate its population below
            int pop_below_vtx = pop[vtx]; // Start with the vertex's own population
            

            // Add population below from each child
            // again using anonymous lambdas
            forest_edges.for_each_neighbor(vtx, edge_index, [&](int const child) {
                // ignore the parent 
                if (child == parent) {
                    return;
                }
                // sum up the pop below 
                pop_below_vtx += pops_below_vertex[child];
            });

            pops_below_vertex[vtx] = pop_below_vtx;

            // Check if any cut can be made
            // If pop below is too small we need to keep going up
            if (pop_below_vtx < smallest_lower_bound ||
                total_region_pop - pop_below_vtx > biggest_upper_bound) {
                continue;
            } else if (pop_below_vtx > biggest_upper_bound ||
                       total_region_pop - pop_below_vtx < smallest_lower_bound) {
                no_valid_edges_vertices[parent] = true;
                continue;
                // Recall pop below is only increasing for the parent so we can skip this entire
                // lineage if we want
            }

            // See if any valid edge cuts can be made with this edge
            // if yes,then the function will add them
            get_all_valid_edge_cuts_from_edge(
                existing_cuts,
                root, vtx, parent, total_region_size, pops_below_vertex[vtx],
                total_region_pop - pops_below_vertex[vtx], lower, upper,
                smaller_cut_sizes_to_try);
        }
    }

}

/*
 * Calls the function `emit_valid_cut` on every potentially splittable edge 
 * in the tree made by joining region 1 and region 2's trees together at the 
 * edge (region1_root, region2_root). 
 * 
 * It proceeds by first going through all internal edges in region 1's tree
 * then all internal edges in region 2's tree then finally the edge across
 * (region1_root, region2_root)
 *  
 */
template <typename EmitValidCut>
void for_each_valid_edge_cut_in_joined_packed_tree(
    MapParams const &map_params,
    EdgeBitset const &forest_edges, 
    TreePopStack &stack,
    std::vector<int> &pops_below_vertex,
    std::vector<bool> &no_valid_edges_vertices,
    int const region1_root,
    int const region1_pop,
    int const region2_root,
    int const region2_pop,
    int const min_potential_cut_size,
    int const max_potential_cut_size,
    std::vector<int> const
        &smaller_cut_sizes_to_try,
    int const total_merged_region_size,
    EmitValidCut &&emit_valid_cut
) {
    int const total_merged_region_pop = region1_pop + region2_pop;

    /*
     * Internal edges of region 1.
     *
     * The population of region 2 is treated as lying above the root of
     * region 1 because the joining edge attaches there.
     */
    for_each_valid_edge_cut_in_undirected_packed_tree(
        map_params.graph_edge_index,
        forest_edges,
        region1_root,
        map_params.pop,
        stack,
        pops_below_vertex,
        no_valid_edges_vertices,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_pop,
        total_merged_region_size,
        map_params.lower,
        map_params.upper,
        emit_valid_cut
    );

    /*
     * Internal edges of region 2.
     */
    for_each_valid_edge_cut_in_undirected_packed_tree(
        map_params.graph_edge_index,
        forest_edges,
        region2_root,
        map_params.pop,
        stack,
        pops_below_vertex,
        no_valid_edges_vertices,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_pop,
        total_merged_region_size,
        map_params.lower,
        map_params.upper,
        emit_valid_cut
    );

    /*
     * The newly inserted joining edge itself. Cutting below region2_root
     * recovers region 2, while the above component is region 1.
     */
    for_each_valid_population_split_from_edge(
        total_merged_region_size,
        static_cast<double>(
            region2_pop
        ),
        static_cast<double>(
            region1_pop
        ),
        map_params.lower,
        map_params.upper,
        smaller_cut_sizes_to_try,
        [&](
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
        ) {
            emit_valid_cut(
                region1_root,
                region2_root,
                region1_root,
                cut_below_region_size,
                cut_below_pop,
                cut_above_region_size,
                cut_above_pop
            );
        }
    );
}

// finds all valid edges if you joined the two trees
// with the edge (region1_root, region2_root)
// THIS INCLUDES (region1_root, region2_root) as an edge!!
std::vector<EdgeCut> get_valid_pop_edges_in_joined_packed_tree(
    MapParams const &map_params, EdgeBitset const &forest_edges, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int region1_root, const int region1_pop, const int region2_root,
    const int region2_pop, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, const int total_merged_region_size) {
    int const total_merged_region_pop = region1_pop + region2_pop;

    // Don't need to reset pops below as it starts from the child nodes 
    // so it overwrites old values 
    // Also don't need to reset no_valid_edges_vertices as we 
    // can just set them to false the first time we see the vertex when 
    // traversing the tree 
    std::vector<EdgeCut> valid_edges;

    for_each_valid_edge_cut_in_joined_packed_tree(
        map_params,
        forest_edges,
        stack,
        pops_below_vertex,
        no_valid_edges_vertices,
        region1_root,
        region1_pop,
        region2_root,
        region2_pop,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_size,
        [&](
            int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
        ) {
            valid_edges.emplace_back(
                root,
                cut_vertex,
                cut_vertex_parent,
                cut_below_region_size,
                cut_below_pop,
                cut_above_region_size,
                cut_above_pop
            );
        }
    );

    if constexpr (FINDING_JOINED_EDGE_CUTS_VERBOSE) {
        pops_below_vertex[region1_root] = region1_pop;
        pops_below_vertex[region2_root] = region2_pop;
        REprintf("Pop below region2_root is %d so above is %d so foound %d\n",
                 pops_below_vertex.at(region2_root),
                 total_merged_region_pop - pops_below_vertex.at(region2_root),
                 (int) valid_edges.size());
    }


    return valid_edges;
}

/*
 * Calls emit_valid_cut once per each edge (skips some edge's we know are unsplittable)
 * Starting from the root
 * NOTE: This adds every root vertex so if being used for region tree boundary 
 * computation make sure there is no edge from the root to the other tree
 *
 * 
 */
template <typename EmitValidCut>
void for_each_valid_edge_cut_in_undirected_packed_tree(
    GraphEdgeIndex const &edge_index,
    EdgeBitset const &forest_edges, 
    int const root,
    std::vector<unsigned int> const &pop,
    TreePopStack &stack,
    std::vector<int> &pops_below_vertex,
    std::vector<bool> &no_valid_edges_vertices,
    int const min_potential_cut_size,
    int const max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    int const total_region_pop,
    int const total_region_size,
    double const lower,
    double const upper,
    EmitValidCut &&emit_valid_cut
){

    // this is the largest size a region can be
    // If the population above is bigger than this you can terminate the serach
    // since pop above only gets larger as you continue down the tree
    double const biggest_upper_bound = upper * max_potential_cut_size;

    // this is the smallest size a region can be
    // If the pop below is below this then you can terminate the search since
    // pop below only gets smaller as you continue along the tree
    double const smallest_lower_bound = lower * min_potential_cut_size;

    // Stack for DFS
    // Elements are: vertex, parent, is_revisiting
    stack.clear();
    no_valid_edges_vertices[root] = false;


    // Start by adding all the roots children to the stack
    forest_edges.for_each_neighbor(root, edge_index, [&](int const root_child) {
        stack.push({root_child, root, false});
    });

    // Loop until the stack is empty
    while (!stack.empty()) {
        // get the top of the stack
        auto popped = stack.pop();

        int const vtx = std::get<0>(popped);
        int const parent = std::get<1>(popped);
        bool const is_revisiting = std::get<2>(popped);

        if (!is_revisiting) { // This is the first time visiting the node
            no_valid_edges_vertices[vtx] = false;

            // Push the vertex back onto the stack as "revisiting"
            stack.push({vtx, parent, true});


            // Push unvisited child vertices onto the stack to get pop below
            forest_edges.for_each_neighbor(vtx, edge_index, [&](int const child_vtx) {
                // if its the parent then skip it
                if (child_vtx == parent)
                    return;
                // else add to the stack
                stack.push({child_vtx, vtx, false});
            });
        } else if (no_valid_edges_vertices[vtx]) {
            // if parent isn't valid then neither is its parent so mark that
            no_valid_edges_vertices[parent] = true;
        } else if (!no_valid_edges_vertices[parent]) {
            // if revisiting it true that means we already visited all the nodes children
            // so we can get pop_below
            // if no valid edges is true we no there's no point in searching up this path
            // anymore

            // All children of this vertex are processed; calculate its population below
            int pop_below_vtx = pop[vtx]; // Start with the vertex's own population
            

            // Add population below from each child
            forest_edges.for_each_neighbor(vtx, edge_index, [&](int const child) {
                // ignore the parent 
                if (child == parent) {
                    return;
                }
                pop_below_vtx += pops_below_vertex[child]; // Add population from child vertices
            });
            pops_below_vertex[vtx] = pop_below_vtx;

            int const pop_above_vtx = total_region_pop - pop_below_vtx;

            // Check if any cut can be made
            // If pop below is too small we need to keep going up
            if (pop_below_vtx < smallest_lower_bound ||
                pop_above_vtx > biggest_upper_bound) {
                continue;
            } else if (pop_below_vtx > biggest_upper_bound ||
                       pop_above_vtx < smallest_lower_bound) {
                no_valid_edges_vertices[parent] = true;
                continue;
                // Recall pop below is only increasing for the parent so we can skip this entire
                // lineage if we want
            }

            // See if any valid edge cuts can be made with this edge
            // if yes,then the function will add them
        for_each_valid_population_split_from_edge(
            total_region_size,
            static_cast<double>(pop_below_vtx),
            static_cast<double>(pop_above_vtx),
            lower, upper,
            smaller_cut_sizes_to_try,
            [&](
                int const cut_below_region_size,
                double const cut_below_pop,
                int const cut_above_region_size,
                double const cut_above_pop
            ) {
                emit_valid_cut(
                    root,
                    vtx,
                    parent,
                    cut_below_region_size,
                    cut_below_pop,
                    cut_above_region_size,
                    cut_above_pop
                );
            }
        );
        }
    }
}



/*
 * Computes, for every vertex v in the packed region tree, the total
 * unnormalized weight contributed by all INTERNAL edges of this tree if the
 * other region tree were attached at v.
 *
 * The tree is stored in `forest_edges`; `edge_index` is used to enumerate the
 * packed forest neighbors of each vertex.
 *
 * On return, for every vertex v in this tree:
 *
 *     reroot_weight[v] = F(v)
 *
 * where F(v) is the sum of the unnormalized weights of all valid cuts on
 * internal edges of this tree when the other region tree is attached at v.
 *
 * This does NOT include the weight contribution from the joining edge between
 * this tree and the other region tree. That contribution must be added
 * separately when evaluating a particular joining edge.
 *
 * The algorithm uses two tree traversals:
 *
 * 1. Postorder traversal:
 *      - compute the subtree population below every nonroot vertex v;
 *      - compute W_out(v), the total valid-cut weight for edge
 *        (parent(v), v) when the other tree is attached outside v's subtree;
 *      - compute W_in(v), the total valid-cut weight for the same edge when
 *        the other tree is attached inside v's subtree;
 *      - store
 *
 *            Delta(v) = W_in(v) - W_out(v)
 *
 *        temporarily in reroot_weight[v];
 *      - accumulate
 *
 *            F(root) = sum_{v != root} W_out(v).
 *
 * 2. Preorder traversal:
 *      - set
 *
 *            reroot_weight[root] = F(root);
 *
 *      - for each parent-child edge, apply the rerooting recurrence
 *
 *            F(child) = F(parent) + Delta(child);
 *
 *        overwriting the previously stored Delta(child) with the final
 *        F(child) value;
 *      - call `emit_finished_vertex(v)` immediately after F(v) has been
 *        finalized.
 *
 * The `emit_finished_vertex` callback allows the caller to perform work as
 * each vertex's final reroot weight becomes available. In particular, when
 * computing effective boundary lengths, one region tree can be processed
 * first and its F(v) values left in `reroot_weight`. While processing the
 * second tree, `emit_finished_vertex(v)` can inspect map-graph neighbors of v
 * and immediately evaluate boundary edges using the already-computed F values
 * from both trees. This avoids storing all boundary-edge endpoints or making
 * an additional traversal of the second tree.
 *
 * ComputeCutWeight must be callable as:
 *
 *     double(
 *         int root,
 *         int cut_vertex,
 *         int cut_vertex_parent,
 *         int cut_below_region_size,
 *         double cut_below_pop,
 *         int cut_above_region_size,
 *         double cut_above_pop
 *     )
 *
 * and must return the unnormalized selection weight of that valid cut.
 *
 * EmitFinishedVertex must be callable as:
 *
 *     void(int vertex)
 *
 * and is called exactly once for every vertex in this region tree, after
 * reroot_weight[vertex] has been set to its final F(vertex) value.
 *
 * `forest_edges` must contain this region's tree and must not contain an edge
 * joining it to the other region tree. Because the packed forest is undirected,
 * the traversal explicitly tracks and skips each vertex's parent.
 */
template <
    typename ComputeCutWeight,
    typename EmitFinishedVertex
>
void compute_tree_path_weights_in_undirected_packed_tree(
    GraphEdgeIndex const &edge_index,
    EdgeBitset const &forest_edges,
    int const root,
    std::vector<unsigned int> const &pop,
    TreePopStack &stack,
    std::vector<int> &pops_below_vertex,
    std::vector<double> &reroot_weight,
    int const min_potential_cut_size,
    int const max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    int const this_region_tree_pop,
    int const other_region_tree_pop,
    int const total_region_size,
    double const lower,
    double const upper,
    ComputeCutWeight &&compute_cut_weight,
    EmitFinishedVertex &&emit_finished_vertex
) {
    /*
     * Broadest possible population bounds for either side of a valid cut.
     * These are used only as a cheap test before entering the cut-size loop.
     */
    double const biggest_upper_bound =
        upper * max_potential_cut_size;

    double const smallest_lower_bound =
        lower * min_potential_cut_size;

    /*
     * ============================================================
     * PASS 1: POSTORDER
     * ============================================================
     *
     * For every non-root vertex v:
     *
     *   1. Compute p_v, the population in the subtree below v.
     *   2. Compute W_out(v), where the other tree is attached
     *      outside the subtree below v.
     *   3. Compute W_in(v), where the other tree is attached
     *      inside the subtree below v.
     *   4. Store
     *
     *          reroot_weight[v] = W_in(v) - W_out(v)
     *
     *      which is Delta(v).
     *
     * At the same time accumulate
     *
     *      F(root) = sum_v W_out(v).
     */
    double root_total_weight = 0.0;

    stack.clear();

    /*
     * The root has no parent edge, so start with its forest children.
     */
    forest_edges.for_each_neighbor(
        root,
        edge_index,
        [&](int const root_child) {
            stack.push({
                root_child,
                root,
                false
            });
        }
    );

    while (!stack.empty()) {
        auto const popped = stack.pop();

        int const vtx =
            std::get<0>(popped);

        int const parent =
            std::get<1>(popped);

        bool const is_revisiting =
            std::get<2>(popped);

        if (!is_revisiting) {
            /*
             * First visit.
             *
             * Push vtx back so it will be processed after all of its
             * children have been processed.
             */
            stack.push({
                vtx,
                parent,
                true
            });

            /*
             * Descend to all forest neighbors except the parent.
             */
            forest_edges.for_each_neighbor(
                vtx,
                edge_index,
                [&](int const child_vtx) {
                    if (child_vtx == parent) {
                        return;
                    }

                    stack.push({
                        child_vtx,
                        vtx,
                        false
                    });
                }
            );

            continue;
        }

        /*
         * ------------------------------------------------------------
         * POSTORDER REVISIT
         * ------------------------------------------------------------
         *
         * Every child has already been processed, so all child
         * pops_below_vertex values are available.
         */

        int pop_below_vtx =
            static_cast<int>(pop[vtx]);

        forest_edges.for_each_neighbor(
            vtx,
            edge_index,
            [&](int const child) {
                if (child == parent) {
                    return;
                }

                pop_below_vtx +=
                    pops_below_vertex[child];
            }
        );

        pops_below_vertex[vtx] =
            pop_below_vtx;

        /*
         * Population of THIS tree outside the subtree below vtx.
         */
        int const this_tree_pop_above_vtx =
            this_region_tree_pop -
            pop_below_vtx;


        /*
         * ============================================================
         * OTHER TREE ATTACHED OUTSIDE S_v
         * ============================================================
         *
         * Cutting (parent, vtx) gives:
         *
         *   below = S_v
         *   above = (this tree \ S_v) + other tree
         */
        int const pop_below_attach_outside =
            pop_below_vtx;

        int const pop_above_attach_outside =
            this_tree_pop_above_vtx +
            other_region_tree_pop;

        bool const out_could_be_valid =
            pop_below_attach_outside >=
                smallest_lower_bound &&
            pop_below_attach_outside <=
                biggest_upper_bound &&
            pop_above_attach_outside >=
                smallest_lower_bound &&
            pop_above_attach_outside <=
                biggest_upper_bound;

        double weight_out = 0.0;

        if (out_could_be_valid) {
            for_each_valid_population_split_from_edge(
                total_region_size,

                /*
                 * IMPORTANT:
                 * for_each_valid_population_split_from_edge takes
                 * BELOW population first, then ABOVE population.
                 */
                static_cast<double>(
                    pop_below_attach_outside
                ),
                static_cast<double>(
                    pop_above_attach_outside
                ),

                lower,
                upper,
                smaller_cut_sizes_to_try,
                [&](
                    int const cut_below_region_size,
                    double const cut_below_pop,
                    int const cut_above_region_size,
                    double const cut_above_pop
                ) {
                    weight_out +=
                        compute_cut_weight(
                            root,
                            vtx,
                            parent,
                            cut_below_region_size,
                            cut_below_pop,
                            cut_above_region_size,
                            cut_above_pop
                        );
                }
            );
        }


        /*
         * ============================================================
         * OTHER TREE ATTACHED INSIDE S_v
         * ============================================================
         *
         * Cutting (parent, vtx) gives:
         *
         *   below = S_v + other tree
         *   above = this tree \ S_v
         */
        int const pop_below_attach_inside =
            pop_below_vtx +
            other_region_tree_pop;

        int const pop_above_attach_inside =
            this_tree_pop_above_vtx;

        bool const in_could_be_valid =
            pop_below_attach_inside >=
                smallest_lower_bound &&
            pop_below_attach_inside <=
                biggest_upper_bound &&
            pop_above_attach_inside >=
                smallest_lower_bound &&
            pop_above_attach_inside <=
                biggest_upper_bound;

        double weight_in = 0.0;

        if (in_could_be_valid) {
            for_each_valid_population_split_from_edge(
                total_region_size,
                static_cast<double>(
                    pop_below_attach_inside
                ),
                static_cast<double>(
                    pop_above_attach_inside
                ),
                lower,
                upper,
                smaller_cut_sizes_to_try,
                [&](
                    int const cut_below_region_size,
                    double const cut_below_pop,
                    int const cut_above_region_size,
                    double const cut_above_pop
                ) {
                    weight_in +=
                        compute_cut_weight(
                            root,
                            vtx,
                            parent,
                            cut_below_region_size,
                            cut_below_pop,
                            cut_above_region_size,
                            cut_above_pop
                        );
                }
            );
        }

        /*
         * Temporarily use reroot_weight[vtx] to store
         *
         *     Delta(v) = W_in(v) - W_out(v).
         *
         * Pass 2 will overwrite this with F(v).
         */
        reroot_weight[vtx] =
            weight_in -
            weight_out;

        /*
         * When the other tree is attached at root, it is outside every
         * proper rooted subtree S_v. Therefore every internal edge is in
         * its OUT state:
         *
         *     F(root) = sum_{v != root} W_out(v).
         */
        root_total_weight +=
            weight_out;
    }


    /*
     * ============================================================
     * PASS 2: PREORDER
     * ============================================================
     *
     * Convert
     *
     *     reroot_weight[v] = Delta(v)
     *
     * into
     *
     *     reroot_weight[v] = F(v)
     *
     * using
     *
     *     F(child) = F(parent) + Delta(child).
     */


    /*
     * The root value is known directly from pass 1.
     */
    reroot_weight[root] = root_total_weight;

    /*
     * The root's F value is now final.
     *
     * This callback lets the caller immediately use F(root), e.g.
     * inspect map-graph neighbors and calculate effective boundary
     * contributions.
     */
    emit_finished_vertex(root);

    stack.clear();

    /*
     * Start preorder traversal with the children of root.
     */
    forest_edges.for_each_neighbor(
        root,
        edge_index,
        [&](int const root_child) {
            stack.push({root_child, root, false});
        }
    );

    while (!stack.empty()) {
        auto const popped = stack.pop();

        int const vtx = std::get<0>(popped);
        int const parent = std::get<1>(popped);

        /*
         * Before this assignment:
         *
         *     reroot_weight[parent] = F(parent)
         *     reroot_weight[vtx]    = Delta(vtx)
         *
         * because the parent has already been processed but vtx has not.
         */
        double const delta_vtx = reroot_weight[vtx];

        /*
         * Apply the reroot recurrence:
         *
         *     F(vtx) = F(parent) + Delta(vtx).
         */
        reroot_weight[vtx] = reroot_weight[parent] + delta_vtx;

        /*
         * reroot_weight[vtx] is now the FINAL F(vtx) value.
         *
         * Call the user-supplied hook before moving on.
         */
        emit_finished_vertex(vtx);

        /*
         * Continue preorder traversal to children.
         */
        forest_edges.for_each_neighbor(
            vtx,
            edge_index,
            [&](int const child) {
                if (child == parent) {
                    return;
                }

                stack.push({child, vtx, false});
            }
        );
    }
}



/*
 * Appends all valid edge cuts in the tree to existing_cuts
 *
 *
 * Returns a vector of all the valid edge cuts (ie an edge and regions for
 * the two cuts) where at least one of the regions is between
 * `min_potential_cut_size` and `max_potential_cut_size` inclusive. Returns
 * the edges as EdgeCut objects. An empty vector means there are no
 * valid edges.
 *
 *
 *
 * @param ust A directed edge spanning tree.
 * @param root The root vertex of the spanning tree.
 * @param cut_below_pops The population corresponding to cutting below each vertex.
 * So `cut_below_pops[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param min_potential_cut_size The smallest potential region size at least one of
 * the regions cut must be
 * @param max_potential_cut_size The largest potential region size at least one of
 * the regions cut must be. Setting this to 1 will result in only 1 district splits.
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details Adds valid cuts to existing_cuts
 *
 * @return A vector of EdgeCut objects
 *
 */



// This sums up the total unnormalized probability of all possible edge cuts in the 
// tree. In the retroactive splitting probability calculation this is the denominator 
// NOTE this assumes the splitting weight can be written as a of just the population and
// vertices induced by the removed edge 
double get_unnormed_weight_sum_in_joined_packed_tree(
    MapParams const &map_params,
    EdgeBitset const &forest_edges,
    TreePopStack &stack,
    std::vector<int> &pops_below_vertex,
    std::vector<bool> &no_valid_edges_vertices,
    TreeSplitter const &tree_splitter,
    int const region1_root,
    int const region1_pop,
    int const region2_root,
    int const region2_pop,
    int const min_potential_cut_size,
    int const max_potential_cut_size,
    std::vector<int> const
        &smaller_cut_sizes_to_try,
    int const total_merged_region_size
) {
    double unnormed_weight_sum = 0.0;

    for_each_valid_edge_cut_in_joined_packed_tree(
        map_params,
        forest_edges,
        stack,
        pops_below_vertex,
        no_valid_edges_vertices,
        region1_root,
        region1_pop,
        region2_root,
        region2_pop,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_size,
        [&](
            int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
        ) {
            unnormed_weight_sum +=
                (
                    tree_splitter
                        .get_unnormed_selection_prob(
                            root,
                            cut_vertex,
                            cut_vertex_parent,
                            cut_below_region_size,
                            cut_below_pop,
                            cut_above_region_size,
                            cut_above_pop
                        )
                );
        }
    );

    return unnormed_weight_sum;
}

double compute_signed_pop_deviance(
    double const target, 
    int const region_pop,
    int const region_size
) {
    // get the target populations for the regions
    double region_target = target * region_size;
    // get the deviation
    return (static_cast<double>(region_pop) - region_target) / region_target;
}

double compute_absolute_pop_deviance(
    double const target, 
    int const region_pop,
    int const region_size
) {
    // get the target populations for the regions
    double region_target = target * region_size;
    // get the deviation
    return std::fabs(compute_signed_pop_deviance(target, region_pop, region_size));
}

arma::vec compute_almost_best_weights_on_smaller_dev_edges(std::vector<EdgeCut> &valid_edges,
                                                           double epsilon, double target) {

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    // find the maximum value
    double global_min = 42.0;

    for (size_t i = 0; i < valid_edges.size(); i++) {
        std::array<double, 2> devs = valid_edges.at(i).compute_abs_pop_deviances(target);
        double smaller_dev = std::min(devs.at(0), devs.at(1));
        unnormalized_wgts(i) = smaller_dev;
        // Rprintf("Bigger abs dev = %.3f, Computed weight %.3f\n",
        //     smaller_dev, unnormalized_wgts(i));

        global_min = std::min(global_min, smaller_dev);

        // Rprintf("devs are (%.6f,%.6f),  Best so far is %.6f\n",
        //     devs.at(0), devs.at(1), global_min);
    }
    // Rprintf("\n\n");

    for (size_t i = 0; i < valid_edges.size(); i++) {
        // make 1 if eqaul to the max, epsilon otherwise
        // REprintf("Set Weight %d, dev %f to %f \n",
        //     (int) i, unnormalized_wgts(i),
        //     (unnormalized_wgts(i) == global_min) ? 1.0 : epsilon);
        unnormalized_wgts(i) = (unnormalized_wgts(i) == global_min) ? 1.0 : epsilon;
    }

    return unnormalized_wgts;
}

}



constexpr bool MERGED_TREE_SPLITTING_VERBOSE = false; // Compile-time constant

std::vector<EdgeCut> TreeSplitter::get_all_valid_pop_edge_cuts_in_directed_tree(
    const MapParams &map_params, Tree const &ust, // FlatGraph const &ust, 
    const int root, 
    std::vector<int> &pops_below_vertex, 
    int const region_population, int const region_size, const int min_potential_cut_size,
    const int max_potential_cut_size, std::vector<int> const &smaller_cut_sizes_to_try) const {

    
    // Don't need to reset pops below as it starts from the child nodes 
    // so it overwrites old values 
    // Also don't need to reset no_valid_edges_vertices as we 
    // can just set them to false the first time we see the vertex when 
    // traversing the tree 

    // reset pops_below_vertex and valid edges thing
    std::vector<EdgeCut> valid_edges = get_all_valid_edges_in_directed_tree(
        ust, root, map_params.pop, stack, pops_below_vertex, no_valid_edges_vertices,
        min_potential_cut_size, max_potential_cut_size, smaller_cut_sizes_to_try,
        region_population, region_size, map_params.lower, map_params.upper, map_params.target);

    return valid_edges;
}

std::pair<bool, EdgeCut> TreeSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
    Plan const &plan, int const split_region1, int const split_region2, Tree const &ust, // FlatGraph const &ust,
    const int root, std::vector<int> &pops_below_vertex,
    int const region_population,
    int const region_size, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, bool save_selection_prob) {
    // get all the valid edges
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, ust, root, pops_below_vertex, 
        region_population, region_size, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if (num_valid_edges == 0) {
        return std::make_pair(false, EdgeCut());
    } else { // else have derived class choose according to its rule
        return select_edge_to_cut(rng_state, valid_edges,
                                  save_selection_prob);
    }
}

// returns edge cut and log probability it was chosen
std::pair<bool, EdgeCut>
TreeSplitter::select_edge_to_cut(
                                 RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
                                 bool save_selection_prob) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately
    if (num_valid_edges == 1) {
        // if only 1 just return that
        // selection prob is just 1 so don't touch
        // if(save_selection_prob){
        //     Rprintf("Save true: %d valid, only 1 edge, log prob is %f \n",
        //         num_valid_edges, valid_edges[0].log_prob);
        // }
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights
    arma::vec unnormalized_wgts(num_valid_edges);

    for (size_t i = 0; i < num_valid_edges; i++) {
        unnormalized_wgts(i) = compute_unnormalized_edge_cut_weight(valid_edges[i]);
    }

    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    if (save_selection_prob) {
        selected_edge_cut.log_prob =
            std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
        // Rprintf("Save, %d valid, log prob is %f and %f\n", num_valid_edges,
        // selected_edge_cut.log_prob,
        //     std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts)));
    }

    return std::make_pair(true, selected_edge_cut);
}

// Takes a vector of valid edge cuts and returns the log probability
// the one an index idx would have been chosen
double TreeSplitter::get_log_selection_prob(std::vector<EdgeCut> &valid_edges, int idx) const {
    auto num_valid_edges = valid_edges.size();
    // get the weights
    double weight_sum = 0.0;
    // get idx weight
    double idx_weight = compute_unnormalized_edge_cut_weight(valid_edges[idx]);

    // get sum of weights
    for (size_t i = 0; i < num_valid_edges; i++) {
        weight_sum += compute_unnormalized_edge_cut_weight(valid_edges[i]);
    }

    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return std::log(idx_weight) - std::log(weight_sum);
}

double TreeSplitter::get_log_retroactive_splitting_prob_for_joined_packed_tree(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    EdgeBitset const &forest_edges, 
    std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
    Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try) {
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];

    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size + region2_size;

    // Get the total unnormed selection probability in the whole tree 
    // made by joining region1_root to region2_root
    long double const denominator =
        get_unnormed_weight_sum_in_joined_packed_tree(
            map_params,
            forest_edges,
            stack,
            pops_below_vertex,
            visited,
            *this,
            region1_root,
            region1_population,
            region2_root,
            region2_population,
            min_potential_cut_size,
            max_potential_cut_size,
            smaller_cut_sizes_to_try,
            total_merged_region_size
        );

    // get the unnormed selection prob for (region1_root, region2_root)
    double const numerator =
        get_unnormed_selection_prob(
            region1_root,
            region2_root,
            region1_root,
            region2_size,
            region2_population,
            region1_size,
            region1_population
        );

    return std::log(numerator) - std::log(denominator);
}


double TreeSplitter::get_log_eff_boundary_len_for_adj_region_pair(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    EdgeBitset const &forest_edges, 
    std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
    Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    bool const count_edges_across
) {
    /*
     * Get the region IDs corresponding to the supplied representative
     * boundary-edge endpoints.
     */
    RegionID const region1_id = plan.region_ids[region1_root];
    RegionID const region2_id = plan.region_ids[region2_root];

    /*
     * Region-level information.
     */
    int const region1_population = plan.region_pops[region1_id];
    int const region2_population = plan.region_pops[region2_id];
    int const region1_size = plan.region_sizes[region1_id];
    int const region2_size = plan.region_sizes[region2_id];
    int const total_merged_region_size = region1_size + region2_size;

   /*
    * Imagine the map edge
    *
    *     region1_root -- region2_root
    *
    * with region1_root in region 1 and region2_root in region 2.
    *
    * We orient the joined tree so that:
    *
    *     tree root         = region1_root
    *     cut vertex        = region2_root
    *     cut vertex parent = region1_root
    *
    * Therefore, cutting the joining edge gives:
    *
    *     below = region 2
    *     above = region 1.
    *
    * A single physical joining edge may correspond to more than
    * one valid size assignment when population bounds are loose,
    * so sum the weights of ALL valid cuts associated with this
    * physical edge.
    */
    double joining_edge_weight = 0.0;

    for_each_valid_population_split_from_edge(
        total_merged_region_size,
        static_cast<double>(region2_population),
        static_cast<double>(region1_population),
        map_params.lower,
        map_params.upper,
        smaller_cut_sizes_to_try,
        [&](int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop) {
            joining_edge_weight +=
                get_unnormed_selection_prob(
                    region1_root,
                    region2_root,
                    region1_root,
                    cut_below_region_size,
                    cut_below_pop,
                    cut_above_region_size,
                    cut_above_pop
                );
        }
    );

    // We now compute the unnormalized selection probability for 
    // the edge (region1_root, region2_root) with the specific sizes
   /*
    * --------------------------------------------------------
    * NUMERATOR
    * --------------------------------------------------------
    *
    * The actual cut whose retroactive probability we want is
    * the cut which recovers the CURRENT two regions:
    *
    *     below = region 2
    *     above = region 1.
    *
    * This is one particular valid size assignment on the
    * joining edge, not the sum over all valid assignments.
    */
    double const numerator = get_unnormed_selection_prob(
        region1_root,
        region2_root,
        region1_root,
        region2_size,
        region2_population,
        region1_size,
        region1_population
    );


    /*
     * This callback gives the unnormalized selection weight of one valid
     * candidate cut.
     *
     * compute_tree_path_weights_in_undirected_packed_tree uses this when
     * computing W_out(v) and W_in(v) for every internal tree edge.
     */
    auto const compute_cut_weight =
        [&](int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop) {

            return get_unnormed_selection_prob(
                root,
                cut_vertex,
                cut_vertex_parent,
                cut_below_region_size,
                cut_below_pop,
                cut_above_region_size,
                cut_above_pop
            );
        };

    /*
     * ================================================================
     * REGION 2
     * ================================================================
     *
     * Compute F_2(v) for every vertex in region 2 first.
     *
     * We do region 2 first because, while computing region 1 below, we
     * will inspect map edges from region-1 vertices to region-2 vertices.
     * At that point we need the region-2 F values to already be available.
     *
     * Since the two regions have disjoint vertex sets, computing region 1
     * afterward does not overwrite any of these region-2 entries in
     * reroot_weight.
     */
    compute_tree_path_weights_in_undirected_packed_tree(
        map_params.graph_edge_index,
        forest_edges,
        region2_root,
        map_params.pop,
        stack,
        pops_below_vertex,
        reroot_weight,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        region2_population,
        region1_population,
        total_merged_region_size,
        map_params.lower,
        map_params.upper,
        compute_cut_weight,

        /*
         * We do not need to do anything as region-2 vertices finish.
         * Their final F_2(v) values simply remain in reroot_weight.
         */
        [](int const) {}
    );

    /*
     * This will accumulate
     *
     *     sum_{(a,b) in allowed boundary edges}
     *         P(select (a,b) | joined at (a,b)).
     *
     * That is the effective tree boundary length between these regions.
     */
    double effective_boundary_len = 0.0;

    /*
     * ================================================================
     * REGION 1
     * ================================================================
     *
     * Compute F_1(v) for every vertex in region 1.
     *
     * As soon as F_1(v) is finalized, inspect all MAP-GRAPH neighbors of
     * v. For every allowed neighbor in region 2 we can immediately compute
     * the retroactive splitting probability associated with joining the
     * two region trees at that map edge.
     */
    compute_tree_path_weights_in_undirected_packed_tree(
        map_params.graph_edge_index,
        forest_edges,
        region1_root,
        map_params.pop,
        stack,
        pops_below_vertex,
        reroot_weight,
        min_potential_cut_size,
        max_potential_cut_size,
        smaller_cut_sizes_to_try,
        region1_population,
        region2_population,
        total_merged_region_size,
        map_params.lower,
        map_params.upper,
        compute_cut_weight,

        [&](int const v) {
            /*
             * At this point:
             *
             *     reroot_weight[v] = F_1(v)
             *
             * and for every vertex u in region 2:
             *
             *     reroot_weight[u] = F_2(u).
             *
             * Now inspect map-graph edges out of v.
             */
            for (int const nbor : map_params.g[v]) {

                /*
                 * Only map edges from region 1 to region 2 contribute to
                 * this particular region pair's effective boundary.
                 */
                if (plan.region_ids[nbor] != region2_id) {
                    continue;
                }

                /*
                 * Under the hierarchical rule, if across-county edges are
                 * not allowed for this region pair, only count boundary
                 * edges whose endpoints are in the same county.
                 */
                if (
                    !count_edges_across &&
                    map_params.counties[v] !=
                        map_params.counties[nbor]
                ) {
                    continue;
                }

                /*
                 * The full denominator for the joined tree is:
                 *
                 *     internal contribution from region 1
                 *   + internal contribution from region 2
                 *   + joining-edge contribution.
                 *
                 * By construction:
                 *
                 *     reroot_weight[v]    = F_1(v)
                 *     reroot_weight[nbor] = F_2(nbor).
                 */
                double const denominator =
                    reroot_weight[v] +
                    reroot_weight[nbor] +
                    joining_edge_weight;

                /*
                 * Probability that the splitter would select the actual
                 * joining cut if the two region trees were joined using
                 * map edge (v, nbor).
                 */
                double const edge_selection_prob = numerator/ denominator;

                /*
                 * Sum these probabilities over all allowed map boundary
                 * edges between region 1 and region 2.
                 *
                 * There is no double counting here because this callback
                 * runs only while traversing region 1. We therefore see
                 * each undirected region1-region2 map edge only from its
                 * region-1 endpoint.
                 */
                effective_boundary_len += edge_selection_prob;
            }
        }
    );

    /*
     * Return the log effective boundary length, matching the existing
     * ForestPlan interface.
     */
    return std::log(effective_boundary_len);
}


double TreeSplitter::get_log_retroactive_splitting_prob_from_valid_pop_cut_list(
        std::vector<EdgeCut> &valid_edges, EdgeCut const actual_cut_edge
){
    if (MERGED_TREE_SPLITTING_VERBOSE) {
        int region1_root = actual_cut_edge.tree_root;
        int region2_root = actual_cut_edge.cut_vertex;
        REprintf("Finding Merge prob for (%d, %d) - %zu valid edges!\n", region1_root,
                 region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    if constexpr (perf_config::bounds_checking){
        if (it == valid_edges.end()) {
            int region1_root = actual_cut_edge.tree_root;
            int region2_root = actual_cut_edge.cut_vertex;
            int region1_size = actual_cut_edge.cut_above_region_size;
            int region1_population = actual_cut_edge.cut_above_pop;
            int region2_size = actual_cut_edge.cut_below_region_size;
            int region2_population = actual_cut_edge.cut_below_pop;
            std::ostringstream oss;

            oss << "Actual cut edge not found in valid_edges in "
                << "get_log_retroactive_splitting_prob_from_valid_pop_cut_list.\n";

            oss << "region1_root=" << region1_root << "\n";
            oss << "region2_root=" << region2_root << "\n";
            oss << "region1_size=" << region1_size << "\n";
            oss << "region2_size=" << region2_size << "\n";
            oss << "region1_population=" << region1_population << "\n";
            oss << "region2_population=" << region2_population << "\n";
            oss << "valid_edges.size()=" << valid_edges.size() << "\n";

            oss << "Actual cut edge: "
                << "tree_root=" << actual_cut_edge.tree_root
                << ", cut_vertex=" << actual_cut_edge.cut_vertex
                << ", cut_vertex_parent=" << actual_cut_edge.cut_vertex_parent
                << ", cut_below_region_size=" << actual_cut_edge.cut_below_region_size
                << ", cut_below_pop=" << actual_cut_edge.cut_below_pop
                << ", cut_above_region_size=" << actual_cut_edge.cut_above_region_size
                << ", cut_above_pop=" << actual_cut_edge.cut_above_pop
                << "\n";

            oss << "Valid edges:\n";
            for (std::size_t k = 0; k < valid_edges.size(); ++k) {
                auto const &e = valid_edges[k];

                oss << "  [" << k << "] "
                    << "tree_root=" << e.tree_root
                    << ", cut_vertex=" << e.cut_vertex
                    << ", cut_vertex_parent=" << e.cut_vertex_parent
                    << ", cut_below_region_size=" << e.cut_below_region_size
                    << ", cut_below_pop=" << e.cut_below_pop
                    << ", cut_above_region_size=" << e.cut_above_region_size
                    << ", cut_above_pop=" << e.cut_above_pop
                    << ", log_prob=" << e.log_prob
                    << "\n";
            }

            throw std::runtime_error(oss.str());
        }
    }

    if (it == valid_edges.end()) {
        std::ostringstream oss;
        oss << "Actual cut edge not found in retroactive "
            << "valid-edge list.\n";
        oss << "valid_edges.size()="
            << valid_edges.size() << "\n";

        throw std::runtime_error(oss.str());
    }

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    
    if constexpr (perf_config::bounds_checking){
        if (actual_cut_edge_index < 0 ||
        actual_cut_edge_index >= static_cast<int>(valid_edges.size())) {
            std::ostringstream oss;
            oss << "actual_cut_edge_index out of bounds. "
                << "actual_cut_edge_index=" << actual_cut_edge_index
                << ", valid_edges.size()=" << valid_edges.size();

            throw std::runtime_error(oss.str());
        }
    }

    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Actual Cut Edge at Index %d and so prob is %f \n", actual_cut_edge_index,
                 get_log_selection_prob(valid_edges, actual_cut_edge_index));
    }

    return get_log_selection_prob(valid_edges, actual_cut_edge_index);
}

void NaiveTopKSplitter::update_single_int_param(int int_param) {
    if (int_param <= 0) {
        std::ostringstream oss;
        oss << "Splitting k must be at least 1.\n";
        oss << "Received int_param=" << int_param << "\n";

        throw std::runtime_error(oss.str());
    }
    k_param = int_param;
}

std::pair<bool, EdgeCut>
NaiveTopKSplitter::select_edge_to_cut(RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
                                      bool save_selection_prob) const {

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if(num_valid_edges > k_param){
    //     REprintf("k was %d but found %d valid edges\n", k_param, num_valid_edges);
    //     // throw Rcpp::exception("K not big enough!\n");
    // }

    int idx = rng_state.r_int(k_param);
    // if we selected k greater than number of edges failure
    if (idx >= num_valid_edges) {
        return std::make_pair(false, EdgeCut());
    } else {
        // we always store selection probability since its so cheap to compute
        EdgeCut selected_edge_cut = valid_edges[idx];
        selected_edge_cut.log_prob = -std::log(k_param);
        return std::make_pair(true, selected_edge_cut);
    }
}

std::pair<bool, EdgeCut> UniformValidSplitter::select_edge_to_cut(
    RNGState &rng_state,
    std::vector<EdgeCut> &valid_edges, bool save_selection_prob) const {
    std::uint32_t num_valid_edges = static_cast<std::uint32_t>(valid_edges.size());
    // if only 1 edge just return that
    if (num_valid_edges == 1)
        return std::make_pair(true, valid_edges[0]);

    // pick one unif at random
    int idx = rng_state.r_int(num_valid_edges);
    // we always store selection probability since its so cheap to compute
    EdgeCut selected_edge_cut = valid_edges[idx];
    selected_edge_cut.log_prob = -std::log(num_valid_edges);

    return std::make_pair(true, selected_edge_cut);
}

double 
UniformValidSplitter::get_unnormed_selection_prob(
            int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
) const {
    // since uniform over edges just return 1 
    return 1.0;
}




/*
 * Assigns weight to a potential split of exp(-alpha*larger abs dev)
 *
 * Assigns weight to a potential split of exp(-alpha*larger dev) where 
 * larger abs dev is the bigger absolute deviation from the target of 
 * the two regions induced by the cut. 
 */
double ExpoWeightedSplitter::get_unnormed_selection_prob(
            int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
) const {
    // get the absolute deviation
    double const above_dev = compute_absolute_pop_deviance(target, cut_above_pop, cut_above_region_size);
    double const below_dev = compute_absolute_pop_deviance(target, cut_below_pop, cut_below_region_size);
    // take the bigger of them 
    double bigger_dev = std::max(above_dev, below_dev);
    // return the value 
    return std::exp(-alpha * bigger_dev);
}

double
ExpoWeightedSplitter::compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const {
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double bigger_dev = std::max(devs.at(0), devs.at(1));
    return std::exp(-alpha * bigger_dev);
}

/*
 * Assigns weight to a potential split of exp(-alpha*smaller abs dev)
 *
 * Assigns weight to a potential split of exp(-alpha*smller dev) where 
 * smaller abs dev is the smaller absolute deviation from the target of 
 * the two regions induced by the cut. 
 */
double ExpoWeightedSmallerDevSplitter::get_unnormed_selection_prob(
            int const root,
            int const cut_vertex,
            int const cut_vertex_parent,
            int const cut_below_region_size,
            double const cut_below_pop,
            int const cut_above_region_size,
            double const cut_above_pop
) const {
    // get the absolute deviation
    double const above_dev = compute_absolute_pop_deviance(target, cut_above_pop, cut_above_region_size);
    double const below_dev = compute_absolute_pop_deviance(target, cut_below_pop, cut_below_region_size);
    // take the smaller of them 
    double smaller_dev = std::min(above_dev, below_dev);
    // return the value 
    return std::exp(-alpha * smaller_dev);
}

double ExpoWeightedSmallerDevSplitter::compute_unnormalized_edge_cut_weight(
    EdgeCut const &edge_cut) const {
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double smaller_dev = std::min(devs.at(0), devs.at(1));
    return std::exp(-alpha * smaller_dev);
}

double PopTemperSplitter::compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const {
    double region1_pop_temper = compute_log_pop_temper(
        target, pop_temper, ndists, edge_cut.cut_above_pop, edge_cut.cut_above_region_size);
    double region2_pop_temper = compute_log_pop_temper(
        target, pop_temper, ndists, edge_cut.cut_below_pop, edge_cut.cut_below_region_size);
    // larger population deviation means bigger pop temper but we want smaller
    // so we add then do exp(- sum)
    return std::exp(-(region1_pop_temper + region2_pop_temper));
}

std::pair<bool, EdgeCut> ExperimentalSplitter::select_edge_to_cut(
    RNGState &rng_state,
    std::vector<EdgeCut> &valid_edges, bool save_selection_prob) const {
    auto num_valid_edges = valid_edges.size();

    // if no valid edges reject immediately
    if (num_valid_edges == 1) {
        // if only 1 just return that
        return std::make_pair(true, valid_edges[0]);
    }

    // get the weights
    arma::vec unnormalized_wgts =
        compute_almost_best_weights_on_smaller_dev_edges(valid_edges, epsilon, target);

    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    if (save_selection_prob) {
        selected_edge_cut.log_prob =
            std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
    }

    return std::make_pair(true, selected_edge_cut);
}

double ExperimentalSplitter::get_log_selection_prob(std::vector<EdgeCut> &valid_edges,
                                                    int idx) const {
    // get the weights
    arma::vec unnormalized_wgts =
        compute_almost_best_weights_on_smaller_dev_edges(valid_edges, epsilon, target);

    // we want log of weight at idx / sum of all weight which is equal to
    // log(prob at idx) - log(sum of all weights)
    return log(unnormalized_wgts(idx)) - log(arma::sum(unnormalized_wgts));
}


std::pair<bool, EdgeCut> ConstraintSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
    Plan const &plan, int const split_region1, int const split_region2, Tree const &ust, // FlatGraph const &ust,
    const int root, std::vector<int> &pops_below_vertex,
    int const region_population,
    int const region_size, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, bool save_selection_prob) {
    // get all the valid edges
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, ust, root, pops_below_vertex, 
        region_population, region_size, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if (num_valid_edges == 0) {
        return std::make_pair(false, EdgeCut());
    } else if (num_valid_edges == 1 && !scoring_function.any_hard_constraints) {
        return std::make_pair(true, valid_edges[0]);
    }


    // copy over the current plan information
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);

    // get the weights
    arma::vec unnormalized_wgts = compute_constraint_edge_cut_weights(
        valid_edges, scoring_function, ust, plan.num_regions + 1, region_ids, region_sizes,
        region_pops, split_region1, split_region2, vertex_queue);

    // check if all filtered
    num_valid_edges = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if (num_valid_edges == 0) {
        return std::make_pair(false, EdgeCut());
    }

    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    if (save_selection_prob) {
        selected_edge_cut.log_prob =
            std::log(unnormalized_wgts[idx]) - std::log(
                std::accumulate(unnormalized_wgts.begin(), unnormalized_wgts.end(), 0.0)
            );
        // REprintf("Selection prob %f\n", selected_edge_cut.log_prob);
    }


    return std::make_pair(true, selected_edge_cut);
}

// assumes two trees in spanning forest have been joined
// at the edge (cut_vertex_root_region_id, cut_vertex_parent)
void assign_region_ids_from_joined_undirected_tree(
    Tree const &forest_graph, PlanVector &region_ids, int const cut_vertex_root,
    int const cut_vertex_root_region_id, int const cut_vertex_parent,
    int const cut_parent_region_id, CircularQueue<std::pair<int, int>> &vertex_queue) {
    // clear the queue
    vertex_queue.clear();

    // Since tree is undirected we first start from cut_vertex_root
    // and just make sure to skip the parent

    // update root and add its children to queue
    region_ids[cut_vertex_root] = cut_vertex_root_region_id;
    for (auto const &child_vertex : forest_graph[cut_vertex_root]) {
        // ignore if its the parent aka the cut edge
        if (child_vertex == cut_vertex_parent)
            continue;

        vertex_queue.push({child_vertex, cut_vertex_root});
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, vtx_parent] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = cut_vertex_root_region_id;
        // add children
        for (auto const &child_vertex : forest_graph[vertex]) {
            // if its the parent then skip it
            if (child_vertex == vtx_parent)
                continue;
            vertex_queue.push({child_vertex, vertex});
        }
    }

    // now we update starting at the cut root vertex
    // update root and add its children to queue
    region_ids[cut_vertex_parent] = cut_parent_region_id;
    for (auto const &child_vertex : forest_graph[cut_vertex_parent]) {
        // ignore if its the child aka the cut edge
        if (child_vertex == cut_vertex_root)
            continue;

        vertex_queue.push({child_vertex, cut_vertex_parent});
    }

    // update all the children
    while (!vertex_queue.empty()) {
        // get and remove head of queue
        auto [vertex, parent_vtx] = vertex_queue.pop();
        // update region ids
        region_ids[vertex] = cut_parent_region_id;
        // add children
        for (auto const &child_vertex : forest_graph[vertex]) {
            // if its the parent then skip it
            if (child_vertex == parent_vtx)
                continue;
            vertex_queue.push({child_vertex, vertex});
        }
    }

    return;
}

double ConstraintSplitter::get_log_retroactive_splitting_prob_for_joined_packed_tree(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    EdgeBitset const &forest_edges, 
    std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
    Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try) {
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];

    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size + region2_size;

    auto const region1_id = plan.region_ids[region1_root];
    auto const region2_id = plan.region_ids[region2_root];

    // Get all the valid edges in the joined tree
    std::vector<EdgeCut> valid_edges = get_valid_pop_edges_in_joined_packed_tree(
        map_params, forest_edges, stack, pops_below_vertex, visited, region1_root,
        region1_population, region2_root, region2_population, min_potential_cut_size,
        max_potential_cut_size, smaller_cut_sizes_to_try, total_merged_region_size);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if only 1 valid edge then its log(1) = 0
    if (num_valid_edges == 1 && !scoring_function.any_hard_constraints) {
        return 0.0;
    }

    // copy the forest over to the dummy forest 
    dummy_forest = forest_edges.get_graph_tree(map_params.graph_edge_index);

    // find the index of the actual edge we cut
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(region1_root, region2_root, region1_root, region2_size,
                            region2_population, region1_size, region1_population);

    // Now return the probability we actually selected that cut edge 
    return custom_get_log_retroactive_splitting_prob_from_valid_pop_cut_list(
        valid_edges, actual_cut_edge, map_params, scoring_function, plan
    );
}


double ConstraintSplitter::custom_get_log_retroactive_splitting_prob_from_valid_pop_cut_list(
    std::vector<EdgeCut> &valid_edges, EdgeCut const actual_cut_edge,
    MapParams const &map_params, ScoringFunction const &scoring_function,
    Plan const &plan
){
    int const region1_root = actual_cut_edge.tree_root;
    int const region2_root = actual_cut_edge.cut_vertex;
    auto const region1_id = plan.region_ids[region1_root];
    auto const region2_id = plan.region_ids[region2_root];
    // copy over the current plan information
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);


    // add the actual removed edge back
    dummy_forest[region1_root].push_back(region2_root);
    dummy_forest[region2_root].push_back(region1_root);

    std::vector<long double> unnormed_wgts;
    unnormed_wgts.reserve(valid_edges.size());

    std::size_t write_idx = 0;

    for (std::size_t read_idx = 0; read_idx < valid_edges.size(); ++read_idx) {
        EdgeCut const &edge_cut = valid_edges[read_idx];

        // Update split information.
        region_sizes[region1_id] = edge_cut.cut_above_region_size;
        region_sizes[region2_id] = edge_cut.cut_below_region_size;

        region_pops[region1_id] = edge_cut.cut_above_pop;
        region_pops[region2_id] = edge_cut.cut_below_pop;

        // Update the region IDs.
        assign_region_ids_from_joined_undirected_tree(
            dummy_forest,
            region_ids,
            edge_cut.cut_vertex,
            region2_id,                 // below side
            edge_cut.cut_vertex_parent,
            region1_id,                 // above side
            vertex_queue
        );

        std::pair<bool, double> const score_result =
            scoring_function.compute_full_split_plan_score(
                plan.num_regions,
                region_ids,
                region_sizes,
                region_pops,
                region1_id,
                region2_id
            );

        // Remove cuts that violate a hard constraint by not copying them
        // into the retained portion of valid_edges.
        if (!score_result.first) {
            continue;
        }

        if (write_idx != read_idx) {
            valid_edges[write_idx] = std::move(valid_edges[read_idx]);
        }

        unnormed_wgts.push_back(
            std::exp(-score_result.second)
        );

        ++write_idx;
    }

    // Remove the rejected edge cuts from the end of the vector.
    valid_edges.resize(write_idx);



    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Finding Merge prob for (%d, %d) - %zu valid edges!\n", region1_root,
                 region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    if (it == valid_edges.end()) {
        // The actual cut was rejected by a hard constraint or was otherwise absent.
        std::ostringstream oss;
        oss << "Error in ConstraintSplitter. In get_log_retroactive_splitting_prob_for_joined_packed_tree "
            << "there were no valid tree cuts which is theoretically impossible unless a non-splittable"
            << " score function was supplied!\n";
        throw std::runtime_error(oss.str());
    }
    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);

    auto sum = std::accumulate(unnormed_wgts.begin(), unnormed_wgts.end(), 0.0);
    auto log_sum = std::log(sum);
    // REprintf("Sum %.30f\n", sum);

    // compute selection probability if needed
    double log_selection_prob = std::log(unnormed_wgts[actual_cut_edge_index]) - log_sum;


    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Actual Cut Edge at Index %d and so prob is %f \n", actual_cut_edge_index,
                 log_selection_prob);
    }

    return log_selection_prob;

    // get the weights
    // arma::vec unnormalized_wgts(unnormed_wgts.size());

    // for (size_t i = 0; i < unnormed_wgts.size(); i++) {
    //     unnormalized_wgts[i] = std::exp(std::log(unnormed_wgts[i]) - log_sum);
    //     REprintf("%.30f\n", unnormalized_wgts[i]);
    // }

    // unnormalized_wgts = arma::cumsum(unnormalized_wgts);

    // REprintf("Weights are:\n");
    // for (auto const v : unnormalized_wgts) {
    //     REprintf("%.20f, ", v);
    // }
    // REprintf("\n");
}