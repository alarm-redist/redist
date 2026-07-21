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

namespace{

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
                                  double const lower, double const target, double const upper,
                                  std::vector<int> const &smaller_cut_sizes_to_try) {

    // iterate over all possible valid sizes of the smaller region
    for (auto const cut_region1_size : smaller_cut_sizes_to_try) {
        int cut_region2_size = total_region_size - cut_region1_size;
        // Get the bounds for region 1
        double cut_region1_lb = lower * cut_region1_size;
        double cut_region1_ub = upper * cut_region1_size;
        // Get the bounds for region 2
        double cut_region2_lb = lower * cut_region2_size;
        double cut_region2_ub = upper * cut_region2_size;

        if constexpr (FINDING_EDGE_CUTS_VERBOSE) {
            REprintf("\tFor (%d, %d): compare %f vs %f vs %f and %f vs %f vs %f \n",
                     cut_region1_size, cut_region2_size, cut_region1_lb, below_pop,
                     cut_region1_ub, cut_region2_lb, above_pop, cut_region2_ub);
        }

        // check if assigning potential_region_size to cut below leads to valid region
        bool cut_below_ok = cut_region1_lb <= below_pop && below_pop <= cut_region1_ub &&
                            cut_region2_lb <= above_pop && above_pop <= cut_region2_ub;

        if (cut_below_ok) {
            // cut region 1 size is cut below
            int cut_below_region_size = cut_region1_size;
            int cut_above_region_size = cut_region2_size;
            existing_cuts.emplace_back(root, cut_vertex, cut_vertex_parent, cut_below_region_size,
                                     below_pop, cut_above_region_size, above_pop);
        }

        // if both sizes are the same then results are symetric so ignore this case
        if (cut_region1_size == cut_region2_size)
            continue;

        if constexpr (FINDING_EDGE_CUTS_VERBOSE) {
            REprintf("\tFor (%d, %d): compare %f vs %f vs %f and %f vs %f vs %f \n",
                     cut_region2_size, cut_region1_size, cut_region2_lb, below_pop,
                     cut_region2_ub, cut_region1_lb, above_pop, cut_region1_ub);
        }

        // check if assigning potential_region_size to cut below leads to valid region
        bool cut_above_ok = cut_region1_lb <= above_pop && above_pop <= cut_region1_ub &&
                            cut_region2_lb <= below_pop && below_pop <= cut_region2_ub;

        if (cut_above_ok) {
            // cut region 1 size is cut above
            int cut_below_region_size = cut_region2_size;
            int cut_above_region_size = cut_region1_size;
            existing_cuts.emplace_back(root, cut_vertex, cut_vertex_parent, cut_below_region_size,
                                     below_pop, cut_above_region_size, above_pop);
        }
    }

    return;
}

/*
 * Pick one of the valid tree edges to split uniformly at random if possible
 *
 *
 * Returns a valid tree edge to split uniformly at random if at least one
 * valid edge to cut is in the tree. If successful returns information on the
 * edge and region sizes associated with the cut.
 *
 * Note even if a successful cut is found it does not
 * update the plan or the tree.
 *
 *
 * It will only attempt to create regions where the size is between
 * min_potential_d and max_potential_d (inclusive). So the one district
 * split case is `min_potential_d=max_potential_d=1`.
 *
 * Valid edge here is defined as an edge and region sizes such that the
 * two induced regions both fall within the population bounds.
 *
 *
 * @param root The root vertex of the spanning tree
 * @param pop_below The population corresponding to cutting below each vertex.
 * So `pop_below[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
 * means the vertex is the root or it is not in the tree.
 * @param min_potential_cut_size The smallest potential region size to try for a cut.
 * @param max_potential_cut_size The largest potential region size it will try for a cut.
 * Setting this to 1 will result in only 1 district splits.
 * @param region_ids A vector mapping 0 indexed vertices to their region id number
 * @param region_id_to_split The id of the region in the plan object we're attempting to split
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details No modifications made
 *
 * @return <True, information on the edge cut> if two valid regions were
 * successfully split, false otherwise
 *
 */

/*
 * Pick a valid tree edges to split with probability ∝ exp(-alpha*larger abs dev)
 *
 *
 * Returns a valid tree edge to split with probability proporitional to
 * exp(-alpha*larger dev) where larger abs dev is the bigger absolute deviation
 * from the target of the two regions induced by the cut. If successful returns
 * information on the edge and region sizes associated with the cut.
 *
 * Note even if a successful cut is found it does not
 * update the plan or the tree.
 *
 *
 * It will only attempt to create regions where the size is between
 * min_potential_d and max_potential_d (inclusive). So the one district
 * split case is `min_potential_d=max_potential_d=1`.
 *
 * Valid edge here is defined as an edge and region sizes such that the
 * two induced regions both fall within the population bounds.
 *
 *
 * @param root The root vertex of the spanning tree
 * @param pop_below The population corresponding to cutting below each vertex.
 * So `pop_below[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
 * means the vertex is the root or it is not in the tree.
 * @param alpha Used in the exp() term. A larger alpha puts more weight on smaller
 * deviations and smaller makes the weight closer to uniform.
 * @param min_potential_cut_size The smallest potential region size to try for a cut.
 * @param max_potential_cut_size The largest potential region size it will try for a cut.
 * Setting this to 1 will result in only 1 district splits.
 * @param region_ids A vector mapping 0 indexed vertices to their region id number
 * @param region_id_to_split The id of the region in the plan object we're attempting to split
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details No modifications made
 *
 * @return <True, information on the edge cut> if two valid regions were
 * successfully split, false otherwise
 *
 */

/*
 * Pick a valid tree edges to split with probability ∝ exp(-alpha*larger abs dev)
 *
 *
 * Returns a valid tree edge to split with probability proporitional to
 * exp(-alpha*larger dev) where larger abs dev is the bigger absolute deviation
 * from the target of the two regions induced by the cut. If successful returns
 * information on the edge and region sizes associated with the cut.
 *
 * Note even if a successful cut is found it does not
 * update the plan or the tree.
 *
 *
 * It will only attempt to create regions where the size is between
 * min_potential_d and max_potential_d (inclusive). So the one district
 * split case is `min_potential_d=max_potential_d=1`.
 *
 * Valid edge here is defined as an edge and region sizes such that the
 * two induced regions both fall within the population bounds.
 *
 *
 * @param root The root vertex of the spanning tree
 * @param pop_below The population corresponding to cutting below each vertex.
 * So `pop_below[v]` is the population associated with the region made by cutting
 * below the vertex `v`
 * @param tree_vertex_parents The parent of each vertex in the tree. A value of -1
 * means the vertex is the root or it is not in the tree.
 * @param alpha Used in the exp() term. A larger alpha puts more weight on smaller
 * deviations and smaller makes the weight closer to uniform.
 * @param min_potential_cut_size The smallest potential region size to try for a cut.
 * @param max_potential_cut_size The largest potential region size it will try for a cut.
 * Setting this to 1 will result in only 1 district splits.
 * @param region_ids A vector mapping 0 indexed vertices to their region id number
 * @param region_id_to_split The id of the region in the plan object we're attempting to split
 * @param total_region_pop The total population of the region being split
 * @param total_region_size The size of the region being split
 * @param lower Acceptable lower bounds on a valid district's population
 * @param upper Acceptable upper bounds on a valid district's population
 * @param target Ideal population of a valid district. This is what deviance is calculated
 * relative to
 *
 * @details No modifications made
 *
 * @return <True, information on the edge cut> if two valid regions were
 * successfully split, false otherwise
 *
 */
arma::vec compute_expo_prob_weights_on_edges(std::vector<EdgeCut> &valid_edges, double alpha,
                                             double target) {

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++) {
        std::array<double, 2> devs = valid_edges.at(i).compute_abs_pop_deviances(target);
        double bigger_dev = std::max(devs.at(0), devs.at(1));
        unnormalized_wgts(i) = std::exp(-alpha * bigger_dev);
        // Rprintf("Bigger abs dev = %.3f, Computed weight %.3f\n",
        //     bigger_dev, unnormalized_wgts(i));

        // Rprintf("devs are (%.3f,%.3f),  Computed weight %.3f\n",
        //     devs.at(0), devs.at(1), unnormalized_wgts(i));
    }
    // Rprintf("\n\n");

    return unnormalized_wgts;
}

arma::vec compute_expo_prob_weights_on_smaller_dev_edges(std::vector<EdgeCut> &valid_edges,
                                                         double alpha, double target) {

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++) {
        std::array<double, 2> devs = valid_edges.at(i).compute_abs_pop_deviances(target);
        double smaller_dev = std::min(devs.at(0), devs.at(1));
        unnormalized_wgts(i) = std::exp(-alpha * smaller_dev);
        // Rprintf("Bigger abs dev = %.3f, Computed weight %.3f\n",
        //     smaller_dev, unnormalized_wgts(i));

        // Rprintf("devs are (%.6f,%.6f),  Computed weight %.6f\n",
        //     devs.at(0), devs.at(1), unnormalized_wgts(i));
    }
    // Rprintf("\n\n");

    return unnormalized_wgts;
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

// takes an uncut tree and an edge in the tree
// and updates the region ids as if the cut edge was really removed
void assign_region_ids_from_uncut_tree(FlatGraph const &ust, PlanVector &region_ids, int const root,
                                       int const root_region_id, int const cut_vertex_root,
                                       int const cut_vertex_parent, int const cut_region_id,
                                       CircularQueue<std::pair<int, int>> &vertex_queue) {
    // clear the queue
    vertex_queue.clear();

    // update root and add its children to queue
    region_ids[root] = root_region_id;
    for (auto const &child_vertex : ust.neighbors(root)) {
        // check if this edge is one to remove in which case we ignore
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
        for (auto const &child_vertex : ust.neighbors(vertex)) {
            // check if this edge is one to remove in which case we ignore
            if (vertex == cut_vertex_parent && child_vertex == cut_vertex_root)
                continue;

            vertex_queue.push({child_vertex, 0});
        }
    }

    // now we update starting at the cut root vertex
    // update root and add its children to queue
    region_ids[cut_vertex_root] = cut_region_id;
    for (auto const &child_vertex : ust.neighbors(cut_vertex_root)) {
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
        for (auto const &child_vertex : ust.neighbors(vertex)) {
            // since tree is directed we don't need to check anything
            vertex_queue.push({child_vertex, 0});
        }
    }

    return;
}

arma::vec compute_soft_constraint_edge_cut_weights(
    std::vector<EdgeCut> &valid_edges, ScoringFunction const &scoring_function, FlatGraph const &ust,
    int const num_regions, PlanVector &region_ids, RegionSizes &region_sizes,
    IntPlanAttribute &region_pops, int const split_region_id1, int const split_region_id2,
    CircularQueue<std::pair<int, int>> &vertex_queue) {
    arma::vec unnormalized_wgts(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++) {
        // update split info
        region_sizes[split_region_id1] = valid_edges[i].cut_above_region_size;
        region_sizes[split_region_id2] = valid_edges[i].cut_below_region_size;

        region_pops[split_region_id1] = valid_edges[i].cut_above_pop;
        region_pops[split_region_id2] = valid_edges[i].cut_below_pop;

        // update the region ids
        assign_region_ids_from_uncut_tree(ust, region_ids, valid_edges[i].tree_root,
                                          split_region_id1, valid_edges[i].cut_vertex,
                                          valid_edges[i].cut_vertex_parent, split_region_id2,
                                          vertex_queue);

        // get the soft score
        double const score = scoring_function.compute_full_split_plan_soft_score(
            num_regions, region_ids, region_sizes, region_pops, split_region_id1,
            split_region_id2);

        // REprintf("Soft score %f, unnormed weight %.20f\n", score, std::exp(-score));

        unnormalized_wgts[i] = std::exp(-score);
    }

    return unnormalized_wgts;
}

std::vector<EdgeCut> get_all_valid_edges_in_directed_tree(
    const FlatGraph &a_ust, const int root, const std::vector<unsigned int> &pop, TreePopStack &stack,
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

    // Start by adding all the roots children to the stack
    for (auto const &root_children : a_ust.neighbors(root)) {
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

            // Push the vertex back onto the stack as "revisiting"
            stack.push({vtx, parent, true});

            // Push unvisited child vertices onto the stack to get pop below
            for (const auto &child_vtx : a_ust.neighbors(vtx)) {
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
            for (const auto &child : a_ust.neighbors(vtx)) {
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
                total_region_pop - pops_below_vertex[vtx], lower, target, upper,
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
                total_region_pop - pops_below_vertex[vtx], lower, target, upper,
                smaller_cut_sizes_to_try);
        }
    }

}

// finds all valid edges if you joined the two trees
// with the edge (region1_root, region2_root)
// THIS INCLUDES (region1_root, region2_root) as an edge!!
std::vector<EdgeCut> get_valid_edges_in_joined_packed_tree(
    MapParams const &map_params, EdgeBitset const &forest_edges, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int region1_root, const int region1_pop, const int region2_root,
    const int region2_pop, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, const int total_merged_region_size) {
    int const total_merged_region_pop = region1_pop + region2_pop;
    // auto func_start = std::chrono::steady_clock::now();
    // reset pops_below_vertex
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    std::fill(no_valid_edges_vertices.begin(), no_valid_edges_vertices.end(), false);

    // create the valid cut list
    std::vector<EdgeCut> edge_across_valid_edge_cuts;

    // find the valid edges in this half of the tree
    get_all_valid_edges_in_undirected_tree(edge_across_valid_edge_cuts,
        map_params.graph_edge_index, forest_edges,
        region1_root, map_params.pop, stack, pops_below_vertex,
        no_valid_edges_vertices, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try, total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target);

    // find the valid edges in this half of the tree
    get_all_valid_edges_in_undirected_tree(edge_across_valid_edge_cuts,
        map_params.graph_edge_index, forest_edges, 
        region2_root, map_params.pop, stack, pops_below_vertex,
        no_valid_edges_vertices, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try, total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target);

    // Now add the joined cut
    // we make region2 the cut vertex and region1 the parent
    
    get_all_valid_edge_cuts_from_edge(edge_across_valid_edge_cuts, 
        region1_root, region2_root, region1_root, total_merged_region_size,
        static_cast<double>(region2_pop), static_cast<double>(region1_pop), map_params.lower,
        map_params.target, map_params.upper, smaller_cut_sizes_to_try);

    if constexpr (FINDING_EDGE_CUTS_VERBOSE) {
        REprintf("Pop below region2_root is %d so above is %d so foound %d\n",
                 pops_below_vertex.at(region2_root),
                 total_merged_region_pop - pops_below_vertex.at(region2_root),
                 (int)edge_across_valid_edge_cuts.size());
    }


    return edge_across_valid_edge_cuts;
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
void get_all_valid_edges_in_undirected_vertex_tree(
    std::vector<EdgeCut> &existing_cuts,
    Tree const &forest_graph, 
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


    // Start by adding all the roots children to the stack
    for (auto const &root_children : forest_graph[root]) {
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

            // Push the vertex back onto the stack as "revisiting"
            stack.push({vtx, parent, true});


            // Push unvisited child vertices onto the stack to get pop below
            for (const auto &child_vtx : forest_graph[vtx]) {
                // if its the parent then skip it
                if (child_vtx == parent)
                    continue;
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
            for (const auto &child : forest_graph[vtx]) {
                // ignore the parent
                if (child == parent)
                    continue;
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
            // if yes,then the function will add them
            get_all_valid_edge_cuts_from_edge(
                existing_cuts,
                root, vtx, parent, total_region_size, pops_below_vertex[vtx],
                total_region_pop - pops_below_vertex[vtx], lower, target, upper,
                smaller_cut_sizes_to_try);
        }
    }

}


// finds all valid edges if you joined the two trees
// with the edge (region1_root, region2_root)
// THIS INCLUDES (region1_root, region2_root) as an edge!!
std::vector<EdgeCut> get_valid_edges_in_joined_vertex_tree(
    MapParams const &map_params, Tree const &forest_graph, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    const int region1_root, const int region1_pop, const int region2_root,
    const int region2_pop, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, const int total_merged_region_size) {
    int const total_merged_region_pop = region1_pop + region2_pop;
    // auto func_start = std::chrono::steady_clock::now();
    // reset pops_below_vertex
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    std::fill(no_valid_edges_vertices.begin(), no_valid_edges_vertices.end(), false);

    // create the valid cut list
    std::vector<EdgeCut> edge_across_valid_edge_cuts;

    // find the valid edges in this half of the tree
    get_all_valid_edges_in_undirected_vertex_tree(edge_across_valid_edge_cuts,
        forest_graph,
        region1_root, map_params.pop, stack, pops_below_vertex,
        no_valid_edges_vertices, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try, total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target);

    // find the valid edges in this half of the tree
    get_all_valid_edges_in_undirected_vertex_tree(edge_across_valid_edge_cuts,
        forest_graph, 
        region2_root, map_params.pop, stack, pops_below_vertex,
        no_valid_edges_vertices, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try, total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target);

    // Now add the joined cut
    // we make region2 the cut vertex and region1 the parent
    
    get_all_valid_edge_cuts_from_edge(edge_across_valid_edge_cuts, 
        region1_root, region2_root, region1_root, total_merged_region_size,
        static_cast<double>(region2_pop), static_cast<double>(region1_pop), map_params.lower,
        map_params.target, map_params.upper, smaller_cut_sizes_to_try);

    if constexpr (FINDING_EDGE_CUTS_VERBOSE) {
        REprintf("Pop below region2_root is %d so above is %d so foound %d\n",
                 pops_below_vertex.at(region2_root),
                 total_merged_region_pop - pops_below_vertex.at(region2_root),
                 (int)edge_across_valid_edge_cuts.size());
    }


    return edge_across_valid_edge_cuts;
}

}



constexpr bool MERGED_TREE_SPLITTING_VERBOSE = false; // Compile-time constant

std::vector<EdgeCut> TreeSplitter::get_all_valid_pop_edge_cuts_in_directed_tree(
    const MapParams &map_params, FlatGraph const &ust, const int root, TreePopStack &stack,
    std::vector<int> &pops_below_vertex, std::vector<bool> &no_valid_edges_vertices,
    int const region_population, int const region_size, const int min_potential_cut_size,
    const int max_potential_cut_size, std::vector<int> const &smaller_cut_sizes_to_try) const {

    // reset pops_below_vertex and valid edges thing
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);
    std::fill(no_valid_edges_vertices.begin(), no_valid_edges_vertices.end(), false);
    std::vector<EdgeCut> valid_edges = get_all_valid_edges_in_directed_tree(
        ust, root, map_params.pop, stack, pops_below_vertex, no_valid_edges_vertices,
        min_potential_cut_size, max_potential_cut_size, smaller_cut_sizes_to_try,
        region_population, region_size, map_params.lower, map_params.upper, map_params.target);

    return valid_edges;
}

std::pair<bool, EdgeCut> TreeSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
    Plan const &plan, int const split_region1, int const split_region2, FlatGraph const &ust,
    const int root, TreePopStack &stack, std::vector<int> &pops_below_vertex,
    std::vector<bool> &no_valid_edges_vertices, int const region_population,
    int const region_size, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, bool save_selection_prob) {
    // get all the valid edges
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, ust, root, stack, pops_below_vertex, no_valid_edges_vertices,
        region_population, region_size, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if (num_valid_edges == 0) {
        return std::make_pair(false, EdgeCut());
    } else { // else have derived class choose according to its rule
        return select_edge_to_cut(scoring_function, ust, rng_state, valid_edges,
                                  save_selection_prob);
    }
}

// returns edge cut and log probability it was chosen
std::pair<bool, EdgeCut>
TreeSplitter::select_edge_to_cut(ScoringFunction const &scoring_function, FlatGraph const &ust,
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
    EdgeBitset const &forest_edges, TreePopStack &stack, std::vector<bool> &visited,
    std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
    Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try) {
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];

    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size + region2_size;

    // Get all the valid edges in the joined tree
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_packed_tree(
        map_params, forest_edges, stack, pops_below_vertex, visited, region1_root,
        region1_population, region2_root, region2_population, min_potential_cut_size,
        max_potential_cut_size, smaller_cut_sizes_to_try, total_merged_region_size);

    // find the index of the actual edge we cut
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(region1_root, region2_root, region1_root, region2_size,
                            region2_population, region1_size, region1_population);

    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Finding Merge prob for (%d, %d) - %zu valid edges!\n", region1_root,
                 region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    if constexpr (perf_config::bounds_checking){
        if (it == valid_edges.end()) {
            std::ostringstream oss;

            oss << "Actual cut edge not found in valid_edges in "
                << "get_log_retroactive_splitting_prob_for_joined_packed_tree.\n";

            oss << "region1_root=" << region1_root << "\n";
            oss << "region2_root=" << region2_root << "\n";
            oss << "region1_id=" << static_cast<int>(plan.region_ids[region1_root]) << "\n";
            oss << "region2_id=" << static_cast<int>(plan.region_ids[region2_root]) << "\n";
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


double TreeSplitter::get_log_retroactive_splitting_prob_for_joined_vertex_tree(
    MapParams const &map_params, ScoringFunction const &scoring_function,
    Tree const &forest_graph, TreePopStack &stack, std::vector<bool> &visited,
    std::vector<int> &pops_below_vertex, const int region1_root, const int region2_root,
    Plan const &plan, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try) {
    const int region1_population = plan.region_pops[plan.region_ids[region1_root]];
    const int region2_population = plan.region_pops[plan.region_ids[region2_root]];

    const int region1_size = plan.region_sizes[plan.region_ids[region1_root]];
    const int region2_size = plan.region_sizes[plan.region_ids[region2_root]];
    int total_merged_region_size = region1_size + region2_size;

    // Get all the valid edges in the joined tree
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_vertex_tree(
        map_params, forest_graph, stack, pops_below_vertex, visited, region1_root,
        region1_population, region2_root, region2_population, min_potential_cut_size,
        max_potential_cut_size, smaller_cut_sizes_to_try, total_merged_region_size);

    // find the index of the actual edge we cut
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(region1_root, region2_root, region1_root, region2_size,
                            region2_population, region1_size, region1_population);

    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Finding Merge prob for (%d, %d) - %zu valid edges!\n", region1_root,
                 region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    if constexpr (perf_config::bounds_checking){
        if (it == valid_edges.end()) {
            std::ostringstream oss;

            oss << "Actual cut edge not found in valid_edges in "
                << "get_log_retroactive_splitting_prob_for_joined_vertex_tree.\n";

            oss << "region1_root=" << region1_root << "\n";
            oss << "region2_root=" << region2_root << "\n";
            oss << "region1_id=" << static_cast<int>(plan.region_ids[region1_root]) << "\n";
            oss << "region2_id=" << static_cast<int>(plan.region_ids[region2_root]) << "\n";
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
NaiveTopKSplitter::select_edge_to_cut(ScoringFunction const &scoring_function, FlatGraph const &ust,
                                      RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
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
    ScoringFunction const &scoring_function, FlatGraph const &ust, RNGState &rng_state,
    std::vector<EdgeCut> &valid_edges, bool save_selection_prob) const {
    int num_valid_edges = static_cast<int>(valid_edges.size());
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
ExpoWeightedSplitter::compute_unnormalized_edge_cut_weight(EdgeCut const &edge_cut) const {
    std::array<double, 2> devs = edge_cut.compute_abs_pop_deviances(target);
    double bigger_dev = std::max(devs.at(0), devs.at(1));
    return std::exp(-alpha * bigger_dev);
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
    ScoringFunction const &scoring_function, FlatGraph const &ust, RNGState &rng_state,
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

// std::pair<bool, EdgeCut> ConstraintSplitter::select_edge_to_cut(
//     ScoringFunction const &scoring_function, Tree const &ust,
//     RNGState &rng_state, std::vector<EdgeCut> &valid_edges,
//     bool save_selection_prob
// ){
//     int num_valid_edges = static_cast<int>(valid_edges.size());
//     // if only 1 edge just return that
//     if(num_valid_edges == 1) return std::make_pair(true, valid_edges[0]);

//     // get the weights
//     arma::vec unnormalized_wgts = compute_soft_constraint_edge_cut_weights(
//         valid_edges, scoring_function, ust,
//         region_ids, region_sizes, region_pops,
//         int const split_region_id1, int const split_region_id2
//     )

//     // select with prob proportional to the weights
//     int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
//     EdgeCut selected_edge_cut = valid_edges.at(idx);
//     // compute selection probability if needed
//     double log_selection_prob = 0.0;
//     if(save_selection_prob){
//         selected_edge_cut.log_prob = std::log(unnormalized_wgts(idx)) -
//         std::log(arma::sum(unnormalized_wgts));
//     }

//     return std::make_pair(true, selected_edge_cut);
// }

std::pair<bool, EdgeCut> ConstraintSplitter::attempt_to_find_edge_to_cut(
    const MapParams &map_params, ScoringFunction const &scoring_function, RNGState &rng_state,
    Plan const &plan, int const split_region1, int const split_region2, FlatGraph const &ust,
    const int root, TreePopStack &stack, std::vector<int> &pops_below_vertex,
    std::vector<bool> &no_valid_edges_vertices, int const region_population,
    int const region_size, const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try, bool save_selection_prob) {
    // get all the valid edges
    std::vector<EdgeCut> valid_edges = get_all_valid_pop_edge_cuts_in_directed_tree(
        map_params, ust, root, stack, pops_below_vertex, no_valid_edges_vertices,
        region_population, region_size, min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if no valid edges immediately return false
    if (num_valid_edges == 0) {
        return std::make_pair(false, EdgeCut());
    } else if (num_valid_edges == 1) {
        return std::make_pair(true, valid_edges[0]);
    }

    // copy over the current plan information
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);

    // get the weights
    arma::vec unnormalized_wgts = compute_soft_constraint_edge_cut_weights(
        valid_edges, scoring_function, ust, plan.num_regions + 1, region_ids, region_sizes,
        region_pops, split_region1, split_region2, vertex_queue);

    // select with prob proportional to the weights
    int idx = rng_state.r_int_unnormalized_wgt(unnormalized_wgts);
    EdgeCut selected_edge_cut = valid_edges.at(idx);
    // compute selection probability if needed
    if (save_selection_prob) {
        selected_edge_cut.log_prob =
            std::log(unnormalized_wgts(idx)) - std::log(arma::sum(unnormalized_wgts));
        // REprintf("Selection prob %f\n", selected_edge_cut.log_prob);
    }

    return std::make_pair(true, selected_edge_cut);
}

// assumes two trees in spanning forest have been joined
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
    EdgeBitset const &forest_edges, TreePopStack &stack, std::vector<bool> &visited,
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
    std::vector<EdgeCut> valid_edges = get_valid_edges_in_joined_packed_tree(
        map_params, forest_edges, stack, pops_below_vertex, visited, region1_root,
        region1_population, region2_root, region2_population, min_potential_cut_size,
        max_potential_cut_size, smaller_cut_sizes_to_try, total_merged_region_size);

    int num_valid_edges = static_cast<int>(valid_edges.size());
    // if only 1 valid edge then its log(1) = 0
    if (num_valid_edges == 1) {
        return 0.0;
    }

    // find the index of the actual edge we cut
    // where we take region2 root as the cut_vertex
    EdgeCut actual_cut_edge(region1_root, region2_root, region1_root, region2_size,
                            region2_population, region1_size, region1_population);

    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Finding Merge prob for (%d, %d) - %zu valid edges!\n", region1_root,
                 region2_root, valid_edges.size());
    }

    // find the index of the edge we actually removed to get these two regions.
    // it should be 0 if pop bounds are tight but this allows it to work even
    // if not.
    auto it = std::find(valid_edges.begin(), valid_edges.end(), actual_cut_edge);

    // copy over the current plan information
    region_ids.copy(plan.region_ids);
    // copy the region sizes vector
    region_sizes.copy(plan.region_sizes);
    // copy population
    region_pops.copy(plan.region_pops);

    int actual_cut_edge_index = std::distance(valid_edges.begin(), it);
    // copy the forest over
    dummy_forest = forest_edges.get_graph_tree(map_params.graph_edge_index);
    // add the actual removed edge back
    dummy_forest[region1_root].push_back(region2_root);
    dummy_forest[region2_root].push_back(region1_root);

    std::vector<long double> unnormed_wgts;
    unnormed_wgts.reserve(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++) {
        // update split info
        region_sizes[region1_id] = valid_edges[i].cut_above_region_size;
        region_sizes[region2_id] = valid_edges[i].cut_below_region_size;

        region_pops[region1_id] = valid_edges[i].cut_above_pop;
        region_pops[region2_id] = valid_edges[i].cut_below_pop;

        // TODO: Reimplement this 
        std::ostringstream oss;
        oss << "Constraint splitter code is under maintenance right now.\n";
        oss << "This path should not be called.\n";
        throw std::runtime_error(oss.str());
        // update the region ids
        // assign_region_ids_from_joined_undirected_tree(
        //     forest_graph, region_ids, valid_edges[i].cut_vertex, region1_id,
        //     valid_edges[i].cut_vertex_parent, region2_id, vertex_queue);

        // get the soft score
        double const score = scoring_function.compute_full_split_plan_soft_score(
            plan.num_regions, region_ids, region_sizes, region_pops, region1_id, region2_id);

        unnormed_wgts.push_back(std::exp(-score));

        // REprintf("Soft score %f, unnormed weight %.30f vs  %.30f \n", score,
        // unnormed_wgts[i], std::exp(-score));
    }

    auto sum = std::accumulate(unnormed_wgts.begin(), unnormed_wgts.end(), 0.0);
    auto log_sum = std::log(sum);
    // REprintf("Sum %.30f\n", sum);

    // compute selection probability if needed
    double log_selection_prob = std::log(unnormed_wgts[actual_cut_edge_index]) - log_sum;
    // REprintf("Actual Cut Edge at Index %d and so prob is %f \n",
    //     actual_cut_edge_index, log_selection_prob);

    if (MERGED_TREE_SPLITTING_VERBOSE) {
        REprintf("Actual Cut Edge at Index %d and so prob is %f \n", actual_cut_edge_index,
                 log_selection_prob);
    }

    return log_selection_prob;

    // get the weights
    arma::vec unnormalized_wgts(unnormed_wgts.size());

    for (size_t i = 0; i < unnormed_wgts.size(); i++) {
        unnormalized_wgts[i] = std::exp(std::log(unnormed_wgts[i]) - log_sum);
        REprintf("%.30f\n", unnormalized_wgts[i]);
    }

    unnormalized_wgts = arma::cumsum(unnormalized_wgts);

    REprintf("Weights are:\n");
    for (auto const v : unnormalized_wgts) {
        REprintf("%.20f, ", v);
    }
    REprintf("\n");
}