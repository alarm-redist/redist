/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2025/1
* Purpose: Functions for finding an edge in a tree to 
remove (splitting the tree)
********************************************************/

#include "tree_splitting.h"





/*
 * Calculate the deviations for cutting at every edge in a spanning tree.
 * and returns them ordered.
//'
//'
//' For each edge it returns the larger of the two deviations associated with 
//' the best region sizes assignment 
 */
std::vector<double> get_ordered_tree_cut_devs(Tree &ust, int root,
                             std::vector<int> const &cut_below_pop, double const target,
                             const arma::subview_col<arma::uword> &region_ids,
                             int const region_id, int const region_size, int const region_pop,
                             int const min_potential_cut_size, int const max_potential_cut_size,
                             std::vector<int> const &smaller_cut_sizes_to_try
                             ) {
    int V = cut_below_pop.size();
    // compile a list of candidate edges to cut
    std::vector<double> devs; 
    devs.reserve(V); // reserve V which is overkill but ok
    // REprintf("Only looking for region id %d!\n", region_id);
    // REprintf("Starting at %d and %2.f and ", region_pop, static_cast<double>(region_pop)*2.0);
    for (int i = 0; i < V; i++) {
        // ignore vertices not in the region
        if (i == root || region_ids(i) != region_id) continue;


        // start at total pop since deviance will never be more than total_pop/2
        double smallest_dev = static_cast<double>(region_pop)*2.0;
        
        int below_pop = cut_below_pop.at(i);
        int above_pop = region_pop - below_pop;

        // if one of the populations is zero give large deviance and continue
        if(below_pop == 0 || above_pop == 0){
            devs.push_back(smallest_dev);
            continue;
        } 


        // for each possible d value get the deviation above and below
        for(auto const cut_region1_size: smaller_cut_sizes_to_try){
            int cut_region2_size = region_size - cut_region1_size;
            double cut_region1_target = target*cut_region1_size; 
            double cut_region2_target = target*cut_region2_size;


            // REprintf("Size 1=%d-Target %.2f, Size 2=%d-Target %.2f", cut_region1_size, cut_region1_target, cut_region2_size, cut_region2_target);
            // find the larger deviation associated with assigning cut_region1_size to cutting below 
            double max_below_dev = std::max(
                std::fabs(below_pop - cut_region1_target) / cut_region1_target,
                std::fabs(above_pop - cut_region2_target) / cut_region2_target
            );
            // find the larger deviation associated with assigning cut_region1_size to cutting above
            double max_above_dev = std::max(
                std::fabs(above_pop - cut_region1_target) / cut_region1_target,
                std::fabs(below_pop - cut_region2_target) / cut_region2_target
            );
            // REprintf("pair dev %.2f, %.2f", max_below_dev, max_above_dev);

            // take this minimum of this d value and all previous ones 
            double min_dev = std::min(
                max_below_dev,
                max_above_dev
            ) ;
            smallest_dev = std::min(smallest_dev, min_dev);
        }
        // REprintf("Max dev is %.3f\n", smallest_dev);
        devs.push_back(smallest_dev);
    }

    std::sort(devs.begin(), devs.end());

    return devs;
}


/*
 * Return a vector of all valid edge cuts for a particular edge and cut region sizes range
 *
 * 
 * Given an edge in a spanning tree and a range of region sizes to consider for
 * the two cut regions this returns a vector of all the the valid edge cuts 
 * (ie an edge and region sizes for the two cuts) that can be made. For strict
 * population bounds this should be only one but if bounds are loose it can 
 * be multiple. An empty vector means there are no valid edge cuts.
 *
 *
 *
 * 
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
 * @details No modifications made
 *
 * @return A vector of EdgeCut objects
 *
 */ 
inline std::vector<EdgeCut> get_all_valid_edge_cuts_from_edge(
    int const root, int const cut_vertex, int const cut_vertex_parent,
    int const total_region_size, 
    double const below_pop, double const above_pop,
    double const lower, double const target, double const upper,
    std::vector<int> const &smaller_cut_sizes_to_try){

        std::vector<EdgeCut> valid_edges;
        // iterate over all possible valid sizes of the smaller region
        for(auto const cut_region1_size: smaller_cut_sizes_to_try){
            int cut_region2_size = total_region_size - cut_region1_size;
            // Get the bounds for region 1
            double cut_region1_lb = lower*cut_region1_size; 
            double cut_region1_ub = upper*cut_region1_size;
            // Get the bounds for region 2
            double cut_region2_lb = lower*cut_region2_size; 
            double cut_region2_ub = upper*cut_region2_size;

            if(TREE_SPLITTING_DEBUG_VERBOSE){
            REprintf("\tFor (%d, %d): compare %f vs %f vs %f and %f vs %f vs %f \n", 
                cut_region1_size, cut_region2_size,
                cut_region1_lb, below_pop, cut_region1_ub,
                cut_region2_lb, above_pop, cut_region2_ub);
            }

            // check if assigning potential_region_size to cut below leads to valid region
            bool cut_below_ok = cut_region1_lb <= below_pop
                            && below_pop <= cut_region1_ub
                            && cut_region2_lb <= above_pop
                            && above_pop <= cut_region2_ub;

            if(cut_below_ok){
                // cut region 1 size is cut below 
                int cut_below_region_size = cut_region1_size;
                int cut_above_region_size = cut_region2_size;
                valid_edges.emplace_back(root, cut_vertex, cut_vertex_parent,
                    cut_below_region_size, below_pop,
                    cut_above_region_size, above_pop);
                }

            // if both sizes are the same then results are symetric so ignore this case
            if(cut_region1_size == cut_region2_size) continue;


            if(TREE_SPLITTING_DEBUG_VERBOSE){
            REprintf("\tFor (%d, %d): compare %f vs %f vs %f and %f vs %f vs %f \n", 
                cut_region2_size, cut_region1_size,
                cut_region2_lb, below_pop, cut_region2_ub,
                cut_region1_lb, above_pop, cut_region1_ub);
            }

            // check if assigning potential_region_size to cut below leads to valid region
            bool cut_above_ok = cut_region1_lb <= above_pop
                            && above_pop <= cut_region1_ub
                            && cut_region2_lb <= below_pop
                            && below_pop <= cut_region2_ub;

            if(cut_above_ok){
                // cut region 1 size is cut above 
                int cut_below_region_size = cut_region2_size;
                int cut_above_region_size = cut_region1_size;
                valid_edges.emplace_back(root, cut_vertex, cut_vertex_parent,
                    cut_below_region_size, below_pop,
                    cut_above_region_size, above_pop);
                }
        }

    return valid_edges;
}


/*
 * Return a vector of all valid edge cuts in the tree
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
 * @details No modifications made
 *
 * @return A vector of EdgeCut objects
 *
 */ 
std::vector<EdgeCut> get_all_valid_edges_in_directed_tree( 
                     const Tree &ust, const int root,
                     const std::vector<int> &cut_below_pops,
                     const int min_potential_cut_size, const int max_potential_cut_size,
                     std::vector<int> const &smaller_cut_sizes_to_try,
                      const int total_region_pop, const int total_region_size,
                     const double lower, const double upper, const double target){

    // this is the largest size a region can be
    // If the population above is bigger than this you can terminate the serach
    // since pop above only gets larger as you continue down the tree 
    double biggest_upper_bound = upper * max_potential_cut_size;

    // this is the smallest size a region can be 
    // If the pop below is below this then you can terminate the search since 
    // pop below only gets smaller as you continue along the tree 
    double smallest_lower_bound = lower * min_potential_cut_size;

    if(TREE_SPLITTING_DEBUG_VERBOSE){
    REprintf("Region pop=%d, Region size=%d\n",total_region_pop, total_region_size);
    REprintf("lower=%.4f, big upper=%.4f, min_mul = %d, max_mul = %d \n", 
        smallest_lower_bound, biggest_upper_bound, 
        min_potential_cut_size, max_potential_cut_size);
    }

    // vector of valid edge cuts
    std::vector<EdgeCut> valid_edges;

    int V = static_cast<int>(ust.size());
    valid_edges.reserve(std::ceil(V * .01)+1); // heuristic to reserve, no rhyme or reason

    // create a queue of vertices to visit 
    std::queue<std::array<int,2>> vertex_queue;

    // add the children of the root to the queue
    for(const auto &child : ust.at(root)){
        vertex_queue.push({child, root});
    }

    // keep running until the queue isn't empty 
    while(!vertex_queue.empty()){
        // get the vertex at the head of the queue and remove it
        std::array<int,2> queue_element = vertex_queue.front();
        vertex_queue.pop();
        int cut_vertex = queue_element.at(0);
        // get the population below and above it
        double pop_below = static_cast<double>(cut_below_pops.at(cut_vertex));
        double pop_above = static_cast<double>(total_region_pop) - pop_below;

        // check if we should terminate going down this branch. We stop if 
        // - pop_below is less than lower since pop_below only gets smaller
        // - pop_above is bigger than biggest_upper_bound since pop_above
        //   only gets bigger 
        // these conditions should be symmetric but too lazy to prove it now
        if(pop_below < smallest_lower_bound || pop_above > biggest_upper_bound){
            if(TREE_SPLITTING_DEBUG_VERBOSE){
            REprintf("Terminated!! v=%d: Pop above = %f, Pop below = %f, parent = %d\n", 
            cut_vertex, pop_above, pop_below, queue_element.at(1));
            }
            continue;
        }else if(pop_below <= biggest_upper_bound && pop_above >= smallest_lower_bound){
            // At this point its only possible to split if 
            //  - pop_above is >= the smallest lower bound or 
            //  - pop_below is <= than biggest_upper_bound
            if(TREE_SPLITTING_DEBUG_VERBOSE){
            REprintf("Trying: v=%d: Pop above = %f, Pop below = %f, parent = %d\n", 
                cut_vertex, pop_above, pop_below, queue_element.at(1));
            }
            // get vertex parent 
            int cut_vertex_parent = queue_element.at(1);

            // See if any valid edge cuts can be made with this edge 
            std::vector<EdgeCut> new_valid_edges = get_all_valid_edge_cuts_from_edge(
                root, cut_vertex, cut_vertex_parent,
                total_region_size,
                pop_below, pop_above,
                lower, target, upper,
                smaller_cut_sizes_to_try);

            // if yes then add them
            if(new_valid_edges.size() > 0){
                // REprintf("Added v=%d: Pop above = %f, Pop below = %f, parent = %d\n", 
                // cut_vertex, pop_above, pop_below, cut_vertex_parent);
                valid_edges.insert(
                    valid_edges.end(),
                    new_valid_edges.begin(),
                    new_valid_edges.end()
                    );
            }
        }else{
            if(TREE_SPLITTING_DEBUG_VERBOSE){
            REprintf("Skipping: v=%d: Pop above = %f, Pop below = %f, parent = %d\n", 
                cut_vertex, pop_above, pop_below, queue_element.at(1));
            }
        }

        // now iterate over the nodes children 
        for(const auto &child_vertex : ust.at(cut_vertex)){
            // since directed we don't need to worry about visiting twice
            // so just add this vertices childrend to the queue 
            vertex_queue.push({child_vertex, cut_vertex});
        }
    }
    
    return valid_edges;
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
arma::vec compute_expo_prob_weights_on_edges(
        std::vector<EdgeCut> valid_edges, double alpha, double target){

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++)
    {
        std::array<double, 2> devs = valid_edges.at(i).compute_abs_pop_deviances(target);
        double bigger_dev = std::max(devs.at(0), devs.at(1));
        unnormalized_wgts(i) = std::exp(-alpha*bigger_dev);
        // Rprintf("Bigger abs dev = %.3f, Computed weight %.3f\n", 
        //     bigger_dev, unnormalized_wgts(i));

        // Rprintf("devs are (%.3f,%.3f),  Computed weight %.3f\n", 
        //     devs.at(0), devs.at(1), unnormalized_wgts(i));
    }
    // Rprintf("\n\n");
    

    return unnormalized_wgts;

}



arma::vec compute_expo_prob_weights_on_smaller_dev_edges(
        std::vector<EdgeCut> valid_edges, double alpha, double target){

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    for (size_t i = 0; i < valid_edges.size(); i++)
    {
        std::array<double, 2> devs = valid_edges.at(i).compute_abs_pop_deviances(target);
        double smaller_dev = std::min(devs.at(0), devs.at(1));
        unnormalized_wgts(i) = std::exp(-alpha*smaller_dev);
        // Rprintf("Bigger abs dev = %.3f, Computed weight %.3f\n", 
        //     smaller_dev, unnormalized_wgts(i));

        // Rprintf("devs are (%.6f,%.6f),  Computed weight %.6f\n", 
        //     devs.at(0), devs.at(1), unnormalized_wgts(i));
    }
    // Rprintf("\n\n");
    

    return unnormalized_wgts;

}


arma::vec compute_almost_best_weights_on_smaller_dev_edges(
        std::vector<EdgeCut> valid_edges, double epsilon, double target){

    // get the weights vector
    arma::vec unnormalized_wgts(valid_edges.size());

    // find the maximum value 
    double global_min = 42.0;

    for (size_t i = 0; i < valid_edges.size(); i++)
    {
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

    for (size_t i = 0; i < valid_edges.size(); i++){
        // make 1 if eqaul to the max, epsilon otherwise
        // REprintf("Set Weight %d, dev %f to %f \n", 
        //     (int) i, unnormalized_wgts(i), 
        //     (unnormalized_wgts(i) == global_min) ? 1.0 : epsilon);
        unnormalized_wgts(i) = (unnormalized_wgts(i) == global_min) ? 1.0 : epsilon;
    }
    

    return unnormalized_wgts;

}




// finds all valid edges if you joined the two trees
// with the edge (region1_root, region2_root)
// THIS INCLUDES (region1_root, region2_root) as an edge!!
std::vector<EdgeCut> get_valid_edges_in_joined_tree(
    MapParams const &map_params,
    Graph const &forest_graph, Tree &ust,
    std::vector<bool> &visited, std::vector<int> &pops_below_vertex,
    const int region1_id, const int region1_root,
    const int region2_id, const int region2_root,
    const int min_potential_cut_size, const int max_potential_cut_size,
    std::vector<int> const &smaller_cut_sizes_to_try,
    const int total_merged_region_pop, const int total_merged_region_size
){
    // clear the tree
    clear_tree(ust);
    // clear the pop below and visited 
    std::fill(visited.begin(), visited.end(), false);
    std::fill(pops_below_vertex.begin(), pops_below_vertex.end(), 0);


    // auto t1 = std::chrono::high_resolution_clock::now();
    // build the tree starting from root 1
    build_directed_tree_and_get_pops_below(
        forest_graph, 
        ust, region1_root, visited, 
        map_params.pop, pops_below_vertex);
    // auto t2 = std::chrono::high_resolution_clock::now();
    //     /* Getting number of milliseconds as a double. */
    // std::chrono::duration<double, std::milli> ms_double = t2 - t1; 
    // Rcout << "Building Tree " << ms_double.count() << " ms\n";

    // auto t1f = std::chrono::high_resolution_clock::now();
    // find the valid edges in this half of the tree 
    std::vector<EdgeCut> valid_tree1_edges = get_all_valid_edges_in_directed_tree(
        ust, region1_root, 
        pops_below_vertex,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target
    );
    // auto t2f = std::chrono::high_resolution_clock::now();
    //     /* Getting number of milliseconds as a double. */
    // ms_double = t2f - t1f; 
    // Rcout << "Getting Edges in Tree " << ms_double.count() << " ms\n";

    // build the tree starting from root 2
    build_directed_tree_and_get_pops_below(
        forest_graph, 
        ust, region2_root, visited, 
        map_params.pop, pops_below_vertex);
    // find the valid edges in this half of the tree 
    std::vector<EdgeCut> valid_tree2_edges = get_all_valid_edges_in_directed_tree(
        ust, region2_root, 
        pops_below_vertex,
        min_potential_cut_size, max_potential_cut_size,
        smaller_cut_sizes_to_try,
        total_merged_region_pop, total_merged_region_size,
        map_params.lower, map_params.upper, map_params.target
    );

    // Now add the joined cut
    // we make region2 the cut vertex and region1 the parent
    std::vector<EdgeCut> edge_across_valid_edge_cuts = get_all_valid_edge_cuts_from_edge(
                region1_root, region2_root, region1_root,
                total_merged_region_size,
                static_cast<double>(pops_below_vertex.at(region2_root)), 
                static_cast<double>(total_merged_region_pop - pops_below_vertex.at(region2_root)),
                map_params.lower, map_params.target, map_params.upper,
                smaller_cut_sizes_to_try);

    // REprintf("Pop below region2_root is %d so above is %d so foound %d\n",
    //     pops_below_vertex.at(region2_root), 
    //     total_merged_region_pop -pops_below_vertex.at(region2_root),
    //     (int) edge_across_valid_edge_cuts.size());

    // now add the edges from the two trees
    edge_across_valid_edge_cuts.insert(
        edge_across_valid_edge_cuts.end(),
        valid_tree1_edges.begin(),
        valid_tree1_edges.end()
    );

    edge_across_valid_edge_cuts.insert(
        edge_across_valid_edge_cuts.end(),
        valid_tree2_edges.begin(),
        valid_tree2_edges.end()
    );

    return edge_across_valid_edge_cuts;
}



