#include "gsmc.h"


//' Attempts to cut one region into two from spanning tree
//'
//' Takes a spanning tree `ust` drawn on a specific region and attempts to cut
//' it to produce two new regions using the generalized splitting procedure
//' outlined <PAPER HERE>. This function is based on `cut_districts` in `smc.cpp`
//'
//' @title Attempt Generalized Region split
//'
//' @param ust A directed spanning tree passed by reference
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param root The root vertex of the spanning tree
//' @param pop A vector of the population associated with each vertex in `g`
//' @param plan A plan object
//' @param region_to_split The label of the region in the plan object we're attempting to split
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param new_region_labels A vector that will be updated by reference to contain the names of
//' the two new split regions if function is successful.
//'
//' @details Modifications
//'    - If two new valid regions are split then the plan object is updated accordingly
//'    - If two new valid regions are split then the new_region_labels is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'
//' @return True if two valid regions were split off false otherwise
//'
bool cut_regions(Tree &ust, int k_param, int root,
                     const uvec &pop,
                     Plan &plan, const std::string region_to_split,
                     double lower, double upper, double target,
                     std::vector<std::string> &new_region_labels){
    // Get population of region being split
    double total_pop = plan.r_to_pop_map.at(region_to_split);
    // Get the number of final districts in region being split
    int num_final_districts = plan.r_to_d_map.at(region_to_split);

    // Rcout << "For " << region_to_split << " Total pop is " << total_pop << " and d_nk is " << num_final_districts << "\n";

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(plan.V, 0);
    std::vector<int> parent(plan.V);
    parent[root] = -1;
    tree_pop(ust, root, pop, pop_below, parent);

    // compile a list of: things for each edge in tree
    std::vector<int> candidates; // candidate edges to cut,
    std::vector<double> deviances; // how far from target pop.
    std::vector<bool> is_ok; // whether they meet constraints
    std::vector<int> new_d_val; // the new d_n,k value

    // Reserve at least as much space for top k_param of them
    candidates.reserve(k_param);
    deviances.reserve(k_param);
    is_ok.reserve(k_param);
    new_d_val.reserve(k_param);

    if(plan.region_labels.at(root) != region_to_split){
        Rcout << "Root vertex is not in region to split!";
    }

    // Now loop over all valid edges to cut
    for (int i = 1; i <= plan.V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (plan.region_labels.at(i-1) != region_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        double below = pop_below.at(i - 1);
        // Get the population of the other side
        double above = total_pop - below;

        // vectors to keep track of value for each d_{n,k}
        std::vector<double> devs(num_final_districts-1);
        std::vector<bool> dev2_bigger(num_final_districts-1);
        std::vector<bool> are_ok(num_final_districts-1);


        // Now try each potential d_nk value from 1 up to d_n-1,k -1
        // remember to correct for 0 indexing
        for(int potential_d = 1; potential_d < num_final_districts; potential_d++){

            double dev1 = std::fabs(below - target * potential_d);
            double dev2 = std::fabs(above - target * potential_d);

            /*
             Rcout << "\nit is " << dev1 << " and " << dev2 << " for " << potential_d << "\n";
             Rcout << std::printf(
             "We are on d = %d and upper is %f while lower is %f \n",
             potential_d, above, below);
             */

            // If dev1 is smaller then we assign d_nk to below
            if (dev1 < dev2) {
                devs[potential_d-1] = dev1;
                dev2_bigger[potential_d-1] = true;
                // check in bounds
                are_ok[potential_d-1] =    lower * potential_d <= below
                                        && below <= upper * potential_d
                                        && lower * (num_final_districts - potential_d) <= above
                                        && above <= upper * (num_final_districts - potential_d);
            } else { // Else if dev2 is smaller we assign d_nk to above
                devs[potential_d-1] = dev2;
                dev2_bigger[potential_d-1] = false;
                are_ok[potential_d-1] =    lower * potential_d <= above
                                        && above <= upper * potential_d
                                        && lower * (num_final_districts - potential_d) <= below
                                        && below <= upper * (num_final_districts - potential_d);
            }

        }

        /*
         Rcout << "The full vector is ";
         for (auto i: devs)
         Rcout << i << ' ';
         Rcout << "\n";
         */

        // Now find the value of d_{n,k} that has the smallest deviation
        std::vector<double>::iterator result = std::min_element(devs.begin(), devs.end());
        int best_potential_d = std::distance(devs.begin(), result);

        /*
         Rcout << "min element has value " << *result << " and index ["
               << best_potential_d << "]\n";
         */

        if(dev2_bigger[best_potential_d]){
            candidates.push_back(i);
        } else{
            candidates.push_back(-i);
        }

        deviances.push_back(devs[best_potential_d]);
        is_ok.push_back(are_ok[best_potential_d]);
        new_d_val.push_back(best_potential_d+1);

    }


    // if less than k_param candidates immediately reject
    if((int) candidates.size() < k_param){
        return false;
    }

    int idx = r_int(k_param);
    idx = select_k(deviances, idx + 1);
    int cut_at = std::fabs(candidates[idx]) - 1;


    // Rcout <<  deviances[idx] << "\n";
    // Rcout << is_ok[idx] << " howdy \n";
    // Rcout << "top k element has value " << deviances[idx] << "\n";


    // reject sample if not ok
    if (!is_ok[idx]){
        return false;
    }

    /*
    Rcout << "Testing this " << region_to_split + ".R2" << " \n";
    Rcout << "We have Final d =" << new_d_val[idx] << " and " << num_final_districts - new_d_val[idx] << "\n";
     */

    // find index of node to cut at
    std::vector<int> *siblings = &ust[parent[cut_at]];
    int length = siblings->size();
    int j;
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_at) break;
    }

    siblings->erase(siblings->begin()+j); // remove edge
    parent[cut_at] = -1;


    // update plan with new regions
    plan.num_regions++; // increase region count by 1
    plan.num_multidistricts--; // Decrease by one to avoid double counting later

    // Create info for two new districts
    int new_region1_d = new_d_val[idx];
    int new_region2_d = num_final_districts - new_region1_d;

    std::string new_region_label1;
    std::string new_region_label2;

    // Set label and count depending on if district or multi district
    if(new_region1_d == 1){
        plan.num_districts++;
        new_region_label1 = region_to_split + "." + std::to_string(plan.num_districts);
    }else{
        new_region_label1 = region_to_split + ".R1";
        plan.num_multidistricts++;
    }

    // Now do it for second region
    if(new_region2_d == 1){
        plan.num_districts++;
        new_region_label2 = region_to_split + "." + std::to_string(plan.num_districts);
    }else{
        new_region_label2 = region_to_split + ".R2";
        plan.num_multidistricts++;
    }

    // Check whether the new regions are districts or multidistricts


    // Vertex to start traversing tree for updating later
    int tree_vertex1;
    int tree_vertex2;

    // population and d_nk of new regions
    double new_region1_pop;
    double new_region2_pop;


    if (candidates[idx] > 0) { // Means cut below so first vertex is cut_at
        tree_vertex1 = cut_at;
        tree_vertex2 = root;

        // Set the new populations
        new_region1_pop = pop_below.at(cut_at);
        new_region2_pop = total_pop - new_region1_pop;

    //     return pop_below.at(cut_at);
    } else { // Means cut above so first vertex is root
        tree_vertex1 = root;
        tree_vertex2 = cut_at;

        // Set the new populations
        new_region2_pop = pop_below.at(cut_at);
        new_region1_pop = total_pop - new_region2_pop;

    }


    // update the function parameter with names of new regions
    if(new_region_labels.size() != 2){
        Rcout << "For some reason the new_region_labels vector in cut_regions is not size 2!\n";
    }
    new_region_labels[0] = new_region_label1;
    new_region_labels[1] = new_region_label2;

    // Get the regions current integer id
    int old_region_num_id = plan.str_label_to_num_id_map[region_to_split];
    // make the first new region have the same integer id
    int new_region_num_id1 = old_region_num_id;
    // Second new region has id of the new number of regions minus 1
    int new_region_num_id2 = plan.num_regions - 1;

    // Now update the two cut portions
    assign_region(ust, plan, tree_vertex1, new_region_label1, new_region_num_id1, new_region1_pop, new_region1_d);
    assign_region(ust, plan, tree_vertex2, new_region_label2, new_region_num_id2, new_region2_pop, new_region2_d);

    // Now update the map objects in the plan
    plan.r_to_d_map.erase(region_to_split);
    plan.r_to_pop_map.erase(region_to_split);
    plan.str_label_to_num_id_map.erase(region_to_split);
    plan.num_id_to_str_label_map.erase(old_region_num_id);

    // Add the new region 1
    plan.r_to_d_map[new_region_label1] = new_region1_d;
    plan.r_to_pop_map[new_region_label1] = new_region1_pop;
    plan.str_label_to_num_id_map[new_region_label1] = new_region_num_id1;
    plan.num_id_to_str_label_map[new_region_num_id1] = new_region_label1;

    // Add the new region 2
    plan.r_to_d_map[new_region_label2] = new_region2_d;
    plan.r_to_pop_map[new_region_label2] = new_region2_pop;
    plan.str_label_to_num_id_map[new_region_label2] = new_region_num_id2;
    plan.num_id_to_str_label_map[new_region_num_id2] = new_region_label2;

    return true;

};



//' Attempts to split a multi-district within a plan into two new regions with
//' valid population bounds
//'
//' Given a plan this attempts to split a multi-district in it into two new
//' regions where both regions have valid population bounds. Does this by
//' drawing a spanning tree uniformly at random then calling `cut_regions` on
//' that. If the split it successful it returns true and modifies `plan` and
//' `new_region_labels` accordingly. This is based on the `split_map` function
//' in smc.cpp
//'
//'
//' @title Attempt Generalized Region split of a multi-district within a plan
//'
//' @param g A graph (adjacency list) passed by reference
//' @param ust A directed tree object (this will be cleared in the function so
//' whatever it was before doesn't matter)
//' @param counties Vector of county labels of each vertex in `g`
//' @param cg multigraph object (not sure why this is needed)
//' @param plan A plan object
//' @param region_to_split The label of the region in the plan object we're attempting to split
//' @param new_region_labels A vector that will be updated by reference to contain the names of
//' the two new split regions if function is successful.
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' Target population (probably Total population of map/Num districts)
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//'
//' @details Modifications
//'    - If two new valid regions are split then the plan object is updated accordingly
//'    - If two new valid regions are split then the new_region_labels is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'
//' @return True if two valid regions were split off false otherwise
//'
bool attempt_region_split(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                 Plan &plan, const std::string region_to_split,
                 std::vector<std::string> &new_region_labels,
                 std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                 double &lower, double upper, double target, int k_param) {

    int V = g.size();

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_labels[i] != region_to_split;
    }

    // Get a uniform spanning tree drawn on that region
    int root;
    clear_tree(ust);
    // Get a tree
    int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
    // Return unsuccessful if tree not drawn
    if (result != 0) return false;

    // Now try to cut the tree and return that result
    return cut_regions(ust, k_param, root, pop,
                       plan, region_to_split,
                       lower, upper, target,
                       new_region_labels);

}





/*
 * Split off a piece from each map in `districts`,
 * keeping deviation between `lower` and `upper`
 */

// This should just update ancestor, lag, and give log prob reason was picked.
// Everything else can happen outside.
// Actually we compute incremental weight in the function.


//' Splits a multidistrict in all of the plans
//'
//' Using the procedure outlined in <PAPER HERE> this function attempts to split
//' a multidistrict in a previous steps plan until M successful splits have been made. This
//' is based on the `split_maps` function in smc.cpp
//'
//' @title Split all the maps
//'
//' @param g A graph (adjacency list) passed by reference
//' @param counties Vector of county labels of each vertex in `g`
//' @param cg County level multigraph
//' @param pop A vector of the population associated with each vertex in `g`
//' @param old_plans_vec A vector of plans from the previous step
//' @param new_plans_vec A vector which will be filled with plans that had a
//' multidistrict split to make them
//' @param original_ancestor_vec A vector used to track which original ancestor
//' the new plans descended from. The value  of `original_ancestor_vec[i]`
//' is the index of the original ancestor the new plan `new_plans_vec[i]` is
//' descended from.
//' @param parent_vec A vector used to track the index of the previous plan
//' sampled that was successfully split. The value of `parent_vec[i]` is the
//' index of the old plan from which the new plan `new_plans_vec[i]` was
//' successfully split from. In other words `new_plans_vec[i]` is equal to
//' `attempt_region_split(old_plans_vec[parent_vec[i]], ...)`
//' @param prev_ancestor_vec A vector used to track the index of the original
//' ancestor of the previous plans. The value of `prev_ancestor_vec[i]` is the
//' index of the original ancestor of `old_plans_vec[i]`
//' @param new_log_incremental_weights A vector of the new log incremental weights
//' computed for the new plans. The value of `new_log_incremental_weights[i]` is
//' the log incremental weight for `new_plans_vec[i]`
//' @param unnormalized_sampling_weights A vector of weights used to sample indices
//' of the `old_plans_vec`. The value of `unnormalized_sampling_weights[i]` is
//' the unnormalized probability that index i is selected
//' @param normalized_weights_to_fill_in A vector which will be filled with the
//' normalized weights the index sampler uses. The value of
//' `normalized_weights_to_fill_in[i]` is the probability that index i is selected
//' @param draw_tries_vec A vector used to keep track of how many plan split
//' attempts were made for index i. The value `draw_tries_vec[i]` represents how
//' many split attempts were made for the i-th new plan (including the successful
//' split). For example, `draw_tries_vec[i] = 1` means that the first split
//' attempt was successful.
//' @param accept_rate The number of accepted splits over the total number of
//' attempted splits. This is equal to `sum(draw_tries_vec)/M`
//' @param n_unique_parent_indices The number of unique parent indices, ie the
//' number of previous plans that had at least one descendant amongst the new
//' plans. This is equal to `unique(parent_vec)`
//' @param n_unique_original_ancestors The number of unique original ancestors,
//' in the new plans. This is equal to `unique(original_ancestor_vec)`
//' @param ancestors Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lags Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
//' WHAT IT IS DOING
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param k_param The top edges to pick parameter for the region splitting
//' algorithm
//' @param pool A threadpool for multithreading
//' @param verbosity A parameter controlling the amount of detail printed out
//' during the algorithms running
//'
//' @details Modifications
//'    - The `new_plans_vec` is updated with all the newly split plans
//'    - The `old_plans_vec` is updated with all the newly split plans as well.
//'    Note that the reason both this and `new_plans_vec` are updated is because
//'    of the nature of the code you need both vectors and so both are passed by
//'    reference to save memory.
//'    - The `original_ancestor_vec` is updated to contain the indices of the
//'    original ancestors of the new plans
//'    - The `parent_vec` is updated to contain the indices of the parents of the
 //'    new plans
//'    - If two new valid regions are split then the new_region_labels is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'    - The `new_log_incremental_weights` is updated to contain the incremental
//'    weights of the new plans NOTE: In the future this will be moved to its
//'    own function
//'    - The `normalized_weights_to_fill_in` is updated to contain the normalized
//'    probabilities the index sampler used. This is only collected for diagnostics
//'    at this point and should just be equal to `unnormalized_sampling_weights`
//'    divided by `sum(unnormalized_sampling_weights)`
//'    - The `draw_tries_vec` is updated to contain the number of tries for each
//'    of the new plans
//'    - The `accept_rate` is updated to contain the average acceptance rate for
//'    this iteration
//'    - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
//'    with the unique number of parents and original ancestors for all the new
//'    plans respectively
//'    - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
//'    I DO NOT KNOW WHAT IT MEANS
//'
//' @return nothing
//'
void generalized_split_maps(
        const Graph &g, const uvec &counties, Multigraph &cg, const uvec &pop,
        std::vector<Plan> &old_plans_vec, std::vector<Plan> &new_plans_vec,
        std::vector<int> &original_ancestor_vec,
        std::vector<int> &parent_vec,
        const std::vector<int> &prev_ancestor_vec,
        std::vector<double> &new_log_incremental_weights,
        const std::vector<double> &unnormalized_sampling_weights,
        std::vector<double> &normalized_weights_to_fill_in,
        std::vector<int> &draw_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        int &n_unique_original_ancestors,
        umat &ancestors, const std::vector<int> &lags,
        double lower, double upper, double target,
        int k_param,
        RcppThread::ThreadPool &pool,
        int verbosity
                ) {
    // important constants
    const int V = g.size();
    const int M = old_plans_vec.size();


    uvec iters(M, fill::zeros); // how many actual iterations, (used to compute acceptance rate)
    uvec original_ancestor_uniques(M); // used to compute unique original ancestors
    uvec parent_index_uniques(M); // used to compute unique parent indicies

    // PREVIOUS SMC CODE I DONT KNOW WHAT IT DOES
    const int dist_ctr = old_plans_vec.at(0).num_regions;
    const int n_lags = lags.size();
    umat ancestors_new(M, n_lags); // lags/ancestor thing



    // Create the obj which will sample from the index with probability
    // proportional to the weights
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> index_sampler(
            unnormalized_sampling_weights.begin(),
            unnormalized_sampling_weights.end()
        );

    // extract and record the normalized weights the sampler is using
    std::vector<double> p = index_sampler.probabilities();
    int nw_index = 0;
    bool print_weights = M < 12 && verbosity > 1;
    if(print_weights){
    Rprintf("Unnormalized weights are: ");
    for (auto w : unnormalized_sampling_weights){
        Rprintf("%.4f, ", w);
    }

    Rprintf("\n");
    Rprintf("Normalized weights are: ");
    }
    for (auto prob : p){
        if(print_weights) Rprintf("%.4f, ", prob);
        normalized_weights_to_fill_in.at(nw_index) = prob;
        nw_index++;
    }
    if(print_weights) Rprintf("\n");

    // Because of multithreading we have to add specific checks for if the user
    // wants to quit the program
    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50; // check for interrupts every _ iterations


    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        int reject_ct = 0;
        bool ok = false;
        int idx;
        std::string region_to_split;
        std::vector<std::string> new_region_labels(2, "INIT");

        Tree ust = init_tree(V);
        std::vector<bool> visited(V);
        std::vector<bool> ignore(V);
        while (!ok) {
            // increase the iters count by one
            iters[i]++;
            // use weights to sample previous plan
            idx = index_sampler(gen);
            Plan proposed_new_plan = old_plans_vec[idx];

            // pick a region to try to split
            choose_multidistrict_to_split(
                old_plans_vec[idx], region_to_split);

            // Now try to split that region
            ok = attempt_region_split(g, ust, counties, cg,
                                      proposed_new_plan, region_to_split,
                                      new_region_labels,
                                      visited, ignore, pop,
                                      lower, upper, target, k_param);

            // bad sample; try again
            if (!ok) {
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);
                continue;
            }

            // else update the new plan and leave the while loop
            new_plans_vec[i] = proposed_new_plan;

        }

        // compute the log incremental weight for the new plan
        new_log_incremental_weights[i] = compute_log_incremental_weight(g, new_plans_vec[i]);
        // Record how many tries needed to create i-th new plan
        draw_tries_vec[i] = static_cast<int>(iters[i]);
        // Make the new plans original ancestor the same as its parent
        original_ancestor_uniques[i] = prev_ancestor_vec[idx];
        // record index of new plan's parent
        parent_index_uniques[i] = idx;
        // clear the spanning tree
        clear_tree(ust);

        // ORIGINAL SMC CODE I DONT KNOW WHAT THIS DOES
        // save ancestors/lags
        for (int j = 0; j < n_lags; j++) {
            if (dist_ctr <= lags[j]) {
                ancestors_new(i, j) = i;
            } else {
                ancestors_new(i, j) = ancestors(idx, j);
            }
        }


        // update this particles ancestor to be the ancestor of its previous one
        parent_vec[i] = idx;
        original_ancestor_vec[i] = prev_ancestor_vec[idx];

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();



    // now replace the old plans with the new ones
    for(int i=0; i < M; i++){
        old_plans_vec[i] = new_plans_vec[i];
    }


    // now compute acceptance rate and unique parents and original ancestors
    accept_rate = M / (1.0 * sum(iters));
    n_unique_original_ancestors = ((uvec) find_unique(original_ancestor_uniques)).n_elem;
    n_unique_parent_indices = ((uvec) find_unique(parent_index_uniques)).n_elem;
    if (verbosity >= 3) {
        Rprintf("%.2f acceptance rate, %d unique parent indices sampled, and %d unique original ancestors!\n",
                100.0 * accept_rate, (int) n_unique_parent_indices , (int) n_unique_original_ancestors);
    }

    if(print_weights){
        Rprintf("Log Incremental weights are: ");
        for (auto w : new_log_incremental_weights){
            Rprintf("%.4f, ", w);
        }
        Rprintf("\n");
    }

    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;

}



//' Uses gsmc method to generate a sample of `M` plans in `c++`
//'
//' Using the procedure outlined in <PAPER HERE> this function uses Sequential
//' Monte Carlo (SMC) methods to generate a sample of `M` plans
//'
//' @title Run redist gsmc
//'
//' @param N The number of districts the final plans will have
//' @param adj_list A 0-indexed adjacency list representing the undirected graph
//' which represents the underlying map the plans are to be drawn on
//' @param counties Vector of county labels of each vertex in `g`
//' @param pop A vector of the population associated with each vertex in `g`
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param M The number of plans (samples) to draw
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param control Named list of additional parameters.
//' @param num_threads The number of threads the threadpool should use
//' @param verbosity What level of detail to print out while the algorithm is
//' running <ADD OPTIONS>
//' @export
List gsmc_plans(
        int N, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        double target, double lower, double upper,
        int M, int k_param, // M is Number of particles aka number of different plans
        List control,
        int num_threads, int verbosity, bool diagnostic_mode){

    // set number of threads
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 1) num_threads = 0;

    // lags thing (copied from original smc code, don't understand what its doing)
    std::vector<int> lags = as<std::vector<int>>(control["lags"]);
    umat ancestors(M, lags.size(), fill::zeros);

    // Create map level graph and county level multigraph
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);

    int V = g.size();
    double total_pop = sum(pop);

    // Loading Info
    if (verbosity >= 1) {
        Rcout.imbue(std::locale(""));
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "GENERALIZED SEQUENTIAL MONTE CARLO\n";
        Rcout << "Sampling " << M << " " << V << "-unit ";
        Rcout << "maps with " << N << " districts and population between "
              << lower << " and " << upper << " using " << num_threads << " threads.\n";
        if (cg.size() > 1){
            Rcout << "Ensuring no more than " << N - 1 << " splits of the "
                  << cg.size() << " administrative units.\n";
        }
    }



    std::vector<Plan> plans_vec(M, Plan(V, N, total_pop));
    std::vector<Plan> new_plans_vec(M, Plan(V, N, total_pop)); // New plans


    // Define output variables that must always be created

    // This is N-1 by M where [i][j] is the index of the parent of particle j on step i
    // ie the index of the previous plan that was sampled and used to create particle j on step i
    std::vector<std::vector<int>> parent_index_mat(N-1, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the index of the original (first) ancestor of particle j on step i
    std::vector<std::vector<int>> original_ancestor_mat(N-1, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the number of tries it took to form particle j on iteration i
    // Inclusive of the final step. ie if succeeds in one try it would be 1
    std::vector<std::vector<int>> draw_tries_mat(N-1, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the log incremental weight of particle j on step i
    std::vector<std::vector<double>> log_incremental_weights_mat(N-1, std::vector<double> (M, -1.0));

    // This is N-1 by M where [i][j] is the normalized weight of particle j on step i
    std::vector<std::vector<double>> normalized_weights_mat(N-1, std::vector<double> (M, -1.0));




    // Tracks the acceptance rate - total number of tries over M - for each round
    std::vector<double> acceptance_rates(N-1, -1.0);

    // Tracks the effective sample size for the weights of each round
    std::vector<double> n_eff(N-1, -1.0);

    // Tracks the number of unique parent vectors sampled to create the next round
    std::vector<int> nunique_parents_vec(N-1, -1);

    // Tracks the number of unique ancestors left at each step
    std::vector<int> nunique_original_ancestors_vec(N-1, -1);


    // Declare variables whose size will depend on whether or not we're in
    // diagnostic mode or not
    std::vector<std::vector<std::vector<int>>> plan_region_ids_mat;
    std::vector<std::vector<std::vector<int>>> plan_d_vals_mat;
    std::vector<std::vector<std::string>> final_plan_region_labels;




    // If diagnostic mode track stuff from every round
    if(diagnostic_mode){
        // Create info tracking we will pass out at the end
        // This is N-1 by M by V where for each n=1,...,N-1 and m=1,...M it maps vertices to integer id
        plan_region_ids_mat.resize(N-1, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));
        // For integers make default -1 for error tracking
        // This is N-1 by M by V where for each n=1,...,N-1 and m=1,...M it maps vertices to region d val
        plan_d_vals_mat.resize(N-1, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));

        //  M by V where for each m=1,...M it maps vertices to region labels in a plan
        final_plan_region_labels.resize(
                M, std::vector<std::string> (V, "MISSING")
        );

    }else{ // else only track for final round
        // Create info tracking we will pass out at the end
        // This is M by V where for each m=1,...M it maps vertices to integer id for plan
        plan_region_ids_mat.resize(1, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));

        // For integers make default -1 for error tracking
        // This is M by V where for each m=1,...M it maps vertices to d value for plan
        plan_d_vals_mat.resize(1, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));
    }





    // Start off all the unnormalized weights at 1
    std::vector<double> unnormalized_sampling_weights(M, 1.0);


    // Create a threadpool

    RcppThread::ThreadPool pool(num_threads);

    std::string bar_fmt = "Split [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} | ETA{cli::pb_eta}";
    RObject bar = cli_progress_bar(N-1, cli_config(false, bar_fmt.c_str()));

    // Now for each run through split the map
    try {
    for(int n=0; n<N-1; n++){
        if(verbosity > 1){
            Rprintf("Iteration %d \n", n+1);
        }


        // For the first iteration we need to pass a special previous ancestor thing
        if(n == 0){
        std::vector<int> dummy_prev_ancestors(M, 1);
        // split the map
        generalized_split_maps(
            g, counties, cg, pop,
            plans_vec, new_plans_vec,
            original_ancestor_mat[n],
            parent_index_mat[n],
            dummy_prev_ancestors,
            log_incremental_weights_mat[n],
            unnormalized_sampling_weights,
            normalized_weights_mat[n],
            draw_tries_mat[n],
            acceptance_rates[n],
            nunique_parents_vec[n],
            nunique_original_ancestors_vec[n],
            ancestors, lags,
            lower, upper, target,
            k_param,
            pool,
            verbosity
        );

            // For the first ancestor one make every ancestor themselves
            std::iota (parent_index_mat[0].begin(), parent_index_mat[0].end(), 0);
            std::iota (original_ancestor_mat[0].begin(), original_ancestor_mat[0].end(), 0);
        }else{
        // split the map and we can use the previous original ancestor matrix row
        generalized_split_maps(
            g, counties, cg, pop,
            plans_vec, new_plans_vec,
            original_ancestor_mat[n],
            parent_index_mat[n],
            original_ancestor_mat[n-1],
            log_incremental_weights_mat[n],
            unnormalized_sampling_weights,
            normalized_weights_mat[n],
            draw_tries_mat[n],
            acceptance_rates[n],
            nunique_parents_vec[n],
            nunique_original_ancestors_vec[n],
            ancestors, lags,
            lower, upper, target,
            k_param,
            pool,
            verbosity
        );
        }

        if (verbosity == 1 && CLI_SHOULD_TICK){
            cli_progress_set(bar, n);
        }
        Rcpp::checkUserInterrupt();



        // Now update the weights, region labels, dval column of the matrix
        for(int j=0; j<M; j++){
            // Make the unnormalized weight just the incremental weight from this round
            // In the future make standalone function to compute these weights
            unnormalized_sampling_weights[j] = std::exp(
                log_incremental_weights_mat.at(n).at(j) // + std::log(unnormalized_sampling_weights[j])
            );

            if(diagnostic_mode && n == N-2){ // record if in diagnostic mode and final step
                plan_region_ids_mat.at(n).at(j) = plans_vec[j].region_num_ids;
                plan_d_vals_mat.at(n).at(j) = plans_vec[j].region_dval;
                final_plan_region_labels.at(j) = plans_vec[j].region_labels;
            }else if(diagnostic_mode){ // record if in diagnostic mode but not final step
                plan_region_ids_mat.at(n).at(j) = plans_vec[j].region_num_ids;
                plan_d_vals_mat.at(n).at(j) = plans_vec[j].region_dval;
            }else if(n == N-2){ // else if not only record final step
                plan_region_ids_mat.at(0).at(j) = plans_vec[j].region_num_ids;
                plan_d_vals_mat.at(0).at(j) = plans_vec[j].region_dval;
            }

        }


        n_eff.at(n) = compute_n_eff(log_incremental_weights_mat[n]);



    }
    } catch (Rcpp::internal::InterruptedException e) {
        cli_progress_done(bar);
        return R_NilValue;
    }

    cli_progress_done(bar);


    // make first number of unique original ancestors just M
    nunique_original_ancestors_vec.at(0) = M;

    // Return results
    List out = List::create(
        _["original_ancestors"] = original_ancestor_mat,
        _["parent_index"] = parent_index_mat,
        _["final_region_labs"] = final_plan_region_labels,
        _["region_ids_mat_list"] = plan_region_ids_mat,
        _["region_dvals_mat_list"] = plan_d_vals_mat,
        _["log_incremental_weights_mat"] = log_incremental_weights_mat,
        _["normalized_weights_mat"] = normalized_weights_mat,
        _["draw_tries_mat"] = draw_tries_mat,
        _["acceptance_rates"] = acceptance_rates,
        _["nunique_parent_indices"] = nunique_parents_vec,
        _["nunique_original_ancestors"] = nunique_original_ancestors_vec,
        _["ancestors"] = ancestors,
        _["step_n_eff"] = n_eff
    );

    return out;

}
