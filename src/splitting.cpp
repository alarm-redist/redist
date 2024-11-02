/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/11
* Purpose: Functions for Splitting Trees and Plans
********************************************************/

#include "splitting.h"


//' Selects a multidistrict with probability proportional to its d_nk value and
//' returns the log probability of the selected region
//'
//' Given a plan object with at least one multidistrict this function randomly
//' selects a multidistrict with probability proporitional to its d_nk value
//' (relative to all multidistricts) and returns the log of the probability that
//' region was chosen.
//'
//'
//' @title Choose multidistrict to split
//'
//' @param plan A plan object
//' @param region_to_split an integer that will be updated by reference with the
//' id number of the region selected to split
//'
//' @details No modifications to inputs made
//'
//' @return the region level graph
//'
double choose_multidistrict_to_split(
        Plan const&plan, int &region_id_to_split){

    if(plan.num_multidistricts < 1){
        Rprintf("ERROR: Trying to find multidistrict to split when there are none!\n");
    }

    // count total
    int total_multi_ds = 0;

    // make vectors with cumulative d value and region label for later
    std::vector<int> multi_d_vals;
    std::vector<int> region_ids;

    // Iterate over all regions
    for(int region_id = 0; region_id < plan.num_regions; region_id++) {

        int d_val = plan.region_dvals.at(region_id);

        // collect info if multidistrict
        if(d_val > 1){
            // Add that regions d value to the total
            total_multi_ds += d_val;
            // add the count and label to vector
            multi_d_vals.push_back(d_val);
            region_ids.push_back(region_id);
        }
    }

    // https://stackoverflow.com/questions/14599057/how-compare-vector-size-with-an-integer
    size_t intendedSize = 1;
    if(region_ids.size() == intendedSize){
        region_id_to_split = region_ids.at(0);
        // probability of picking is 1 and log(1) = 0
        return 0;
    }


    // Now pick an index proportational to d_nk value
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(multi_d_vals.begin(), multi_d_vals.end());

    int idx = d(gen);

    region_id_to_split = region_ids[idx];
    double log_prob = std::log(
        static_cast<double>(multi_d_vals[idx])
        ) - std::log(
                static_cast<double>(total_multi_ds)
        );

    return log_prob;
}


//' Attempts to cut one region into two from a spanning tree and if successful
//' returns information on what the two new regions would be. Does not actually
//' update the plan
//'
//' Takes a spanning tree `ust` drawn on a specific region and attempts to cut
//' it to produce two new regions using the generalized splitting procedure
//' outlined <PAPER HERE>. This function is based on `cut_districts` in `smc.cpp`
//' however the crucial difference is even if a cut is successful it does not
//' update the plan. Instead it just returns the information on the two new
//' regions if successful and the vertices to use to update the plans.
//'
//' Depending on the value of max_potential_d will only attempt to split off
//' a single district or allows for more general splits.
//'
//' By convention the first new region (`new_region1`) will always be the region
//' with the smaller d-value (although they can be equal).
//'
//' @title Attempt to Cut Region Tree into Two New Regions
//'
//' @param ust A directed spanning tree passed by reference
//' @param root The root vertex of the spanning tree
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param max_potential_d The largest potential d value it will try for a cut. Setting this to 
//' 1 will result in only 1 district splits. 
//' @param pop A vector of the population associated with each vertex in `g`
//' @param region_ids A vector mapping 0 indexed vertices to their region id number
//' @param region_id_to_split The id of the region in the plan object we're attempting to split
//' @param total_region_pop The total population of the region being split 
//' @param total_region_dval The dval of the region being split 
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param new_region1_tree_root The index of the root of tree associated with
//' the first new region (if the tree cut was successful)
//' @param new_region1_dval The d-value of the first new region (if the tree cut
//'  was successful)
//' @param new_region1_pop The population of the first new region (if the tree cut
//' was successful)
//' @param new_region2_tree_root The index of the root of tree associated with
//' the second new region (if the tree cut was successful)
//' @param new_region2_dval The d-value of the second new region (if the tree cut
//'  was successful)
//' @param new_region2_pop The population of the second new region (if the tree cut
//' was successful)
//'
//' @details Modifications
//'    - If two new valid regions are split then the tree `ust` is cut into two
//'    distjoint pieces
//'    - If two new valid regions are split then the 6 `new_region` inputs are all
//'    updated by reference with the values associated with the new regions
//'
//' @return True if two valid regions were successfully split, false otherwise
//'
bool get_edge_to_cut(Tree &ust, int root,
                     int k_param, int max_potential_d,
                     const uvec &pop, const std::vector<int> &region_ids,
                     const int region_id_to_split, double total_region_pop, int total_region_dval,
                     const double lower, const double upper, const double target,
                     int &new_region1_tree_root, int &new_region1_dval, double &new_region1_pop,
                     int &new_region2_tree_root, int &new_region2_dval, double &new_region2_pop
                     ){

    // total_region_dval

    int V = static_cast<int>(region_ids.size());

    // Rcout << "For " << region_id_to_split << " Total pop is " << total_pop << " and d_nk is " << num_final_districts << "\n";

    // create list that points to parents & computes population below each vtx
    std::vector<int> pop_below(V, 0);
    std::vector<int> parent(V);
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

    if(region_ids.at(root) != region_id_to_split){
        Rcout << "Root vertex is not in region to split!";
    }

    // Now loop over all valid edges to cut
    for (int i = 1; i <= V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (region_ids.at(i-1) != region_id_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        double below = pop_below.at(i - 1);
        // Get the population of the other side
        double above = total_region_pop - below;

        // vectors to keep track of value for each d_{n,k}
        std::vector<double> devs(max_potential_d);
        std::vector<bool> dev2_bigger(max_potential_d);
        std::vector<bool> are_ok(max_potential_d);


        // Now try each potential d_nk value from 1 up to d_n-1,k -1
        // remember to correct for 0 indexing
        for(int potential_d = 1; potential_d <= max_potential_d; potential_d++){

            double dev1 = std::fabs(below - target * potential_d);
            double dev2 = std::fabs(above - target * potential_d);

            // If dev1 is smaller then we assign d_nk to below
            if (dev1 < dev2) {
                devs[potential_d-1] = dev1;
                dev2_bigger[potential_d-1] = true;
                // check in bounds
                are_ok[potential_d-1] =    lower * potential_d <= below
                                        && below <= upper * potential_d
                                        && lower * (total_region_dval - potential_d) <= above
                                        && above <= upper * (total_region_dval - potential_d);
            } else { // Else if dev2 is smaller we assign d_nk to above
                devs[potential_d-1] = dev2;
                dev2_bigger[potential_d-1] = false;
                are_ok[potential_d-1] =    lower * potential_d <= above
                                        && above <= upper * potential_d
                                        && lower * (total_region_dval - potential_d) <= below
                                        && below <= upper * (total_region_dval - potential_d);
            }

        }

        // Now find the value of d_{n,k} that has the smallest deviation
        std::vector<double>::iterator result = std::min_element(devs.begin(), devs.end());
        int best_potential_d = std::distance(devs.begin(), result);


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


    // reject sample if not ok
    if (!is_ok[idx]){
        return false;
    }


    // find index of node to cut at
    std::vector<int> *siblings = &ust[parent[cut_at]];
    int length = siblings->size();
    int j;
    for (j = 0; j < length; j++) {
        if ((*siblings)[j] == cut_at) break;
    }

    // remove edge from tree
    siblings->erase(siblings->begin()+j);
    parent[cut_at] = -1;

    // Get dval for two new regions
    new_region1_dval = new_d_val[idx];
    new_region2_dval = total_region_dval - new_region1_dval;


    if (candidates[idx] > 0) { // Means cut below so first vertex is cut_at
        new_region1_tree_root = cut_at;
        new_region2_tree_root = root;

        // Set the new populations
        new_region1_pop = pop_below.at(cut_at);
        new_region2_pop = total_region_pop - new_region1_pop;

    //     return pop_below.at(cut_at);
    } else { // Means cut above so first vertex is root
        new_region1_tree_root = root;
        new_region2_tree_root = cut_at;

        // Set the new populations
        new_region2_pop = pop_below.at(cut_at);
        new_region1_pop = total_region_pop - new_region2_pop;

    }

    // Now make it so that the first region is always the smallest
    // So if the dval of region 1 is bigger swap the two
    if(new_region1_dval > new_region2_dval){
        std::swap(new_region1_dval, new_region2_dval);
        std::swap(new_region1_pop, new_region2_pop);
        std::swap(new_region1_tree_root, new_region2_tree_root);
    }

    return true;

};






//' Updates a `Plan` object using a cut tree
//'
//' Takes a cut spanning tree `ust` and variables on the two new regions
//' induced by the cuts and updates `plan` with information on those two
//' new regions. Assumes that the plan attributes already have the correct
//' size and accessing either of the region ids won't create issues.
//'
//' It also sets `plan.remainder_region` equal to `new_region2_id` if 
//' split_district_only is true. 
//'
//'
//' @title Update plan regions from cut tree
//'
//' @param ust A cut (ie has two partition pieces) directed spanning tree
//' passed by reference
//' @param plan A plan object
//' @param split_district_only Whether or not this was split according to a 
//' one district split scheme (as in does the remainder need to be updated)
//' @param new_region1_tree_root The vertex of the root of one piece of the cut
//' tree. This always corresponds to the region with the smaller dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region1_dval The dval associated with the new region 1
//' @param new_region1_pop The population associated with the new region 1
//' @param new_region2_tree_root The vertex of the root of other piece of the cut
//' tree. This always corresponds to the region with the bigger dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region2_dval The dval associated with the new region 2
//' @param new_region2_pop The population associated with the new region 2
//' @param new_region1_id The id the new region 1 was assigned in the plan
//' @param new_region2_id The id the new region 2 was assigned in the plan
//'
//' @details Modifications
//'    - `plan` is updated in place with the two new regions
//'
void update_plan_from_cut(
        Tree &ust, Plan &plan, bool split_district_only,
        const int new_region1_tree_root, const int new_region1_dval, const double new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const double new_region2_pop,
        const int new_region1_id,  const int new_region2_id
){

    // make the first new region get max plus one
    plan.region_order_max++;
    int new_region1_order_added_num = plan.region_order_max;

    // Now make the second new region max plus one again
    plan.region_order_max++;
    int new_region2_order_added_num = plan.region_order_max;


    // Now update the two cut portions
    assign_region(ust, plan, new_region1_tree_root, new_region1_id);
    assign_region(ust, plan, new_region2_tree_root, new_region2_id);

    // Now update the region level information

    // updates the new region 1
    plan.region_dvals.at(new_region1_id) = new_region1_dval;
    plan.region_added_order.at(new_region1_id) = new_region1_order_added_num;
    plan.region_pops.at(new_region1_id) = new_region1_pop;

    // updates the new region 2
    // New region 2's id is the highest id number so push back
    plan.region_dvals.at(new_region2_id) = new_region2_dval;
    plan.region_added_order.at(new_region2_id) = new_region2_order_added_num;
    plan.region_pops.at(new_region2_id) = new_region2_pop;

    // If district split only then set the remainder region as well
    if(split_district_only){
        // update the remainder region value if needed
        plan.remainder_region = new_region2_id;
    }

    // done
    return;

};


//' Creates new regions and updates the `Plan` object using a cut tree
//'
//' Takes a cut spanning tree `ust` and variables on the two new regions
//' induced by the cuts and creates space/updates the information on those
//' two new regions in the `plan` object. This function increases the number
//' of regions aspect by 1 and updates the region level information and all
//' other variables changed by adding a new region. 
//'
//' It also sets `plan.remainder_region` equal to `new_region2_id` if 
//' split_district_only is true. 
//'
//'
//' @title Create and update new plan regions from cut tree
//'
//' @param ust A cut (ie has two partition pieces) directed spanning tree
//' passed by reference
//' @param plan A plan object
//' @param split_district_only Whether or not this was split according to a 
//' one district split scheme (as in does the remainder need to be updated)
//' @param old_split_region_id The id of the region that was split into the two
//' new ones 
//' @param new_region1_tree_root The vertex of the root of one piece of the cut
//' tree. This always corresponds to the region with the smaller dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region1_dval The dval associated with the new region 1
//' @param new_region1_pop The population associated with the new region 1
//' @param new_region2_tree_root The vertex of the root of other piece of the cut
//' tree. This always corresponds to the region with the bigger dval (allowing
//' for the possiblity the dvals are equal).
//' @param new_region2_dval The dval associated with the new region 2
//' @param new_region2_pop The population associated with the new region 2
//' @param new_region1_id The id the new region 1 was assigned in the plan
//' @param new_region2_id The id the new region 2 was assigned in the plan
//'
//' @details Modifications
//'    - `plan` is updated in place with the two new regions
//'    - `new_region1_id` is set to the id new region1 was assigned
//'    which is just the `old_split_region_id`
//'    - `new_region2_id` is set to the id new region2 was assigned 
//'    which is just `plan.num_regions-1`
//'
void add_new_regions_to_plan_from_cut(
        Tree &ust, Plan &plan, bool split_district_only,
        const int old_split_region_id,
        const int new_region1_tree_root, const int new_region1_dval, const double new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const double new_region2_pop,
        int &new_region1_id,  int &new_region2_id
){
    // if successful then update the plan
    // update plan with new regions
    plan.num_regions++; // increase region count by 1
    plan.num_multidistricts--; // Decrease by one to avoid double counting later

    // Create info for two new districts

    // Set label and count depending on if district or multi district
    if(new_region1_dval == 1){
        plan.num_districts++;
    }else{
        plan.num_multidistricts++;
    }

    // Now do it for second region
    if(new_region2_dval == 1){
        plan.num_districts++;
    }else{
        plan.num_multidistricts++;
    }

    // make the first new region have the same integer id as the split region
    new_region1_id = old_split_region_id;
    // Second new region has id of the new number of regions minus 1
    new_region2_id = plan.num_regions - 1;

    // Now resize the region level attributes 
    plan.region_dvals.resize(plan.num_regions, -1);
    plan.region_added_order.resize(plan.num_regions, -1);
    plan.region_pops.resize(plan.num_regions, -1.0);

    // now update using tree
    update_plan_from_cut(
            ust, plan, split_district_only,
            new_region1_tree_root, new_region1_dval,  new_region1_pop,
            new_region2_tree_root, new_region2_dval, new_region2_pop,
            new_region1_id, new_region2_id
        );

}



//' Attempts to split a multi-district within a plan into two new regions with
//' valid population bounds (has option for one district splits only)
//'
//' Given a plan this attempts to split a multi-district in it into two new
//' regions where both regions have valid population bounds. (If
//' `split_district_only` is true it will only attempt to split off a district
//' and a remainder). It does this by drawing a spanning tree uniformly at random
//' then calling `get_edge_to_cut` on  that. If the a valid cut is found it
//' then calls `update_plan_from_cut` to update the plan accordingly. If not
//' successful returns false and does nothing.
//'
//' This is based on the `split_map` function in `smc.cpp`
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
//' @param region_id_to_split The label of the region in the plan object we're attempting to split
//' @param new_region_ids A vector that will be updated by reference to contain the names of
//' the two new split regions if function is successful.
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' Target population (probably Total population of map/Num districts)
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param split_district_only Whether or not to only allow for one district
//' split. If `true` then only splits off districts.
//'
//' @details Modifications
//'    - If two new valid regions are split then the plan object is updated accordingly
//'    - If two new valid regions are split then the new_region_ids is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'
//' @return True if two valid regions were split off false otherwise
//'
bool attempt_region_split(const Graph &g, Tree &ust, const uvec &counties, Multigraph &cg,
                 Plan &plan, const int region_id_to_split,
                 std::vector<int> &new_region_ids,
                 std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                 double &lower, double upper, double target,
                 int k_param, bool split_district_only) {

    int V = g.size();

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_ids.at(i) != region_id_to_split;
    }

    // Get a uniform spanning tree drawn on that region
    int root;
    clear_tree(ust);
    // Get a tree
    int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
    // Return unsuccessful if tree not drawn
    if (result != 0) return false;


    // splitting related params
    int new_region1_tree_root, new_region2_tree_root;
    int new_region1_dval, new_region2_dval;
    double new_region1_pop, new_region2_pop;


    bool successful_edge_found;

    int max_potential_d;

    if(split_district_only){
        max_potential_d = 1;
    }else{
        max_potential_d = plan.region_dvals.at(region_id_to_split) - 1;
    }

    // try to get an edge to cut
    successful_edge_found = get_edge_to_cut(ust, root,
                    k_param, max_potential_d,
                    pop, plan.region_ids, 
                    region_id_to_split, plan.region_pops.at(region_id_to_split),
                    plan.region_dvals.at(region_id_to_split),
                    lower, upper, target,
                    new_region1_tree_root, new_region1_dval, new_region1_pop,
                    new_region2_tree_root, new_region2_dval, new_region2_pop);


    if(!successful_edge_found){
        return false;
    }


    // now update the plan with the two new cut regions 
    add_new_regions_to_plan_from_cut(
        ust, plan, split_district_only,
        region_id_to_split,
        new_region1_tree_root, new_region1_dval,  new_region1_pop,
        new_region2_tree_root, new_region2_dval, new_region2_pop,
        new_region_ids[1], new_region_ids[2]
    );

    return true;

}







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
//' @param parent_unsuccessful_tries_vec A vector used to keep track of how many times the
//' previous rounds plans were sampled and unsuccessfully split. The value
//' `parent_unsuccessful_tries_vec[i]` represents how many times `old_plans_vec[i]` was sampled
//' and then unsuccessfully split while creating all `M` of the new plans.
//' THIS MAY NOT BE THREAD SAFE
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
//' @param split_district_only Whether or not to only allow for single district
//' splits. If set to `true` will only attempt to split off one district at a
//' time
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
//'    - If two new valid regions are split then the new_region_ids is updated so the
//'    first entry is the first new region and the second entry is the second new region
//'    - The `normalized_weights_to_fill_in` is updated to contain the normalized
//'    probabilities the index sampler used. This is only collected for diagnostics
//'    at this point and should just be equal to `unnormalized_sampling_weights`
//'    divided by `sum(unnormalized_sampling_weights)`
//'    - The `draw_tries_vec` is updated to contain the number of tries for each
//'    of the new plans
//'    - The `parent_unsuccessful_tries_vec` is updated to contain the number of unsuccessful
//'    samples of the old plans
//'    - The `accept_rate` is updated to contain the average acceptance rate for
//'    this iteration
//'    - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
//'    with the unique number of parents and original ancestors for all the new
//'    plans respectively
//'    - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
//'    I DO NOT KNOW WHAT IT MEANS
//'
//' @keywords internal
void generalized_split_maps(
        const Graph &g, const uvec &counties, Multigraph &cg, const uvec &pop,
        std::vector<Plan> &old_plans_vec, std::vector<Plan> &new_plans_vec,
        std::vector<int> &original_ancestor_vec,
        std::vector<int> &parent_vec,
        const std::vector<int> &prev_ancestor_vec,
        const std::vector<double> &unnormalized_sampling_weights,
        std::vector<double> &normalized_weights_to_fill_in,
        std::vector<int> &draw_tries_vec,
        std::vector<int> &parent_unsuccessful_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        int &n_unique_original_ancestors,
        umat &ancestors, const std::vector<int> &lags,
        double lower, double upper, double target,
        int k_param, bool split_district_only,
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
        int region_id_to_split;
        std::vector<int> new_region_ids(2, -1);

        Tree ust = init_tree(V);
        std::vector<bool> visited(V);
        std::vector<bool> ignore(V);
        while (!ok) {
            // increase the iters count by one
            iters[i]++;
            // use weights to sample previous plan
            idx = index_sampler(gen);
            Plan proposed_new_plan = old_plans_vec[idx];

            if(split_district_only){
                // if just doing district splits just use remainder region
                region_id_to_split = proposed_new_plan.remainder_region;
            }else{
                // if generalized split pick a region to try to split
                choose_multidistrict_to_split(
                    old_plans_vec[idx], region_id_to_split);
            }

            // Now try to split that region
            ok = attempt_region_split(g, ust, counties, cg,
                                      proposed_new_plan, region_id_to_split,
                                      new_region_ids,
                                      visited, ignore, pop,
                                      lower, upper, target,
                                      k_param, split_district_only);


            // bad sample; try again
            if (!ok) {
                // THIS MAY NOT BE THREAD SAFE
                parent_unsuccessful_tries_vec[idx]++; // update unsuccessful try
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);
                continue;
            }

            // else update the new plan and leave the while loop
            new_plans_vec[i] = proposed_new_plan;

        }

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

    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;

}
