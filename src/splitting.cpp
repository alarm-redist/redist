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
        Plan const&plan, int &region_id_to_split, int min_region_cut_size
        ){


    if(plan.num_multidistricts < 1){
        throw Rcpp::exception("ERROR: Trying to find multidistrict to split when there are none!\n");
    }

    // count total
    int total_multi_ds = 0;

    // make vectors with cumulative d value and region label for later
    std::vector<int> multi_d_vals;
    std::vector<int> region_ids;
    

    // REprintf("[ ");
    // Iterate over all regions
    for(int region_id = 0; region_id < plan.num_regions; region_id++) {

        int d_val = plan.region_dvals(region_id);

        // collect info if multidistrict
        if(d_val > min_region_cut_size){
            // Add that regions d value to the total
            total_multi_ds += d_val;
            // add the count and label to vector
            multi_d_vals.push_back(d_val);
            region_ids.push_back(region_id);
        }
    }
    // REprintf("]\n");

    

    // https://stackoverflow.com/questions/14599057/how-compare-vector-size-with-an-integer
    size_t intendedSize = 1;
    if(region_ids.size() == intendedSize){
        region_id_to_split = region_ids.at(0);
        // probability of picking is 1 and log(1) = 0
        return 0;
    }

    // REprintf("Here in split 2!\n");

    // Now pick an index proportational to d_nk value
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(multi_d_vals.begin(), multi_d_vals.end());

    int idx = d(gen);

    // REprintf("Here in split 3 - picked idx %d !", idx);

    // REprintf("[ ");
    // for (auto i: region_ids){
    //     REprintf("%d ", i);
    // }
    // REprintf("]\n");
    

    region_id_to_split = region_ids.at(idx);

    // REprintf("Here in split 3.5 - picked idx %d !\n", idx);

    double log_prob = std::log(
        static_cast<double>(multi_d_vals.at(idx))
        ) - std::log(
                static_cast<double>(total_multi_ds)
        );

    // REprintf("Here in split 4\n");

    return log_prob;
}


//' Attempts to cut one region into two from a spanning tree and if successful
//' cuts the tree object and returns information on what the two new regions 
//' would be. Does not actually update the plan vertex list.
//'
//' Takes a spanning tree `ust` drawn on a specific region and attempts to cut
//' it to produce two new regions using the generalized splitting procedure
//' outlined <PAPER HERE>. This function is based on `cut_districts` in `smc.cpp`
//' however the crucial difference is even if a cut is successful it does not
//' update the plan. Instead it just returns the information on the two new
//' regions if successful and the cut tree.
//'
//' It will only attempt to create regions where the size is between
//' min_potential_d and max_potential_d (inclusive). So the one district
//' split case is `min_potential_d=max_potential_d=1`.
//'
//' By convention the first new region (`new_region1`) will always be the region
//' with the smaller d-value (although they can be equal).
//'
//' @title Attempt to Find a Valid Spanning Tree Edge to Cut into Two New Regions 
//'
//' @param ust A directed spanning tree passed by reference
//' @param root The root vertex of the spanning tree
//' @param k_param The k parameter from the SMC algorithm, you choose among the top k_param edges
//' @param min_potential_d The smallest potential d value it will try for a cut. 
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
                     int k_param, int min_potential_d, int max_potential_d,
                     const uvec &pop, const arma::subview_col<arma::uword> &region_ids,
                     const int region_id_to_split, int total_region_pop, int total_region_dval,
                     const double lower, const double upper, const double target,
                     int &new_region1_tree_root, int &new_region1_dval, int &new_region1_pop,
                     int &new_region2_tree_root, int &new_region2_dval, int &new_region2_pop
                     ){

    int V = static_cast<int>(region_ids.n_elem);

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

    if(region_ids(root) != region_id_to_split){
        REprintf("Root vertex %u is not in region to split %d!\n", region_ids(root), region_id_to_split);
        throw Rcpp::exception("Root vertex is not in region to split!");
    }


    // Now loop over all valid edges to cut
    for (int i = 1; i <= V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (region_ids(i-1) != region_id_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        int below = pop_below.at(i - 1);
        // Get the population of the other side
        int above = total_region_pop - below;

        // keeps track of the information for the best d_val cut
        double smallest_dev = target * max_potential_d * max_potential_d; // start off with fake maximal value
        bool is_dev2_bigger = false;
        bool cut_is_ok = false;
        int cut_dval;

        // Now try each potential d_nk value from min_potential_d up to max_potential_d
        for(int potential_d = min_potential_d; potential_d <= max_potential_d; potential_d++){

            double dev1 = std::fabs(below - target * potential_d) / (target * potential_d);
            double dev2 = std::fabs(above - target * potential_d) / (target * potential_d);

            // if dev1 is smaller of the two and beats the last one then update
            if(dev1 <= dev2 && dev1 < smallest_dev){
                is_dev2_bigger = true;
                smallest_dev = dev1;
                // check if the split induced fits population bounds
                cut_is_ok = lower * potential_d <= below
                            && below <= upper * potential_d
                            && lower * (total_region_dval - potential_d) <= above
                            && above <= upper * (total_region_dval - potential_d);

                cut_dval = potential_d;
            }else if(dev2 < dev1 && dev2 < smallest_dev){
                is_dev2_bigger = false;
                smallest_dev = dev2;
                // check if the split induced fits population bounds
                cut_is_ok = lower * potential_d <= above
                            && above <= upper * potential_d
                            && lower * (total_region_dval - potential_d) <= below
                            && below <= upper * (total_region_dval - potential_d);

                cut_dval = potential_d;
            }
        }

        if(is_dev2_bigger){
            candidates.push_back(i);
        } else{
            candidates.push_back(-i);
        }

        deviances.push_back(smallest_dev);
        is_ok.push_back(cut_is_ok);
        new_d_val.push_back(cut_dval);
    }

    // if less than k_param candidates immediately reject
    if((int) candidates.size() < k_param){
        return false;
    }

    int idx = r_int(k_param);
    idx = select_k(deviances, idx + 1);
    int cut_at = std::fabs(candidates.at(idx)) - 1;


    // reject sample if not ok
    if (!is_ok.at(idx)){
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
    new_region1_dval = new_d_val.at(idx);
    new_region2_dval = total_region_dval - new_region1_dval;

    if (candidates.at(idx) > 0) { // Means cut below so first vertex is cut_at
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
        const int new_region1_tree_root, const int new_region1_dval, const int new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const int new_region2_pop,
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
    plan.region_dvals(new_region1_id) = new_region1_dval;
    plan.region_added_order.at(new_region1_id) = new_region1_order_added_num;
    plan.region_pops.at(new_region1_id) = new_region1_pop;

    // updates the new region 2
    // New region 2's id is the highest id number so push back
    plan.region_dvals(new_region2_id) = new_region2_dval;
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
//' new ones. Region1 will be set to this id
//' @param new_region_id The id that region2 will be set 
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
        const int old_split_region_id, const int new_region_id,
        const int new_region1_tree_root, const int new_region1_dval, const int new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const int new_region2_pop
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
    int new_region1_id = old_split_region_id;
    // Second new region is new_region_id
    int new_region2_id = new_region_id;

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
//' @param new_region_id The id of the larger of the two new split reasons if successful
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
                 const int new_region_id,
                 std::vector<bool> &visited, std::vector<bool> &ignore, const uvec &pop,
                 const double lower, const double upper, const double target,
                 int k_param, int const min_region_cut_size, int const max_region_cut_size, 
                 bool split_district_only) {

    int V = g.size();

    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_ids(i) != region_id_to_split;
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
    int new_region1_pop, new_region2_pop;


    bool successful_edge_found;

    int max_potential_d;

    if(split_district_only){
        max_potential_d = 1;
    }else{
        max_potential_d = plan.region_dvals(region_id_to_split) - 1;
    }

    int min_potential_d = 1;

    // REprintf("Attempting to split\n");
    // try to get an edge to cut
    successful_edge_found = get_edge_to_cut(ust, root,
                    k_param, min_region_cut_size, max_region_cut_size,
                    pop, plan.region_ids, 
                    region_id_to_split, plan.region_pops.at(region_id_to_split),
                    plan.region_dvals(region_id_to_split),
                    lower, upper, target,
                    new_region1_tree_root, new_region1_dval, new_region1_pop,
                    new_region2_tree_root, new_region2_dval, new_region2_pop);


    if(!successful_edge_found){
        return false;
    }

    // REprintf("Old was %d and new is (%d,%d)\n", max_potential_d+1, new_region1_dval, new_region2_dval);


    // now update the plan with the two new cut regions 
    add_new_regions_to_plan_from_cut(
        ust, plan, split_district_only,
        region_id_to_split, new_region_id,
        new_region1_tree_root, new_region1_dval,  new_region1_pop,
        new_region2_tree_root, new_region2_dval, new_region2_pop
    );


    return true;

}




/*
 * Choose k and multiplier for efficient, accurate sampling
 */
void estimate_cut_k(const Graph &g, int &k, int const last_k, 
                      const std::vector<double> &unnormalized_weights, double thresh,
                      double tol, std::vector<Plan> const &plans_vec, 
                      const uvec &counties,
                      Multigraph &cg, const uvec &pop,
                      int const min_region_cut_size, int const max_region_cut_size,
                      bool split_district_only,
                      double const target, int const verbosity) {
    // sample some spanning trees and compute deviances
    int V = g.size();
    int k_max = std::min(10 + (int) (2.0 * V * tol), last_k + 4); // heuristic
    int N_max = plans_vec.size();
    int N_adapt = std::min(60 + (int) std::floor(5000.0 / sqrt((double)V)), N_max);

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);

    std::vector<std::vector<double>> devs;
    devs.reserve(N_adapt);
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    int idx = 0;
    int max_V = 0;
    Tree ust = init_tree(V);

    // REprintf("max_ok starting at %d\n", max_ok);


    for (int i = 0; i < N_max && idx < N_adapt; i++, idx++) {
        if (unnormalized_weights.at(i) == 0) { // skip if not valid
            idx--;
            continue;
        }

        int n_vtx = V;

        //auto r = plan.region_dvals(arma::span(0, plan.num_regions - 1));

        // arma::span(0,  plan.num_regions - 1);

        // Get the index of the region with the largest dval
        int biggest_region_id; int biggest_dval; int max_valid_dval;

        // if split district only just do remainder 
        if(split_district_only){
            biggest_region_id = plans_vec.at(i).remainder_region;
            biggest_dval = plans_vec.at(i).region_dvals(biggest_region_id);
            max_valid_dval = 1; // max valid dval is 1
        }else{
            biggest_region_id = plans_vec.at(i).region_dvals.head(plans_vec.at(i).num_regions).index_max();
            biggest_dval = plans_vec.at(i).region_dvals(biggest_region_id);
            max_valid_dval = biggest_dval-1;
        }

        double biggest_dval_region_pop = plans_vec.at(i).region_pops.at(biggest_region_id);

        for (int j = 0; j < V; j++) {
            // if not the biggest region mark as ignore
            if (plans_vec.at(i).region_ids(j) != biggest_region_id) {
                ignore[j] = true;
                n_vtx--;
            }
        }
        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                    pop, lower, upper, counties, cg);
        if (result != 0) {
            idx--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, biggest_dval_region_pop, target, min_region_cut_size, max_region_cut_size));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            if (devs.at(idx).at(j) <= tol) { // sorted
                n_ok++;
            } else {
                break;
            }
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max){
            max_ok = n_ok;
            // REprintf("max_ok now %d\n", max_ok);
        }
            
        Rcpp::checkUserInterrupt();
    }

    if (idx < N_adapt) N_adapt = idx; // if rejected too many in last step
    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    for (k = 1; k <= k_max; k++) {
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            double dev = devs.at(i).at(r_int(k));
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k >= k_max) {
        if (verbosity >= 3) {
            Rcout << " [maximum hit; falling back to naive k estimator]";
        }
        k = std::max(max_ok, k_max); // NOTE: used to be max_ok but that seemed too small??
    }

    if (last_k < k_max && k < last_k * 0.6) k = (int) (0.5*k + 0.5*last_k);

    k = std::min(std::max(max_ok + 1, k) + 1 - (distr_ok(k) > 0.99) + (thresh == 1),
                 max_V - 1);
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
//' @param parent_index_vec A vector used to track the index of the previous plan
//' sampled that was successfully split. The value of `parent_index_vec[i]` is the
//' index of the old plan from which the new plan `new_plans_vec[i]` was
//' successfully split from. In other words `new_plans_vec[i]` is equal to
//' `attempt_region_split(old_plans_vec[parent_index_vec[i]], ...)`
//' @param unnormalized_sampling_weights A vector of weights used to sample indices
//' of the `old_plans_vec`. The value of `unnormalized_sampling_weights[i]` is
//' the unnormalized probability that index i is selected
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
//' plans. This is equal to `unique(parent_index_vec)`
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
//'    - The `parent_index_vec` is updated to contain the indices of the parents of the
 //'    new plans
//'    - If two new valid regions are split then the new_region_ids is updated so the
//'    first entry is the first new region and the second entry is the second new region
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
        Rcpp::IntegerMatrix::Column parent_index_vec,
        const std::vector<double> &unnormalized_sampling_weights,
        Rcpp::IntegerMatrix::Column draw_tries_vec,
        Rcpp::IntegerMatrix::Column parent_unsuccessful_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        umat &ancestors, const std::vector<int> &lags,
        double const lower, double const upper, double const target,
        int const k_param, int const min_region_cut_size, int const max_region_cut_size, 
        bool const split_district_only,
        RcppThread::ThreadPool &pool,
        int const verbosity, int const diagnostic_level
                ) {
    // important constants
    const int V = g.size();
    const int M = old_plans_vec.size();

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

    bool print_weights = M < 12 && verbosity > 1;
    if(print_weights){
    std::vector<double> p = index_sampler.probabilities();
    Rprintf("Unnormalized weights are: ");
    for (auto w : unnormalized_sampling_weights){
        Rprintf("%.4f, ", w);
    }

    Rprintf("\n");
    Rprintf("Normalized weights are: ");
    
    for (auto prob : p){
        if(print_weights) Rprintf("%.4f, ", prob);

    }
    Rprintf("\n");
    }

    // Because of multithreading we have to add specific checks for if the user
    // wants to quit the program
    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50; // check for interrupts every _ iterations

    // The new region in the split plans is the number of regions in a split plan minus
    // one so the number of regions in a presplit plan
    int new_region_id = old_plans_vec.at(0).num_regions;

    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // REprintf("Plan %d\n", i);
        int reject_ct = 0;
        bool ok = false;
        int idx;
        int region_id_to_split;

        Tree ust = init_tree(V);
        std::vector<bool> visited(V);
        std::vector<bool> ignore(V);
        while (!ok) {
            // increase the number of tries for particle i by 1
            draw_tries_vec[i]++;
            //draw_tries_vec(i)++;
            // REprintf("Plan %d, iter %d\n", i, iters[i]);
            // use weights to sample previous plan
            idx = index_sampler(gen);
            // REprintf("Picked idx %d, from vector of size %u\n", i, old_plans_vec.size());

            // We want the data from the old plans vec but we can't modify it since they
            // point to the same matrix so we do shallow copy             
            new_plans_vec.at(i).shallow_copy(old_plans_vec.at(idx));

            // REprintf("Now picking region to split\n");
            if(split_district_only){
                // if just doing district splits just use remainder region
                region_id_to_split = new_plans_vec.at(i).remainder_region;
            }else{
                // if generalized split pick a region to try to split
                choose_multidistrict_to_split(
                    old_plans_vec.at(idx), region_id_to_split, min_region_cut_size);
            }

            auto plan_specific_max_region_cut_size = std::min(
                max_region_cut_size, // max_region_cut_size 
                (int) old_plans_vec.at(idx).region_dvals(region_id_to_split)
            );

            // REprintf("Now splitting\n");

            // Now try to split that region
            ok = attempt_region_split(g, ust, counties, cg,
                                      new_plans_vec.at(i), region_id_to_split,
                                      new_region_id,
                                      visited, ignore, pop,
                                      lower, upper, target,
                                      k_param, min_region_cut_size, 
                                      plan_specific_max_region_cut_size,
                                      split_district_only);


            // bad sample; try again
            if (!ok) {
                 // update unsuccessful try
                RcppThread::checkUserInterrupt(++reject_ct % reject_check_int == 0);

                // if diagnostic level 2 or higher get unsuccessful count 
                if(diagnostic_level >= 2){
                    // not atomic so technically not thread safe but doesn't seem to differ in practice
                    parent_unsuccessful_tries_vec[idx]++;
                }
                continue;
            }

            // else leave the while loop


        }

        // record index of new plan's parent
        parent_index_vec[i] = idx;


        // ORIGINAL SMC CODE I DONT KNOW WHAT THIS DOES
        // save ancestors/lags
        for (int j = 0; j < n_lags; j++) {
            if (dist_ctr <= lags[j]) {
                ancestors_new(i, j) = i;
            } else {
                ancestors_new(i, j) = ancestors(idx, j);
            }
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();

    // REprintf("Copying now!\n");

    // now replace the old plans with the new ones
    for(int i=0; i < M; i++){
        old_plans_vec.at(i).shallow_copy(new_plans_vec.at(i));
    }


    // now compute acceptance rate and unique parents and original ancestors
    accept_rate = M / (1.0 * sum(draw_tries_vec));

    // Get number of unique parents
    std::set<int> unique_parents(parent_index_vec.begin(), parent_index_vec.end());
    n_unique_parent_indices = unique_parents.size();
    if (verbosity >= 3) {
        Rprintf("%.2f acceptance rate and %d unique parent indices sampled!\n",
                100.0 * accept_rate, (int) n_unique_parent_indices);
    }

    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;

}
