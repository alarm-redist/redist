/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Run SMC with MCMC merge-split steps mixed in
********************************************************/

#include "smc_and_mcmc.h"


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
//' Depending on the value of split_district_only will only attempt to split off
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
//' @param split_district_only If true then only tries to split a district, if false allows for
//' arbitrary region splits
//' @param pop A vector of the population associated with each vertex in `g`
//' @param plan A plan object
//' @param region_id_to_split The id of the region in the plan object we're attempting to split
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
//'
//'
bool get_edge_to_cut(Tree &ust, int root,
                     int k_param, bool split_district_only,
                     const uvec &pop, const Plan &plan, const int region_id_to_split,
                     const double lower, const double upper, const double target,
                     int &new_region1_tree_root, int &new_region1_dval, double &new_region1_pop,
                     int &new_region2_tree_root, int &new_region2_dval, double &new_region2_pop
                     ){
    // Get population of region being split
    double total_pop = plan.region_pops.at(region_id_to_split);

    // Get the number of final districts in region being split
    int num_final_districts = plan.region_dvals.at(region_id_to_split);

    // Get the largest possible d value we can try
    int max_potential_d;

    if(split_district_only){
        // if only splitting districts its just 1 bc then it will only try
        // a potential d of 1 later on
        max_potential_d = 1;
    }else{
        // if arbitrary splits then d value of region minus 1
        max_potential_d = num_final_districts - 1;
    }

    // Rcout << "For " << region_id_to_split << " Total pop is " << total_pop << " and d_nk is " << num_final_districts << "\n";

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

    if(plan.region_num_ids.at(root) != region_id_to_split){
        Rcout << "Root vertex is not in region to split!";
    }

    // Now loop over all valid edges to cut
    for (int i = 1; i <= plan.V; i++) { // 1-indexing here
        // Ignore any vertex not in this region or the root vertex as we wont be cutting those
        if (plan.region_num_ids.at(i-1) != region_id_to_split || i - 1 == root) continue;

        // Get the population of one side of the partition removing that edge would cause
        double below = pop_below.at(i - 1);
        // Get the population of the other side
        double above = total_pop - below;

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
    new_region2_dval = num_final_districts - new_region1_dval;



    if (candidates[idx] > 0) { // Means cut below so first vertex is cut_at
        new_region1_tree_root = cut_at;
        new_region2_tree_root = root;

        // Set the new populations
        new_region1_pop = pop_below.at(cut_at);
        new_region2_pop = total_pop - new_region1_pop;

    //     return pop_below.at(cut_at);
    } else { // Means cut above so first vertex is root
        new_region1_tree_root = root;
        new_region2_tree_root = cut_at;

        // Set the new populations
        new_region2_pop = pop_below.at(cut_at);
        new_region1_pop = total_pop - new_region2_pop;

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




void update_plan_from_cut(
        Tree &ust, Plan &plan,
        const int old_region_id,
        const int new_region1_tree_root, const int new_region1_dval, const double new_region1_pop,
        const int new_region2_tree_root, const int new_region2_dval, const double new_region2_pop,
        int &new_region1_id,  int &new_region2_id
){

    // update plan with new regions
    plan.num_regions++; // increase region count by 1
    plan.num_multidistricts--; // Decrease by one to avoid double counting later

    // Create info for two new districts
    std::string new_region_label1;
    std::string new_region_label2;

    std::string region_to_split = plan.region_str_labels.at(old_region_id);

    // Set label and count depending on if district or multi district
    if(new_region1_dval == 1){
        plan.num_districts++;
        // if district then string label just adds district number
        new_region_label1 = region_to_split + "." + std::to_string(plan.num_districts);
    }else{
        plan.num_multidistricts++;
        // if region then just add current region number
        new_region_label1 = region_to_split + ".R" + std::to_string(old_region_id);
    }

    // Now do it for second region
    if(new_region2_dval == 1){
        plan.num_districts++;
        // if district then string label just adds district number
        new_region_label2 = region_to_split + "." + std::to_string(plan.num_districts);
    }else{
        plan.num_multidistricts++;
        // if region then just add current region number
        new_region_label2 = region_to_split + ".R" + std::to_string(plan.num_regions - 1);
    }

    // TODO: Figure out how to do this to preserve the order districts were added

    // make the first new region have the same integer id
    new_region1_id = old_region_id;
    // Second new region has id of the new number of regions minus 1
    new_region2_id = plan.num_regions - 1;


    // Now update the two cut portions
    assign_region(ust, plan, new_region1_tree_root, new_region1_id);
    assign_region(ust, plan, new_region2_tree_root, new_region2_id);

    // Now update the region level information

    // Add the new region 1
    // New region 1 has the same id number as old region so update that
    plan.region_dvals.at(new_region1_id) = new_region1_dval;
    plan.region_str_labels.at(new_region1_id) = new_region_label1;
    plan.region_pops.at(new_region1_id) = new_region1_pop;

    // Add the new region 2
    // New region 2's id is the highest id number so push back
    plan.region_dvals.push_back(new_region2_dval);
    plan.region_str_labels.push_back(new_region_label2);
    plan.region_pops.push_back(new_region2_pop);

    // done
    return;

};
