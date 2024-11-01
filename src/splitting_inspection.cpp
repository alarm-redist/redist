/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes various aspects of the splitting procedure
*   to R to allow for inspection for diagnostic purposes.
********************************************************/

#include "splitting_inspection.h"



List perform_a_valid_region_split(
        List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int k_param, int region_id_to_split,
        double target, double lower, double upper,
        int N, int num_regions, int num_districts,
        std::vector<int> region_ids, std::vector<int> region_dvals,
        std::vector<double> region_pops,
        bool split_district_only, bool verbose
){
    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    double total_pop = sum(pop);

    // Create a plan object
    Plan plan = Plan(V, N, total_pop);

    // fill in the plan
    plan.num_regions = num_regions;
    plan.num_districts = num_districts;
    plan.num_multidistricts = plan.num_regions - plan.num_districts;
    plan.region_ids = region_ids;
    plan.region_dvals = region_dvals;
    plan.region_pops = region_pops;


    if(verbose){
        plan.Rprint();
    }

    // Create tree related stuff
    int root;
    Tree ust = init_tree(V);
    Tree pre_split_ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);


    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan.region_ids.at(i) != region_id_to_split;
    }



    // Counts the number of split attempts
    bool successful_split_made = false;
    int try_counter {0};


    // splitting related params
    int new_region1_tree_root, new_region2_tree_root;
    int new_region1_dval, new_region2_dval;
    double new_region1_pop, new_region2_pop;

    // Keep running until done
    while(!successful_split_made){
        clear_tree(ust);
        // Get a tree
        int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
        // Return unsuccessful if tree not drawn
        if (result != 0){
            try_counter++;
            continue;
        }

        // copy uncut tree
        pre_split_ust = ust;

        // Try to make a cut
        successful_split_made = get_edge_to_cut(ust, root,
                                                k_param,  split_district_only, pop,
                                                plan, region_id_to_split,
                                                lower, upper, target,
                                                new_region1_tree_root, new_region1_dval, new_region1_pop,
                                                new_region2_tree_root, new_region2_dval, new_region2_pop);
        try_counter++;
        // increase the counter by 1
        if(verbose && false){
            Rcout << "Attempt " << try_counter << "\n";
        }

    }


    if(verbose){
        plan.Rprint();
    }

    // update
    // TODO: make this optional
    int new_region1_id, new_region2_id;

    update_plan_from_cut(
        ust, plan, split_district_only,
        region_id_to_split,
        new_region1_tree_root, new_region1_dval,  new_region1_pop,
        new_region2_tree_root, new_region2_dval, new_region2_pop,
        new_region1_id, new_region2_id
    );

    if(verbose){
        plan.Rprint();
    }

    List out = List::create(
        _["num_attempts"] = try_counter,
        _["region_id_that_was_split"] = region_id_to_split,
        _["region_dvals"] = plan.region_dvals,
        _["plan_vertex_ids"] = plan.region_ids,
        _["pops"] = plan.region_pops,
        _["num_regions"] = plan.num_regions,
        _["num_districts"] = plan.num_districts,
        _["tree"] = ust,
        _["uncut_tree"] = pre_split_ust,
        _["new_region1_id"] = new_region1_id,
        _["new_region1_tree_root"] = new_region1_tree_root,
        _["new_region1_dval"] = new_region1_dval,
        _["new_region1_pop"] = new_region1_pop,
        _["new_region2_id"] = new_region2_id,
        _["new_region2_tree_root"] = new_region2_tree_root,
        _["new_region2_dval"] = new_region2_dval,
        _["new_region2_pop"] = new_region2_pop
    );

    return out;
}


/*
 * New splitter test function. Performs cut tree function on the entire map
 */
List get_successful_proposed_cut(int N, List adj_list, const arma::uvec &counties, const arma::uvec &pop,
                                        double target, double lower, double upper,
                                        bool split_district_only, bool verbose
) {


    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    double total_pop = sum(pop);

    // Create a plan object
    Plan plan = Plan(V, N, total_pop);


    // Create tree related stuff
    Tree ust = init_tree(V);
    Tree pre_split_ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);

    auto region_id_to_split = 0;

    // Get a uniform spanning tree drawn on that region
    int root;




    // set k param for splitter
    auto k_param {6};

    // split the only region
    int region_to_split = 0;


    bool successful_split_made = false;
    int try_counter {0};



    int new_region1_tree_root, new_region2_tree_root;
    int new_region1_dval, new_region2_dval;
    double new_region1_pop, new_region2_pop;

    // Keep running until done
    while(!successful_split_made){
        clear_tree(ust);
        // Get a tree
        int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
        // Return unsuccessful if tree not drawn
        if (result != 0){
            try_counter++;
            continue;
        }

        // copy uncut tree
        pre_split_ust = ust;

        // Try to make a cut
        successful_split_made = get_edge_to_cut(ust, root,
                                                k_param,  split_district_only, pop,
                                                plan, region_id_to_split,
                                                lower, upper, target,
                                                new_region1_tree_root, new_region1_dval, new_region1_pop,
                                                new_region2_tree_root, new_region2_dval, new_region2_pop);
        try_counter++;
        // increase the counter by 1
        if(verbose){
            Rcout << "Attempt " << try_counter << "\n";
        }

    }


    if(verbose){
        plan.Rprint();
    }



    List out = List::create(
        _["num_attempts"] = try_counter,
        _["region_dvals"] = plan.region_dvals,
        _["plan"] = plan.region_ids,
        _["tree"] = ust,
        _["uncut_tree"] = pre_split_ust,
        _["pops"] = plan.region_pops,
        _["new_region1_tree_root"] = new_region1_tree_root,
        _["new_region1_dval"] = new_region1_dval,
        _["new_region1_pop"] = new_region1_pop,
        _["new_region2_tree_root"] = new_region2_tree_root,
        _["new_region2_dval"] = new_region2_dval,
        _["new_region2_pop"] = new_region2_pop
    );

    return out;
}
