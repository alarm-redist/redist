/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes various aspects of the splitting procedure
*   to R to allow for manual splitting of plans.
********************************************************/

#include "manual_splitting.h"


//' Draws a spanning tree uniformly at random on a region and returns it
//'
//' Draws a spanning tree uniformly at random on a region of a plan using
//' Wilson's algorithm. 
//'
//' @title Draw a uniformly random spanning tree on a region of a plan
//'
//'
//' @param adj_list A 0-indexed adjacency list representing the undirected graph
//' which represents the underlying map the plans are to be drawn on
//' @param counties Vector of county labels of each vertex in `g`
//' @param pop A vector of the population associated with each vertex in `g`
//' @param ndists The number of districts the final plans will have
//' @param num_regions The number of regions in the inputted plan
//' @param num_districts The number of districts in the inputted plan
//' @param region_id_to_draw_tree_on The id of the region in the plan to draw
//' the tree on. 
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param region_ids A V by 1 matrix with the region ids of each vertex
//' @param region_dvals A ndists by 1 matrix with the sizes of each regions 
//' @param verbose Whether or not to print out the inputted plan before
//' attemping to draw a tree. 
//'
//' @returns A list with the following 
//'     - `uncut_tree`: The spanning tree drawn on the region stored as a
//'     0-indexed directed edge adjacency graph.
//'     - `num_attempts`: The number of attempts it took to draw the tree.
//'     - `root`: The root vertex of the tree (0-indexed)
//' @export
List draw_a_tree_on_a_region(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int ndists, int num_regions, int num_districts,
    int region_id_to_draw_tree_on,
    double lower, double upper,
    arma::umat region_ids, arma::umat region_sizes,
    bool verbose
){
    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();

    
    double total_pop = sum(pop);
    double target = total_pop / ndists;

    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);

    bool split_district_only = region_sizes.max() == 1;

    // Create a plan object
    Plan * plan = new GraphPlan(
        region_ids.col(0), region_sizes.col(0), 
        ndists, total_pop, split_district_only,
        num_regions, num_districts, pop);


    if(verbose){
        Rprintf("Drawing Tree on Region %d of Plan: ", region_id_to_draw_tree_on);
        plan->Rprint();
    }

    // Create tree related stuff
    int root;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);



    // Counts the how many attempts it took to draw the tree
    bool successful_split_made = false;
    int num_attempts = 0;
    

    // Keep running until a tree is successfully drawn
    while(!successful_split_made){
        // try to draw a tree 
        successful_split_made = plan->draw_tree_on_region(map_params, region_id_to_draw_tree_on,
        ust, visited, ignore, root);
        num_attempts++;
    }


    List out = List::create(
        _["uncut_tree"] = ust,
        _["root"] = root,
        _["num_attempts"] = num_attempts
    );

    return out;
}




//' @title Split a multidistrict into two regions
//' 
//' Splits a multidistrict into Two New regions within population bounds
//'
//' Splits a multidistrict into two new valid regions by drawing spanning
//' trees uniformly at random and attempting to find an edge to cut until
//' a successful cut is made.
//'
//'
//' @inheritParams gsmc_plans
//' @inheritParams get_edge_to_cut
//'
//' @returns A list with the following
//' \itemize{
//'   \item{uncut_tree}{ - The spanning tree drawn stored as a 0-indexed directed
//'   adjacency list.}
//'   \item{root}{ - The 0-indexed root of the tree.}
//'   \item{num_attempts}{ - The number of attempts it took to draw the tree.}
//' }
List perform_a_valid_multidistrict_split(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int ndists, int num_regions, int num_districts,
    int region_id_to_split,
    double target, double lower, double upper,
    arma::umat region_ids, arma::umat region_sizes,
    int split_dval_min, int split_dval_max, bool split_district_only,
    bool verbose, int k_param
){
    if(split_dval_min > split_dval_max) throw Rcpp::exception("Split min must be less than split max!\n");
    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    double total_pop = sum(pop);

    Rprintf("It is %d\n", (int) split_district_only);

    // Create a plan object
    Plan * plan = new GraphPlan(
        region_ids.col(0), region_sizes.col(0), 
        ndists, total_pop, split_district_only,
        num_regions, num_districts, pop);


    if(verbose){
        Rprintf("Splitting Plan: ");
        plan->Rprint();
    }

    // Create tree related stuff
    int root;
    Tree ust = init_tree(V);
    Tree pre_split_ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);


    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan->region_ids(i) != region_id_to_split;
    }

    // Counts the number of split attempts
    bool successful_split_made = false;
    int try_counter {0};


    // splitting related params
    int new_region1_tree_root, new_region2_tree_root;
    int new_region1_dval, new_region2_dval;
    int new_region1_pop, new_region2_pop;

    int uncut_tree_root;

    // Keep running until done
    while(!successful_split_made){
        clear_tree(ust);
        // Get a tree
        int result = sample_sub_ust(g, ust, V, root, visited, ignore, pop, lower, upper, counties, cg);
        uncut_tree_root = root;
        // Return unsuccessful if tree not drawn
        if (result != 0){
            try_counter++;
            continue;
        }

        // copy uncut tree
        pre_split_ust = ust;

        // Try to make a cut
        int max_potential_d;

        split_dval_max = std::min(
            split_dval_max, 
            (int) plan->region_dvals(region_id_to_split) - 1
            );


        // try to get an edge to cut
        successful_split_made = get_edge_to_cut(ust, root,
                        k_param, split_dval_min, split_dval_max,
                        pop, plan->region_ids, 
                        region_id_to_split, plan->region_pops.at(region_id_to_split),
                        plan->region_dvals(region_id_to_split),
                        lower, upper, target,
                        new_region1_tree_root, new_region1_dval, new_region1_pop,
                        new_region2_tree_root, new_region2_dval, new_region2_pop);

        try_counter++;
        // increase the counter by 1
        if(verbose){
            Rcout << "Attempt " << try_counter << "\n";
        }

    }


    // make new region ids the split one and num_regions
    int new_region1_id = region_id_to_split; 
    int new_region2_id = plan->num_regions;

    // now update things with the new region ids 
    add_new_regions_to_plan_from_cut(
        ust, *plan, split_district_only,
        region_id_to_split, plan->num_regions,
        new_region1_tree_root, new_region1_dval,  new_region1_pop,
        new_region2_tree_root, new_region2_dval, new_region2_pop
    );


    if(verbose){
        plan->Rprint();
    }

    List out = List::create(
        _["num_attempts"] = try_counter,
        _["region_id_that_was_split"] = region_id_to_split,
        _["region_sizes"] = plan->region_dvals,
        _["partial_plan_labels"] = plan->region_ids,
        _["region_pops"] = plan->region_pops,
        _["num_regions"] = plan->num_regions,
        _["num_districts"] = plan->num_districts,
        _["uncut_tree"] = pre_split_ust,
        _["uncut_tree_root"] = uncut_tree_root,
        _["cut_tree"] = ust,
        _["new_region1_id"] = new_region1_id,
        _["new_region1_tree_root"] = new_region1_tree_root,
        _["new_region1_size"] = new_region1_dval,
        _["new_region1_pop"] = new_region1_pop,
        _["new_region2_id"] = new_region2_id,
        _["new_region2_tree_root"] = new_region2_tree_root,
        _["new_region2_size"] = new_region2_dval,
        _["new_region2_pop"] = new_region2_pop
    );

    return out;
}


// Does a preset number of merge split steps
List perform_merge_split_steps(
        List adj_list, const arma::uvec &counties, const arma::uvec &pop,
        int k_param, 
        double target, double lower, double upper,
        int ndists, int num_regions, int num_districts,
        arma::umat region_ids, arma::umat region_dvals,
        std::vector<int> region_pops,
        bool split_district_only, int num_merge_split_steps,
        bool verbose
){
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();
    double total_pop = sum(pop);

    auto dummy_region_ids = region_ids; auto dummy_region_dvals = region_dvals;

    // Create a plan object
    Plan *plan = new GraphPlan(region_ids.col(0), region_dvals.col(0), V, ndists, total_pop);
    Plan *new_plan = new GraphPlan(dummy_region_ids.col(0), dummy_region_dvals.col(0), V, ndists, total_pop);

    // fill in the plan
    plan->num_regions = num_regions;
    plan->num_districts = num_districts;
    plan->num_multidistricts = plan->num_regions - plan->num_districts;
    plan->region_ids = region_ids.col(0);
    plan->region_dvals = region_dvals.col(0);
    plan->region_pops = region_pops;

    // TODO FIX THIS but need to create this
    std::iota(plan->region_added_order.begin(), plan->region_added_order.end(), -1);


    if(verbose){
        plan->Rprint();
    }

    // Create tree related stuff
    int root;
    Tree ust = init_tree(V);
    Tree pre_split_ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);


    // now do merge split 
    int num_successes = run_merge_split_step_on_a_plan(
        map_params,
        split_district_only, "uniform",
        k_param,
        *plan, *new_plan, num_merge_split_steps
    );

    List out = List::create(
        _["region_dvals"] = plan->region_dvals,
        _["plan_vertex_ids"] = plan->region_ids,
        _["pops"] = plan->region_pops,
        _["num_regions"] = plan->num_regions,
        _["num_districts"] = plan->num_districts,
        _["num_success"] = num_successes
    );

    return out;
}

