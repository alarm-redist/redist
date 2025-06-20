/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Exposes various aspects of the splitting procedure
*   to R to allow for manual splitting of plans.
********************************************************/

#include "manual_splitting.h"


std::pair<int, int> TEMP_get_potential_region_size_for_loop_bounds(
    const int total_region_size,
    const int min_potential_cut_size, const int max_potential_cut_size
){
    
    // if biggest possible cut size is leq half just return the same bounds
    if(max_potential_cut_size <= total_region_size/2){
        return std::pair<int, int>(min_potential_cut_size, max_potential_cut_size);
    }else if(min_potential_cut_size > total_region_size/2){
        // else if smallest possible size if more than half just subtract 
        // total_region_size and flip 
        return std::pair<int, int>(
            total_region_size - max_potential_cut_size,
            total_region_size - min_potential_cut_size
            );
    }else{ // else must be true that 
    // min_potential_cut_size <= total_region_size/2 < max_potential_cut_size

    return std::pair<int, int>(
            std::min(min_potential_cut_size, total_region_size - max_potential_cut_size),
            total_region_size/2
            );

    }
    
}

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
//' @param region_sizes A ndists by 1 matrix with the sizes of each regions 
//' @param verbose Whether or not to print out the inputted plan before
//' attemping to draw a tree. 
//'
//' @returns A list with the following 
//'     - `uncut_tree`: The spanning tree drawn on the region stored as a
//'     0-indexed directed edge adjacency graph.
//'     - `num_attempts`: The number of attempts it took to draw the tree.
//'     - `root`: The root vertex of the tree (0-indexed)
//'     - `pop_below`: The population below each vertex in `uncut_tree` ie
//'     the population induced by removing the edge terminating in that vertex
//'     - `uncut_tree_vertex_parents`: The parents of each vertex in the tree
//' @export
List draw_a_tree_on_a_region(
    List adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int ndists, int num_regions, int num_districts,
    int region_id_to_draw_tree_on,
    double lower, double upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    bool verbose
){
    // unpack control params
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);
    int V = g.size();

    
    double total_pop = sum(pop);
    double target = total_pop / ndists;

    MapParams map_params(
        adj_list, counties, pop, 
        ndists, ndists, std::vector<int>{},
        lower, target, upper);


    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);
    // Create a plan object via size 1 ensemble
    RcppThread::ThreadPool pool(0);
    PlanEnsemble plan_ensemble(
        map_params, num_regions,
        1, SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );


    if(verbose){
        Rprintf("Drawing Tree on Region %d of Plan: ", region_id_to_draw_tree_on);
        plan_ensemble.plan_ptr_vec[0]->Rprint();
    }

    // Create tree related stuff
    int root;
    Tree ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V, -2);

    // RNGState rng_state();
    RNGState rng_state;


    // Counts the how many attempts it took to draw the tree
    bool successful_split_made = false;
    int num_attempts = 0;
    

    // Keep running until a tree is successfully drawn
    while(!successful_split_made){
        // try to draw a tree 

        // Mark it as ignore if its not in the region to split
        for (int i = 0; i < V; i++){
            ignore[i] = plan_ensemble.plan_ptr_vec[0]->region_ids[i] != region_id_to_draw_tree_on;
        }

        // clear the tree
        clear_tree(ust);
        // Get a uniform spanning tree drawn on that region
        int result = sample_sub_ust(map_params.g, ust, 
            V, root, visited, ignore, 
            map_params.pop, 
            map_params.lower, map_params.upper, 
            map_params.counties, map_params.cg,
            rng_state);

        // result == 0 means it was successful
        successful_split_made = result == 0;

        num_attempts++;
    }

    // computes population below each vtx and parent of each vertex
    tree_vertex_parents.at(root) = -1;
    tree_pop(ust, root, pop, pop_below, tree_vertex_parents);


    List out = List::create(
        _["uncut_tree"] = ust,
        _["root"] = root,
        _["num_attempts"] = num_attempts,
        _["pop_below"] = pop_below,
        _["uncut_tree_vertex_parents"] = tree_vertex_parents
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
//' @inheritParams run_redist_gsmc
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
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    int split_dval_min, int split_dval_max, bool split_district_only,
    bool verbose, int k_param
){
    throw Rcpp::exception("Not working right now!");
    if(split_dval_min > split_dval_max) throw Rcpp::exception("Split min must be less than split max!\n");
    MapParams map_params(
        adj_list, counties, pop,
        ndists, ndists, std::vector<int>{},
        lower, target, upper);
    // unpack control params
    int V = map_params.V;
    

    Rprintf("It is %d\n", (int) split_district_only);

    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);
    // Create a plan object via size 1 ensemble
    RcppThread::ThreadPool pool(0);
    PlanEnsemble plan_ensemble(
        map_params, num_regions,
        1, SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );

    // Create tree splitter 
    TreeSplitter * tree_splitter = new ExperimentalSplitter(map_params.V, .0001, map_params.target);


    if(verbose){
        Rprintf("Splitting Plan: ");
        plan_ensemble.plan_ptr_vec[0]->Rprint();
    }

    auto splitting_schedule_ptr = std::make_unique<PureMSSplittingSchedule>(ndists, ndists, std::vector<int>{1});

    // Create tree related stuff
    int uncut_tree_root;
    USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
    Tree pre_split_ust = init_tree(V);
    std::vector<bool> visited(V);
    std::vector<bool> ignore(V, false);
    std::vector<int> pop_below(V, 0);
    std::vector<int> tree_vertex_parents(V, -2);


    // Mark it as ignore if its not in the region to split
    for (int i = 0; i < V; i++){
        ignore[i] = plan_ensemble.plan_ptr_vec[0]->region_ids[i] != region_id_to_split;
    }

    // Counts the number of split attempts
    bool successful_split_made = false;
    int try_counter = 1;

    int region_to_split_size = static_cast<int>(plan_ensemble.plan_ptr_vec[0]->region_ids[region_id_to_split]);

    std::pair<int, int> loop_bounds = TEMP_get_potential_region_size_for_loop_bounds(
    region_to_split_size,
    split_dval_min, split_dval_max
    );

    std::vector<int> smaller_cut_sizes_to_try;
    for (int i = loop_bounds.first; i <= loop_bounds.second; i++)
    {
        smaller_cut_sizes_to_try.push_back(i);
    }


    RNGState rng_state;

    // make new region ids the split one and num_regions
    int new_region1_id = region_id_to_split; 
    int new_region2_id = plan_ensemble.plan_ptr_vec[0]->num_regions;

    EdgeCut cut_edge;
    // Keep running until done
    while(!successful_split_made){
        if(verbose){
            Rcout << "Attempt " << try_counter << "\n";
        }
        clear_tree(ust_sampler.ust);
 
        // Try to draw a tree on region
        // bool tree_drawn = plan_ensemble.plan_ptr_vec[0]->draw_tree_on_region(map_params, region_id_to_split,
        //     ust_sampler.ust, visited, ignore, uncut_tree_root, rng_state);
        bool tree_drawn = false;

        // Try again and increase counter if tree not drawn
        if (!tree_drawn){
            try_counter++;
            continue;
        }

        split_dval_max = std::min(
            split_dval_max, 
            (int) plan_ensemble.plan_ptr_vec[0]->region_sizes[region_id_to_split] - 1
            );



        // Now try to select an edge to cut
        std::pair<bool, EdgeCut> edge_search_result = tree_splitter->attempt_to_find_edge_to_cut(
            map_params, rng_state, 
            ust_sampler.ust, uncut_tree_root,
            pop_below, visited, 
            plan_ensemble.plan_ptr_vec[0]->region_pops[region_id_to_split],
            plan_ensemble.plan_ptr_vec[0]->region_sizes[region_id_to_split], 
            split_dval_min, split_dval_max, 
            smaller_cut_sizes_to_try //,
            //bool save_selection_prob = false
        );


        // Try again and increase counter if no good edge found
        if (!std::get<0>(edge_search_result)){
            try_counter++;
            continue;
        }

        // copy uncut tree before splitting 
        pre_split_ust = ust_sampler.ust;

        // If successful extract the edge cut info
        cut_edge = std::get<1>(edge_search_result);
        // Now erase the cut edge in the tree
        erase_tree_edge(ust_sampler.ust, cut_edge);


        // now update the region level information from the edge cut
        plan_ensemble.plan_ptr_vec[0]->update_region_info_from_cut(
            cut_edge,
            region_id_to_split, new_region2_id, 
            true
        );

        // Now update the vertex level information
        plan_ensemble.plan_ptr_vec[0]->update_vertex_and_plan_specific_info_from_cut(
            *tree_splitter,
            ust_sampler, cut_edge, 
            region_id_to_split, new_region2_id,
            true
        );
    }


    if(verbose){
        plan_ensemble.plan_ptr_vec[0]->Rprint();
    }

    // splitting related params
    int new_region1_tree_root, new_region2_tree_root;
    int new_region1_dval, new_region2_dval;
    int new_region1_pop, new_region2_pop;


    cut_edge.get_split_regions_info(
        new_region1_tree_root, new_region1_dval, new_region1_pop,
        new_region2_tree_root, new_region2_dval, new_region2_pop
    );

    List out = List::create(
        _["num_attempts"] = try_counter,
        _["region_id_that_was_split"] = region_id_to_split,
        _["region_sizes"] = plan_ensemble.flattened_all_region_sizes,
        _["partial_plan_labels"] = plan_ensemble.flattened_all_plans,
        _["region_pops"] = plan_ensemble.flattened_all_region_pops,
        _["num_regions"] = plan_ensemble.plan_ptr_vec[0]->num_regions,
        _["num_districts"] = plan_ensemble.plan_ptr_vec[0]->get_num_district_and_multidistricts().first,
        _["uncut_tree"] = pre_split_ust,
        _["uncut_tree_root"] = uncut_tree_root,
        _["cut_tree"] = ust_sampler.ust,
        _["pop_below"] = pop_below,
        _["uncut_tree_vertex_parents"] = tree_vertex_parents,
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


// // Does a preset number of merge split steps
// List perform_merge_split_steps(
//         List adj_list, const arma::uvec &counties, const arma::uvec &pop,
//         int k_param, 
//         double target, double lower, double upper,
//         int ndists, int num_regions, int num_districts,
//         arma::umat region_ids, arma::umat region_sizes,
//         std::vector<int> region_pops,
//         bool split_district_only, int num_merge_split_steps,
//         bool verbose
// ){
//     throw Rcpp::exception("Not support rn!");
//     double rho = 1; bool is_final = false;
//     MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
//     // unpack control params
//     Graph g = list_to_graph(adj_list);
//     Multigraph cg = county_graph(g, counties);
//     int V = g.size();
//     ScoringFunction scoring_function(map_params, adj_list, 0);

//     auto dummy_region_ids = region_ids; auto dummy_region_dvals = region_sizes;

//     // Create a plan object
//     std::unique_ptr<Plan> plan = std::make_unique<GraphPlan>(region_ids.col(0), region_sizes.col(0), ndists, num_regions, pop, split_district_only);
//     std::unique_ptr<Plan> new_plan = std::make_unique<GraphPlan>(dummy_region_ids.col(0), dummy_region_dvals.col(0), ndists, num_regions, pop, split_district_only);

//     // create splitter
//     TreeSplitter *tree_splitter = new NaiveTopKSplitter(map_params.V, k_param);

//     // fill in the plan
//     plan->num_regions = num_regions;
//     // plan->region_ids = std::vector< region_ids.col(0);
//     // plan->region_sizes = region_sizes.col(0);
//     plan->region_pops = region_pops;

//     // TODO FIX THIS but need to create this
//     std::iota(plan->region_added_order.begin(), plan->region_added_order.end(), -1);


//     if(verbose){
//         plan->Rprint();
//     }

//     // Create tree related stuff
//     int root;
//     Tree ust = init_tree(V);
//     Tree pre_split_ust = init_tree(V);
//     std::vector<bool> visited(V);
//     std::vector<bool> ignore(V, false);

//     Rcpp::List fake_control;
//     SplittingSizeScheduleType splitting_type = get_splitting_size_regime("FAKE");

//     auto splitting_schedule_ptr = get_splitting_schedule(
//         1, ndists, splitting_type, fake_control
//     );

//     int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
//     RNGState rng_state(global_rng_seed);
//     SamplingSpace sampling_space = get_sampling_space("graph_space");
//     USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
//     std::vector<int> tree_sizes(ndists, 0);
//     std::vector<int> succesful_tree_sizes(ndists, 0);

//     // now do merge split 
//     int num_successes = run_merge_split_steps(
//         map_params, *splitting_schedule_ptr, scoring_function,
//         rng_state, sampling_space,
//         *plan, *new_plan, 
//         ust_sampler, *tree_splitter,
//         "uniform", 
//         rho, is_final, 
//         num_merge_split_steps,
//         tree_sizes, succesful_tree_sizes
//     );

//     List out = List::create(
//         _["region_sizes"] = plan->region_sizes,
//         _["plan_vertex_ids"] = plan->region_ids,
//         _["pops"] = plan->region_pops,
//         _["num_regions"] = plan->num_regions,
//         _["num_districts"] = plan->get_num_district_and_multidistricts().first,
//         _["num_success"] = num_successes
//     );

//     return out;
// }


// TODO: Add support for multimember districts 
// Draws num_trees number of trees on a region
List draw_trees_on_a_region(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int const ndists,
    int const region_id_to_draw_tree_on, int const region_size,
    double const lower, double const target, double const upper,
    arma::uvec const &region_ids, 
    int const num_tree, int num_threads,
    bool const verbose
){

    // Create adj params 
    MapParams map_params(adj_list, counties, pop, 
        ndists, ndists, std::vector<int>{},
        lower, target, upper);
    // count how many times we had to call sample_sub_ust

    // create thread pool
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    RcppThread::ThreadPool pool(num_threads);


    // create list of trees to return 
    std::vector<std::vector<Graph>> thread_undirected_trees(num_threads == 0 ? 1 : num_threads); 
    std::vector<int> thread_attempts(num_threads == 0 ? 1 : num_threads,0);
    std::atomic<int> thread_id_counter{0};


    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    int num_rng_states = num_threads > 0 ? num_threads : 1;
    std::vector<RNGState> rng_states;rng_states.reserve(num_rng_states);
    for (size_t i = 1; i <= num_rng_states; i++)
    {
        // same seed with i*3 long_jumps for state
        rng_states.emplace_back(global_rng_seed, i*3);
    }

    RcppThread::ProgressBar bar(num_tree, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, num_tree, [&] (int ree) {

        static thread_local int thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
        // Stuff for drawing tree
        int root;
        static thread_local Tree ust = init_tree(map_params.V); 
        static thread_local std::vector<bool> visited(map_params.V);
        static thread_local std::vector<bool> ignore(map_params.V, false);
        // Tree ust = init_tree(map_params.V); 
        // std::vector<bool> visited(map_params.V);
        // std::vector<bool> ignore(map_params.V, false);

        // reset result 
        int result = 1;
        while(result != 0){
            // clear tree 
            clear_tree(ust);
            for (size_t i = 0; i < map_params.V; i++)
            {
                ignore[i] = region_ids(i) != region_id_to_draw_tree_on;
            }
            
            // sample until successful
            result = sample_sub_ust(
                map_params.g, ust, map_params.V, root,
                visited, ignore, map_params.pop,
                lower, upper,
                map_params.counties, map_params.cg,
                rng_states[thread_id]
            );
            ++thread_attempts[thread_id];
        }

        // go through the tree from the root and add the backwards edge and sort 
        std::queue<std::pair<int,int>> vertex_queue;
        // add roots children to queue 
        for(auto const &child_vertex: ust[root]){
            vertex_queue.push({child_vertex, root});
        }
        // sort the children 
        std::sort(ust[root].begin(), ust[root].end());

        // update all the children
        while(!vertex_queue.empty()){
            // get and remove head of queue 
            auto queue_pair = vertex_queue.front();
            int vertex = queue_pair.first;
            int parent_vertex = queue_pair.second;
            vertex_queue.pop();
            // add children to the queue 
            for(auto const &child_vertex: ust[vertex]){
                // add children to queue
                vertex_queue.push({child_vertex, vertex});
            }
            // add the edge from vertex to parent 
            ust[vertex].push_back(parent_vertex);
            // now sort edges
            std::sort(ust[vertex].begin(), ust[vertex].end());
        }
        // REprintf("about to copy! %d\n", thread_id);
        thread_undirected_trees[thread_id].push_back(ust);
        // REprintf("Copied! %d\n", thread_id);
        ++bar;
    });

    pool.wait();

    std::vector<Graph> undirected_trees; undirected_trees.reserve(num_tree);
    int num_attempts = 0;

    for (size_t i = 0; i < num_threads; i++)
    {
        // move don't copy 
        // https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
        undirected_trees.insert(
            undirected_trees.end(),
            std::make_move_iterator(thread_undirected_trees[i].begin()),
            std::make_move_iterator(thread_undirected_trees[i].end())
            );
        num_attempts += thread_attempts[i];
    }
    
    

    List out = List::create(
        _["trees_list"] = undirected_trees,
        _["num_attempts"] = num_attempts
    );

    return out;
}






// Draws num_plans number of plans on a region
// if unsuccessful then just returns the unsplit plan
List attempt_splits_on_a_region(
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    int const ndists, int const init_num_regions,
    int const region_id_to_split,
    double const lower, double const target, double const upper,
    Rcpp::IntegerMatrix const &region_ids, 
    Rcpp::IntegerMatrix const &region_sizes,
    std::string const &splitting_schedule_str, int const k_param,
    int const num_plans, int num_threads,
    bool const verbose
){

    // Create adj params 
    MapParams map_params(
        adj_list, counties, pop, 
        ndists, ndists, std::vector<int>{},
        lower, target, upper);
    // count how many times we had to call sample_sub_ust
    Rcpp::List control;

    // get splitting schedule 
    SplittingSizeScheduleType splitting_schedule_type = get_splitting_size_regime(splitting_schedule_str);
    auto splitting_schedule = get_splitting_schedule(
        1, ndists, ndists, std::vector<int>{},
        splitting_schedule_type, control
    );
    
    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    int num_rng_states = num_threads > 0 ? num_threads : 1;
    std::vector<RNGState> rng_states;rng_states.reserve(num_rng_states);
    for (size_t i = 1; i <= num_rng_states; i++)
    {
        // same seed with i*3 long_jumps for state
        rng_states.emplace_back(global_rng_seed, i*3);
    }
    
    // create thread pool
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    RcppThread::ThreadPool pool(num_threads);

    // create the plan
    PlanEnsemble plan_ensemble(
        map_params, init_num_regions,
        1, SamplingSpace::GraphSpace,
        region_ids, 
        region_sizes, rng_states,
        pool 
    );


    splitting_schedule->set_potential_cut_sizes_for_each_valid_size(
        0, plan_ensemble.plan_ptr_vec[0]->num_regions-1
    );


    PlanEnsemble thread_plan_ensemble(
        map_params, arma::sum(pop),
        num_threads, SamplingSpace::GraphSpace,
        pool
    );

    // create the splitter
    NaiveTopKSplitter tree_splitter(map_params.V, k_param);


    // Create the vector of plans to return
    Rcpp::IntegerMatrix saved_plans_mat(map_params.V, num_plans);
    Rcpp::IntegerMatrix saved_region_sizes_mat(ndists, num_plans);


    // create list of trees to return 
    std::vector<std::vector<Graph>> thread_undirected_trees(num_threads == 0 ? 1 : num_threads); 
    std::vector<int> thread_attempts(num_threads == 0 ? 1 : num_threads,0);
    std::atomic<int> thread_id_counter{0};




    std::vector<bool> successful_update(num_plans);

    RcppThread::ProgressBar bar(num_plans, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, num_plans, [&] (int i) {

        static thread_local int thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
        // Stuff for drawing tree
        static thread_local USTSampler ust_sampler(map_params, *splitting_schedule);
        // copy the plan again
        thread_plan_ensemble.plan_ptr_vec[thread_id]->shallow_copy(*plan_ensemble.plan_ptr_vec[0]);


        // keep trying until a tree is drawn
        // ie ignore cases where algorithm fails bc of randomness
        bool tree_successfully_drawn = ust_sampler.attempt_to_draw_tree_on_region(
            rng_states[thread_id], *thread_plan_ensemble.plan_ptr_vec[thread_id], region_id_to_split);
        ++thread_attempts[thread_id];

        while(!tree_successfully_drawn){
            tree_successfully_drawn = ust_sampler.attempt_to_draw_tree_on_region(
                rng_states[thread_id], *thread_plan_ensemble.plan_ptr_vec[thread_id], region_id_to_split);
            ++thread_attempts[thread_id];
        }

        // now draw a tree 
        auto edge_search_result = ust_sampler.try_to_sample_splittable_tree(
            rng_states[thread_id], tree_splitter,
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_pops[region_id_to_split], 
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_sizes[region_id_to_split],
            false
        );

        // if successful then update
        if(std::get<0>(edge_search_result)){
            // now split that region we found on the old one
            thread_plan_ensemble.plan_ptr_vec[thread_id]->update_from_successful_split(
                tree_splitter,
                ust_sampler, std::get<1>(edge_search_result),
                region_id_to_split, thread_plan_ensemble.plan_ptr_vec[thread_id]->num_regions, 
                true
            );
        }
        successful_update[i] = std::get<0>(edge_search_result);
        // Copy the plan into the matrix 
        std::copy(
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_ids.begin(), 
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_ids.end(), 
            saved_plans_mat.column(i).begin() // Start of column in Rcpp::IntegerMatrix
        );
        std::copy(
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_sizes.begin(), 
            thread_plan_ensemble.plan_ptr_vec[thread_id]->region_sizes.end(), 
            saved_region_sizes_mat.column(i).begin() // Start of column in Rcpp::IntegerMatrix
        );

    });

    pool.wait();


    int num_attempts = 0;

    for (size_t i = 0; i < num_threads; i++)
    {
        num_attempts += thread_attempts[i];
    }
    
    

    List out = List::create(
        _["plans_mat"] = saved_plans_mat,
        _["sizes_mat"] = saved_region_sizes_mat,
        _["successful_search"] = successful_update,
        _["num_attempts"] = num_attempts
    );

    return out;
}