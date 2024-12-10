/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Run SMC with MCMC merge-split steps mixed in
********************************************************/

#include "smc_and_mcmc.h"



// Takes a plan and reorders the regions 
// IT IS VERY IMPORTANT THAT THE TWO PLANS NOT POINT TO THE SAME arma::umat or everything breaks
// After reordering it copies that over to the dummy plan as well
void reorder_plan_by_oldest_split(
    Plan &plan, Plan &dummy_plan) {
    // Make dummy plan a shallow copy of the plan
    dummy_plan.shallow_copy(plan);

    // Recall that the region_added_order attribute is a vector that stores info
    // on the relative split order. So if entry i is greater than entry j that 
    // means region i was split more recently than region j

    // For example if you had a plan with 6 regions and a region_added_order vector set to
    // {4,1,7,3,9,2} then that means 
    //  - Region 0 was added 3rd
    //  - Region 1 was added 1st
    //  - Region 2 was added 5th
    //  - Region 3 was added 4rd
    //  - Region 4 was added 6th
    //  - Region 5 was added 2nd 

    // want to return [3, 0, 4, 2, 5, 1]
    // std::vector<int> rel_order = {4,1,7,3,9,2};

    // get the relative order vector of the n regions
    std::vector<int> rel_order(
        plan.region_added_order.begin(), 
        plan.region_added_order.begin() + plan.num_regions
        );


    // Create a vector of indices [0, 1, 2, ..., n-1]
    std::vector<int> indices(rel_order.size());
    std::iota(indices.begin(), indices.end(), 0);

    /*
    Sort the indices based on the values in rel_order so that means 
    indices[i] < indices[j] iff rel_order[i] < rel_order[j] so that means
    indices[0] is the old (ie not updated) label of the oldest split region, 
    In general indices[i] == k  means that the region with the old label of
    k will have a new id of i

    For our example with rel_order = {4,1,7,3,9,2} then that means the indices
    would be sorted to be {1, 5, 3, 0, 2, 4} so we know old region 1 is now zero,
    old region 5 is now 1, etc. Alternatively interpret as indices[i] == k means
    new region i was old region k
    */
    std::sort(indices.begin(), indices.end(), [&rel_order](int a, int b) {
        return rel_order[a] < rel_order[b];
    });


    /*
    Create a vector mapping old region id to the new ordered region id. So 
    `old_region_id_to_new_vec[i]` is the new region id of the region with old
    id i. In other words it means region i should be relabeled as old_region_id_to_new_vec[i]
    */
    std::vector<int> old_region_id_to_new_vec(rel_order.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        old_region_id_to_new_vec[indices[i]] = i;
    }


    // First we relabel all the region vertex ids
    for (size_t i = 0; i < plan.region_ids.n_elem; i++)
    {
        // Send region id i to old_region_id_to_new_vec[i]
        plan.region_ids(i) = old_region_id_to_new_vec.at(dummy_plan.region_ids(i));
    }


    // Relabel the remainder region if needed
    if(plan.remainder_region >= 0){
        plan.remainder_region = old_region_id_to_new_vec.at(dummy_plan.remainder_region);
    }
    


    // Now we reorder the region dvals and population 
    for (size_t i = 0; i < plan.num_regions; i++)
    {
        // Recall indices[i] == k means the old region with id k now has id i
        // so we want to set the value at the new region id `i` to the value it
        // had in the old one which is `indices[i]`
        int old_region_id = indices[i];
        int new_region_id = i;

        plan.region_dvals(new_region_id) = dummy_plan.region_dvals(old_region_id);
        plan.region_pops.at(new_region_id) = dummy_plan.region_pops.at(old_region_id);
        // Since the regions are in order of added reset the order added to just be 1,..., n
        plan.region_added_order.at(i) = i+1;
    }

    // reset the max region counter 
    plan.region_order_max = plan.num_regions + 6;
    
    // copy the dummy plan over
    dummy_plan.shallow_copy(plan);
}


void reorder_all_plans(
    RcppThread::ThreadPool &pool,
    std::vector<Plan> &plans_vec, 
    std::vector<Plan> &dummy_plans_vec){

    int M = (int) plans_vec.size();


    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // reorder every plan
        reorder_plan_by_oldest_split(plans_vec.at(i), dummy_plans_vec.at(i));
    });

    return;
    
}


// Different diagnostic levels
//      - level 0 - Does not capture any ancestry information or retain intermediate weights
//      - level 1 - Saves ancestry information, intermediate weights and the number of tries 
//      - level 2 - Captures the parent tries mat
//      - level 3 - Saves intermediate region dvals and plan ids


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
//' @param control Named list of additional parameters.
//' @param num_threads The number of threads the threadpool should use
//' @param verbosity What level of detail to print out while the algorithm is
//' running <ADD OPTIONS>
//' @export
List optimal_gsmc_with_merge_split_plans(
        int N, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        double target, double lower, double upper,
        int M, // M is Number of particles aka number of different plans
        arma::umat region_id_mat, arma::umat region_dvals_mat,
        List control, // control has pop temper, and k parameter value, and whether only district splits are allowed
        int verbosity, bool diagnostic_mode){

    // unpack control params
    int num_threads = (int) control["num_threads"];

    // set number of threads
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 1) num_threads = 0;

    
    // lags thing (copied from original smc code, don't understand what its doing)
    std::vector<int> lags = as<std::vector<int>>(control["lags"]);
    // k param values to use
    std::vector<int> k_params = as<std::vector<int>>(control["k_params"]);
    // Whether or not to only do district splits
    bool split_district_only = as<bool>(control["split_district_only"]);
    // This is a total_ms_steps by M vector where [s][i] is the number of 
    // successful merge splits performed for plan i on merge split round s
    std::vector<bool> merge_split_step_vec = as<std::vector<bool>>(control["merge_split_step_vec"]);

    double pop_temper = as<double>(control["pop_temper"]);

    int ms_steps_multiplier = as<int>(control["ms_steps_multiplier"]);

    // there are N-1 splits so for now just do it 
    int total_smc_steps = N-1;
    int total_ms_steps = std::count(merge_split_step_vec.begin(), merge_split_step_vec.end(), true);
    
    // total number of steps to run 
    int total_steps = total_smc_steps + total_ms_steps;

    // Merge split related diagnostics

    // For each merge split step this counts the number of attempts that were made
    std::vector<int> num_merge_split_attempts_vec(total_ms_steps, -1);

    // This is a M by total_ms_steps matrix where [i][s] is the number of 
    // successful merge splits performed for plan i on merge split round s
    arma::umat merge_split_successes_mat(M, total_ms_steps, arma::fill::zeros);

    // Make sure first merge split argument isn't true 
    if(merge_split_step_vec.at(0)){
        throw std::invalid_argument("The first entry of merge_split_step_vec cannot be true.");
    };
    
    umat ancestors(M, lags.size(), fill::zeros);

    // Create map level graph and county level multigraph
    Graph g = list_to_graph(adj_list);
    Multigraph cg = county_graph(g, counties);

    int V = g.size();
    double total_pop = sum(pop);

    // Loading Info
    if (verbosity >= 1) {
        Rcout.imbue(std::locale(""));
        Rcout << std::fixed << std::setprecision(4);
        if(!split_district_only){
            Rcout << "GENERALIZED SEQUENTIAL MONTE CARLO";
        }else{
            Rcout << "SEQUENTIAL MONTE CARLO";
        }
        Rcout << " WITH MERGE SPLIT\n";
        Rcout << "Sampling " << M << " " << V << "-unit ";
        Rcout << "maps with " << N << " districts and population between "
              << lower << " and " << upper << " using " << num_threads << " threads, "
              << total_ms_steps << " merge split steps, ";
        if(!split_district_only){
            Rcout << "and generalized region splits.\n";
        }else{
            Rcout << "and only performing 1-district splits.\n";
        }
        if (cg.size() > 1){
            Rcout << "Ensuring no more than " << N - 1 << " splits of the "
                  << cg.size() << " administrative units.\n";
        }
    }


    // Create copies of the matrices
    arma::umat dummy_region_id_mat = region_id_mat;
    arma::umat dummy_region_dvals_mat = region_dvals_mat;

    // Now create the vector of plans
    std::vector<Plan> plans_vec; plans_vec.reserve(M);
    std::vector<Plan> new_plans_vec; new_plans_vec.reserve(M);
    
    // Loop over each column of region_id_mat
    for (size_t i = 0; i < region_id_mat.n_cols; ++i) {
        // Create a Plan object and add it to the vector
        plans_vec.emplace_back(region_id_mat.col(i), region_dvals_mat.col(i), N, total_pop, split_district_only);
        new_plans_vec.emplace_back(dummy_region_id_mat.col(i), dummy_region_dvals_mat.col(i), N, total_pop, split_district_only);
    }


    // Define output variables that must always be created


    // TODO: All of these are actually total_steps instead of N-1 even though right now the 
    // merge split steps don't record anything 

    // This is N-1 by M where [i][j] is the index of the parent of particle j on step i
    // ie the index of the previous plan that was sampled and used to create particle j on step i
    // std::vector<std::vector<int>> parent_index_mat(total_steps, std::vector<int> (M, -1));

    // First column is dummy column
    arma::umat parent_index_mat(M, total_smc_steps, arma::fill::ones);

    

    // Entry [i][s] is the number of tries it took to form particle i on split s
    // arma::umat draw_tries_mat(M, total_steps, arma::fill::ones);

    Rcpp::IntegerMatrix draw_tries_mat(M, total_steps);
    draw_tries_mat.fill(0);

    

    // This is N-1 by M where [i][j] is the number of times particle j from the
    // previous round was sampled and unsuccessfully split on iteration i so this
    // does not count successful sample then split

    arma::umat parent_unsuccessful_tries_mat(M, total_smc_steps, arma::fill::zeros);

    //std::vector<std::vector<std::atomic<int>>> new_parent_unsuccessful_tries_mat(total_steps, std::vector<std::atomic<int>> (M, -1));

    arma::umat da_parent_unsuccessful_tries_mat(M, total_smc_steps, arma::fill::none);

    std::vector<std::atomic<uint>> new_parent_unsuccessful_tries_mat(M);

    // Initialize the vector elements
    for (size_t i = 0; i < M; ++i) {
        new_parent_unsuccessful_tries_mat.at(i).store(static_cast<uint>(0)); // Explicitly store values
    }

    // Tracks the number of unique parent vectors sampled to create the next round
    std::vector<int> nunique_parents_vec(total_steps, -1);


    // These are M by number of smc steps 
    // entry [i][s] is the normalized weight of particle i AFTER split s
    arma::dmat log_incremental_weights_mat(M, total_smc_steps, arma::fill::none);


    // Tracks the acceptance rate - total number of tries over M - for each round
    std::vector<double> acceptance_rates(total_steps, -1.0);

    // Tracks the effective sample size for the weights of each round
    std::vector<double> n_eff(total_smc_steps, -1.0);


    // Declare variables whose size will depend on whether or not we're in
    // diagnostic mode or not
    arma::ucube plan_region_ids_mat;

    
    int plan_dval_list_size = (diagnostic_mode & !split_district_only) ? total_steps-1 : 0;
    Rcpp::List region_dvals_mat_list(plan_dval_list_size);

    // If diagnostic mode track vertex region ids from every round
    if(diagnostic_mode){
        // This is a V by M by N-2 cube (3d array) 
        // plan_region_ids_mat[,j, s] is the jth plan of split s
        // plan_region_ids_mat[v ,j, s] is the region id of the jth plan of split s
        plan_region_ids_mat.set_size(V, M, total_steps-1);

        // if doing generalized region splits then allocate the matrices for the first N-2 splits
        if(!split_district_only){
            for (size_t i = 0; i < total_steps-1; i++)
            {
                region_dvals_mat_list(i) = arma::umat(i+2, M, arma::fill::none);
            }
        }
    }

    

    // Start off all the unnormalized weights at 1
    std::vector<double> unnormalized_sampling_weights(M, 1.0);
    
    // Create a threadpool
    RcppThread::ThreadPool pool(num_threads);

    std::string bar_fmt = "Split [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} | ETA{cli::pb_eta}";
    RObject bar = cli_progress_bar(total_steps, cli_config(false, bar_fmt.c_str()));


    // counts the number of smc steps
    int smc_step_num = 0;
    int merge_split_step_num = 0;

    // Now for each run through split the map
    try {
    for(int step_num=0; step_num < total_steps; step_num++){
        if(verbosity > 1){
            if(merge_split_step_vec[step_num]){
                Rprintf("Iteration %d: Merge Split Step %d \n", step_num+1, merge_split_step_num + 1);
            }else{
                Rprintf("Iteration %d: SMC Step %d \n",  step_num+1, smc_step_num + 1);
            }
        }

        // Check what step type
        if(!merge_split_step_vec[step_num]){

        // split the map
        generalized_split_maps(
            g, counties, cg, pop,
            plans_vec, new_plans_vec,
            parent_index_mat.col(smc_step_num),
            unnormalized_sampling_weights,
            draw_tries_mat.column(step_num),
            parent_unsuccessful_tries_mat.col(smc_step_num),
            new_parent_unsuccessful_tries_mat,
            acceptance_rates.at(step_num),
            nunique_parents_vec.at(smc_step_num),
            ancestors, lags,
            lower, upper, target,
            k_params.at(smc_step_num), split_district_only,
            pool,
            verbosity, diagnostic_mode ? 3 : 0
        );

        // compute log incremental weights and sampling weights for next round
        get_all_plans_log_gsmc_weights(
            pool,
            g,
            plans_vec,
            split_district_only,
            log_incremental_weights_mat.col(smc_step_num),
            unnormalized_sampling_weights,
            target,
            pop_temper
        );

        // compute effective sample size
        n_eff.at(smc_step_num) = compute_n_eff(log_incremental_weights_mat.col(smc_step_num));

        if(step_num == 0){
            // For the first ancestor one make every ancestor themselves
            // regspace<ucolvec>(0,  9);
            std::iota(parent_index_mat.col(0).begin(), parent_index_mat.col(0).end(), 0);
        }

        // only increase if we have smc steps left else it will cause index issues
        // with merge split
        if(smc_step_num < total_smc_steps-1){
            smc_step_num++;
        }

        }else if(merge_split_step_vec[step_num]){ // check if its a merge split step
             // Rprintf("%d!\n", step_num);
            // run merge split 
            // Set the number of steps to run at 1 over previous stage acceptance rate
            // int nsteps_to_run = std::ceil(1/std::pow(acceptance_rates.at(step_num-1),1.5)); // * std::max(merge_split_step_num,1);
            int nsteps_to_run = ms_steps_multiplier * std::ceil(1/acceptance_rates.at(step_num-1)); // * std::max(merge_split_step_num,1);
            num_merge_split_attempts_vec.at(merge_split_step_num) = nsteps_to_run;

            run_merge_split_step_on_all_plans(
                pool,
                g, counties, cg, pop,
                plans_vec, new_plans_vec,
                split_district_only, k_params.at(smc_step_num),
                nsteps_to_run,
                lower, upper, target,
                merge_split_successes_mat.col(merge_split_step_num)
            );

            // set the acceptance rate 
            int total_ms_successes = arma::sum(merge_split_successes_mat.col(merge_split_step_num));

            int total_ms_attempts = M * nsteps_to_run;

            acceptance_rates.at(step_num) = total_ms_successes / static_cast<double>(total_ms_attempts);

            if (verbosity >= 3){
                Rprintf("Ran %d Merge Split Attempts: %d Successes out of %d attempts. Acceptance Rate: %.2f\n", 
                    nsteps_to_run,  
                    total_ms_successes,total_ms_attempts, 100.0*acceptance_rates.at(step_num));
            }


            // draw_tries_mat.column(step_num) = arma::ucolvec(M, arma::fill::value(nsteps_to_run));

            merge_split_step_num++;
        }

        if (verbosity == 1 && CLI_SHOULD_TICK){
            cli_progress_set(bar, step_num);
        }
        Rcpp::checkUserInterrupt();

        

        // Now update the diagnostic info if needed, region labels, dval column of the matrix
        if(diagnostic_mode){ // record if in diagnostic mode and generalized splits

            // reorder the plans by oldest split
            reorder_all_plans(pool, plans_vec, new_plans_vec);
            
            arma::ucolvec data_vec(M);
            for (size_t i = 0; i < M; ++i) {
                // fetch the count
                data_vec(i) = static_cast<arma::uword>(new_parent_unsuccessful_tries_mat[i].load());
                // set it to 
                new_parent_unsuccessful_tries_mat.at(i).store(static_cast<uint>(0));
            }
            da_parent_unsuccessful_tries_mat.col(smc_step_num) = data_vec;

            if(step_num < total_steps -1){
                // make the step_num slice the region ids matrix
                plan_region_ids_mat.slice(step_num) = region_id_mat;
                if(!split_district_only){
                // Get the first number of region rows in the dvals mat since the rest
                // are zero
                region_dvals_mat_list(step_num) = region_dvals_mat.submat(
                    0, 0,
                    plans_vec.at(0).num_regions-1,
                    region_dvals_mat.n_cols-1
                );
                }
            }


        }

    }
    } catch (Rcpp::internal::InterruptedException e) {
        cli_progress_done(bar);
        return R_NilValue;
    }
    

    cli_progress_done(bar);

    // Return results
    List out = List::create(
        _["plans_mat"] = region_id_mat,
        _["parent_index"] = parent_index_mat,
        _["region_ids_mat_list"] = plan_region_ids_mat,
        _["region_dvals_mat_list"] = region_dvals_mat_list,
        _["log_incremental_weights_mat"] = log_incremental_weights_mat,
        _["draw_tries_mat"] = draw_tries_mat,
        _["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat,
        _["acceptance_rates"] = acceptance_rates,
        _["nunique_parent_indices"] = nunique_parents_vec,
        _["ancestors"] = ancestors,
        _["step_n_eff"] = n_eff,
        _["est_k"] = k_params,
        _["merge_split_steps"] = merge_split_step_vec,
        _["merge_split_attempt_counts"] = num_merge_split_attempts_vec,
        _["merge_split_success_mat"] = merge_split_successes_mat,
        _["Daco"] = da_parent_unsuccessful_tries_mat
    );

    return out;

}
