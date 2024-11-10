/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Run SMC with MCMC merge-split steps mixed in
********************************************************/

#include "smc_and_mcmc.h"


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
        List control, // control has pop temper, and k parameter value, and whether only district splits are allowed
        int num_threads, int verbosity, bool diagnostic_mode){

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

    // there are N-1 splits so for now just do it 
    int total_smc_steps = N-1;
    int total_ms_steps = std::count(merge_split_step_vec.begin(), merge_split_step_vec.end(), true);
    
    // total number of steps to run 
    int total_steps = total_smc_steps + total_ms_steps;

    // Merge split related diagnostics

    // For each merge split step this counts the number of attempts that were made
    std::vector<int> num_merge_split_attempts_vec(total_ms_steps, -1);

    // This is a total_ms_steps by M vector where [s][i] is the number of 
    // successful merge splits performed for plan i on merge split round s
    std::vector<std::vector<int>> merge_split_successes_mat(total_ms_steps,
        std::vector<int> (M, -1)
    );


    // return List::create(
    //     _["merge_split_steps"] = merge_split_step_vec,
    //     _["count"] = cnnt
    // );

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


    std::vector<Plan> plans_vec(M, Plan(V, N, total_pop, split_district_only));
    std::vector<Plan> new_plans_vec(M, Plan(V, N, total_pop, split_district_only)); // New plans


    // Define output variables that must always be created


    // TODO: All of these are actually total_steps instead of N-1 even though right now the 
    // merge split steps don't record anything 

    // This is N-1 by M where [i][j] is the index of the parent of particle j on step i
    // ie the index of the previous plan that was sampled and used to create particle j on step i
    std::vector<std::vector<int>> parent_index_mat(total_steps, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the index of the original (first) ancestor of particle j on step i
    std::vector<std::vector<int>> original_ancestor_mat(total_steps, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the number of tries it took to form particle j on iteration i
    // Inclusive of the final step. ie if succeeds in one try it would be 1
    std::vector<std::vector<int>> draw_tries_mat(total_steps, std::vector<int> (M, -1));

    // This is N-1 by M where [i][j] is the number of times particle j from the
    // previous round was sampled and unsuccessfully split on iteration i so this
    // does not count successful sample then split
    std::vector<std::vector<int>> parent_unsuccessful_tries_mat(total_steps, std::vector<int> (M, 0));

    // Tracks the number of unique parent vectors sampled to create the next round
    std::vector<int> nunique_parents_vec(total_steps, -1);

    // Tracks the number of unique ancestors left at each step
    std::vector<int> nunique_original_ancestors_vec(total_steps, -1);


    // This is N-1 by M where [i][j] is the log incremental weight of particle j on step i
    std::vector<std::vector<double>> log_incremental_weights_mat(total_steps, std::vector<double> (M, -1.0));

    // This is N-1 by M where [i][j] is the normalized weight of particle j on step i
    std::vector<std::vector<double>> normalized_weights_mat(total_steps, std::vector<double> (M, -1.0));


    // Tracks the acceptance rate - total number of tries over M - for each round
    std::vector<double> acceptance_rates(total_steps, -1.0);

    // Tracks the effective sample size for the weights of each round
    std::vector<double> n_eff(total_steps, -1.0);


    // Declare variables whose size will depend on whether or not we're in
    // diagnostic mode or not
    std::vector<std::vector<std::vector<int>>> plan_region_ids_mat;
    std::vector<std::vector<std::vector<int>>> plan_d_vals_mat;
    std::vector<std::vector<std::vector<int>>> plan_region_order_added_mat;


    // If diagnostic mode track stuff from every round
    if(diagnostic_mode){
        // Create info tracking we will pass out at the end
        // This is N-1 by M by V where for each n=1,...,N-1 and m=1,...M it maps vertices to integer id
        plan_region_ids_mat.resize(total_steps, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));

        // reserve space for total_steps-1 elements if doing generalized region splits
        if(!split_district_only){
            plan_d_vals_mat.reserve(total_steps-1);
        }
        
        // reserve space for total_steps number of elements
        plan_region_order_added_mat.reserve(total_steps);


        int num_regions = 2;
        for (size_t i = 0; i < total_steps; i++){
            // Rprintf("Iteration %zu and Number of Regions is %d\n", i, num_regions);
            // If doing generalized splits then make num regions by M vector to 
            // track dvals whenever its less than N regions
            if(!split_district_only && num_regions < N){
                plan_d_vals_mat.push_back(
                    std::vector<std::vector<int>>(M, std::vector<int>(num_regions, -1))
                );
            }

            plan_region_order_added_mat.push_back(
                std::vector<std::vector<int>>(M, std::vector<int>(num_regions, -1))
            );

            // if its not an MCMC step increase the number of regions
            if(!merge_split_step_vec.at(i)){
                num_regions++;
            }

        }
        
    }else{ // else only track for final round
        // Create info tracking we will pass out at the end
        // This is M by V where for each m=1,...M it maps vertices to integer id for plan
        plan_region_ids_mat.resize(1, std::vector<std::vector<int>>(
                M, std::vector<int>(V, -1)
        ));
        plan_region_order_added_mat.resize(1, std::vector<std::vector<int>>(
                M, std::vector<int>(N, -1)
        ));
    }

    // return List::create(
    //     _["merge_split_steps"] = merge_split_step_vec,
    //     _["count"] = cnnt,
    //     _["region_dvals_mat_list"] = plan_d_vals_mat,
    //     _["region_order_added_list"] = plan_region_order_added_mat,
    //     _["region_ids_mat_list"] = plan_region_ids_mat
    // );

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
    for(int step_num=0; step_num< total_steps; step_num++){
        if(verbosity > 1){
            if(merge_split_step_vec[step_num]){
                Rprintf("Iteration %d: Merge Split Step %d \n", step_num+1, step_num - smc_step_num + 1);
            }else{
                Rprintf("Iteration %d: SMC Step %d \n",  step_num+1, smc_step_num + 1);
            }
        }

        // For the first iteration we need to pass a special previous ancestor thing
        if(step_num == 0){
        std::vector<int> dummy_prev_ancestors(M, 1);
        // split the map
        generalized_split_maps(
            g, counties, cg, pop,
            plans_vec, new_plans_vec,
            original_ancestor_mat.at(step_num),
            parent_index_mat.at(step_num),
            dummy_prev_ancestors,
            unnormalized_sampling_weights,
            normalized_weights_mat.at(step_num),
            draw_tries_mat.at(step_num),
            parent_unsuccessful_tries_mat.at(step_num),
            acceptance_rates.at(step_num),
            nunique_parents_vec.at(step_num),
            nunique_original_ancestors_vec.at(step_num),
            ancestors, lags,
            lower, upper, target,
            k_params.at(smc_step_num), split_district_only,
            pool,
            verbosity
        );

            // For the first ancestor one make every ancestor themselves
            std::iota (parent_index_mat[0].begin(), parent_index_mat[0].end(), 0);
            std::iota (original_ancestor_mat[0].begin(), original_ancestor_mat[0].end(), 0);
            smc_step_num++;
        }else if(merge_split_step_vec[step_num]){ // check if its a merge split step
             // Rprintf("%d!\n", step_num);
            // run merge split 
            // Set the number of steps to run at 1 over previous stage acceptance rate
            // int nsteps_to_run = std::ceil(1/std::pow(acceptance_rates.at(step_num-1),1.5)); // * std::max(merge_split_step_num,1);
            int nsteps_to_run = std::ceil(1/acceptance_rates.at(step_num-1)); // * std::max(merge_split_step_num,1);
            num_merge_split_attempts_vec.at(merge_split_step_num) = nsteps_to_run;

            run_merge_split_step_on_all_plans( 
                pool,
                g, counties, cg, pop,
                plans_vec,
                split_district_only, k_params.at(smc_step_num),
                nsteps_to_run,
                lower, upper, target,
                merge_split_successes_mat.at(merge_split_step_num)
            );

            // set the acceptance rate 
            int total_ms_successes = std::accumulate(
                merge_split_successes_mat.at(merge_split_step_num).begin(), 
                merge_split_successes_mat.at(merge_split_step_num).end(), 
                0);

            int total_ms_attempts = M * nsteps_to_run;

            acceptance_rates.at(step_num) = total_ms_successes / static_cast<double>(total_ms_attempts);

            if (verbosity >= 3){
                Rprintf("Ran %d Merge Split Attempts: %d Successes out of %d attempts. Acceptance Rate: %.2f\n", 
                    nsteps_to_run,  
                    total_ms_successes,total_ms_attempts, 100.0*acceptance_rates.at(step_num));
            }



            // Copy results from previous step
            original_ancestor_mat.at(step_num) = original_ancestor_mat.at(step_num-1);
            parent_index_mat.at(step_num) = parent_index_mat.at(step_num-1);
            draw_tries_mat.at(step_num) = std::vector<int>(M, nsteps_to_run);
            parent_unsuccessful_tries_mat.at(step_num) = parent_unsuccessful_tries_mat.at(step_num-1);
            nunique_parents_vec.at(step_num) = nunique_parents_vec.at(step_num-1);
            nunique_original_ancestors_vec.at(step_num) = nunique_original_ancestors_vec.at(step_num-1);
            merge_split_step_num++;
        }else{ // else just run a normal smc step 
        // split the map and we can use the previous original ancestor matrix row
            generalized_split_maps(
                g, counties, cg, pop,
                plans_vec, new_plans_vec,
                original_ancestor_mat.at(step_num),
                parent_index_mat.at(step_num),
                original_ancestor_mat.at(step_num-1),
                unnormalized_sampling_weights,
                normalized_weights_mat.at(step_num),
                draw_tries_mat.at(step_num),
                parent_unsuccessful_tries_mat.at(step_num),
                acceptance_rates.at(step_num),
                nunique_parents_vec.at(step_num),
                nunique_original_ancestors_vec.at(step_num),
                ancestors, lags,
                lower, upper, target,
                k_params.at(smc_step_num), split_district_only,
                pool,
                verbosity
            );
            smc_step_num++;
        }

        if (verbosity == 1 && CLI_SHOULD_TICK){
            cli_progress_set(bar, step_num);
        }
        Rcpp::checkUserInterrupt();


        // compute log incremental weights and sampling weights for next round
        if(!merge_split_step_vec[step_num]){
        get_all_plans_log_gsmc_weights(
            pool,
            g,
            plans_vec,
            split_district_only,
            log_incremental_weights_mat.at(step_num),
            unnormalized_sampling_weights,
            target,
            pop_temper
        );
        }else{
            log_incremental_weights_mat.at(step_num) = log_incremental_weights_mat.at(step_num-1);
        }



        // compute effective sample size
        n_eff.at(step_num) = compute_n_eff(log_incremental_weights_mat.at(step_num));

        // Now update the diagnostic info if needed, region labels, dval column of the matrix
        if(diagnostic_mode && step_num < total_steps -1 && !split_district_only){ // record if in diagnostic mode and generalized splits
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(step_num).at(j) = plans_vec.at(j).region_ids;
                plan_region_order_added_mat.at(step_num).at(j) = plans_vec.at(j).region_added_order;
                plan_d_vals_mat.at(step_num).at(j) = plans_vec.at(j).region_dvals;
            }
        }else if(diagnostic_mode){ // record if in diagnostic mode but not generalized splits 
            // or if its the last round and so dval doesn't matter
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(step_num).at(j) = plans_vec.at(j).region_ids;
                plan_region_order_added_mat.at(step_num).at(j) = plans_vec.at(j).region_added_order;
            }
        }else if(step_num == total_steps-1){ // else if not diagnostic only record final step
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(0).at(j) = plans_vec.at(j).region_ids;
                plan_region_order_added_mat.at(0).at(j) = plans_vec.at(j).region_added_order;
            }
        }

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
        _["region_ids_mat_list"] = plan_region_ids_mat,
        _["region_dvals_mat_list"] = plan_d_vals_mat,
        _["region_order_added_list"] = plan_region_order_added_mat,
        _["log_incremental_weights_mat"] = log_incremental_weights_mat,
        _["normalized_weights_mat"] = normalized_weights_mat,
        _["draw_tries_mat"] = draw_tries_mat,
        _["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat,
        _["acceptance_rates"] = acceptance_rates,
        _["nunique_parent_indices"] = nunique_parents_vec,
        _["nunique_original_ancestors"] = nunique_original_ancestors_vec,
        _["ancestors"] = ancestors,
        _["step_n_eff"] = n_eff,
        _["merge_split_steps"] = merge_split_step_vec,
        _["merge_split_attempt_counts"] = num_merge_split_attempts_vec,
        _["merge_split_success_counts"] = merge_split_successes_mat
    );

    return out;

}
