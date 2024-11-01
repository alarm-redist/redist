/********************************************************
* Author: Philip O'Sullivan
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Sequential Monte Carlo redistricting sampler with
*           optimal weights
********************************************************/

#include "optimal_gsmc.h"



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
List optimal_gsmc_plans(
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

    double pop_temper = as<double>(control["pop_temper"]);


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
        if(!split_district_only){
            Rcout << "GENERALIZED SEQUENTIAL MONTE CARLO\n";
        }else{
            Rcout << "SEQUENTIAL MONTE CARLO\n";
        }
        Rcout << "Sampling " << M << " " << V << "-unit ";
        Rcout << "maps with " << N << " districts and population between "
              << lower << " and " << upper << " using " << num_threads << " threads ";
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

    // This is N-1 by M where [i][j] is the number of times particle j from the
    // previous round was sampled and unsuccessfully split on iteration i so this
    // does not count successful sample then split
    std::vector<std::vector<int>> parent_unsuccessful_tries_mat(N-1, std::vector<int> (M, 0));

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

        // reserve space for N-1 elements
        plan_d_vals_mat.reserve(N-1);

        // This is N-2 by M by 1,2,...N where for each n=1,...,N-1, m=1,...,M it maps the region
        // id to the region's d value. So for a given n it is n by M by n+1
        // It stops at N-2 because for N-1 its all 1
        for(int n = 1; n < N-1; n++){
            plan_d_vals_mat.push_back(
                std::vector<std::vector<int>>(M, std::vector<int>(n+1, -1))
            );
        }

        //  M by N-1 where for each m=1,...M it maps region ids to region labels in a plan
        final_plan_region_labels.resize(
                M, std::vector<std::string>(N, "MISSING")
        );

    }else{ // else only track for final round
        // Create info tracking we will pass out at the end
        // This is M by V where for each m=1,...M it maps vertices to integer id for plan
        plan_region_ids_mat.resize(1, std::vector<std::vector<int>>(
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
            unnormalized_sampling_weights,
            normalized_weights_mat[n],
            draw_tries_mat[n],
            parent_unsuccessful_tries_mat.at(n),
            acceptance_rates[n],
            nunique_parents_vec[n],
            nunique_original_ancestors_vec[n],
            ancestors, lags,
            lower, upper, target,
            k_params[n], split_district_only,
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
            unnormalized_sampling_weights,
            normalized_weights_mat[n],
            draw_tries_mat[n],
            parent_unsuccessful_tries_mat.at(n),
            acceptance_rates[n],
            nunique_parents_vec[n],
            nunique_original_ancestors_vec[n],
            ancestors, lags,
            lower, upper, target,
            k_params[n], split_district_only,
            pool,
            verbosity
        );
        }

        if (verbosity == 1 && CLI_SHOULD_TICK){
            cli_progress_set(bar, n);
        }
        Rcpp::checkUserInterrupt();


        // compute log incremental weights and sampling weights for next round
        get_all_plans_log_gsmc_weights(
            pool,
            g,
            new_plans_vec,
            split_district_only,
            log_incremental_weights_mat.at(n),
            unnormalized_sampling_weights,
            target,
            pop_temper
        );


        // compute effective sample size
        n_eff.at(n) = compute_n_eff(log_incremental_weights_mat[n]);

        // Now update the diagnostic info if needed, region labels, dval column of the matrix
        if(diagnostic_mode && n == N-2){ // record if in diagnostic mode and final step
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(n).at(j) = plans_vec[j].region_num_ids;
                final_plan_region_labels.at(j) = plans_vec[j].region_str_labels;
            }
        }else if(diagnostic_mode){ // record if in diagnostic mode but not final step
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(n).at(j) = plans_vec[j].region_num_ids;
                plan_d_vals_mat.at(n).at(j) = plans_vec[j].region_dvals;
            }
        }else if(n == N-2){ // else if not only record final step
            for(int j=0; j<M; j++){
                plan_region_ids_mat.at(0).at(j) = plans_vec[j].region_num_ids;
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
        _["final_region_labs"] = final_plan_region_labels,
        _["region_ids_mat_list"] = plan_region_ids_mat,
        _["region_dvals_mat_list"] = plan_d_vals_mat,
        _["log_incremental_weights_mat"] = log_incremental_weights_mat,
        _["normalized_weights_mat"] = normalized_weights_mat,
        _["draw_tries_mat"] = draw_tries_mat,
        _["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat,
        _["acceptance_rates"] = acceptance_rates,
        _["nunique_parent_indices"] = nunique_parents_vec,
        _["nunique_original_ancestors"] = nunique_original_ancestors_vec,
        _["ancestors"] = ancestors,
        _["step_n_eff"] = n_eff
    );

    return out;

}




