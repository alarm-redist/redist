/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Run SMC with MCMC merge-split steps mixed in
********************************************************/

constexpr bool DEBUG_GSMC_PLANS_VERBOSE = false; // Compile-time constant

#include "smc.h"

/* 
 *  Use SMC Sampler method to split a multidistrict in all of the plans
 * 
 *  Using the procedure outlined in <PAPER HERE> this function attempts to split
 *  a multidistrict in a previous step's plan until `M` successful splits have been made. 
 *  This is based on the `split_maps` function in smc.cpp
 * 
 *  @param map_params Map parameters (adj graph, population, etc.)
 *  @param old_plans_ptr_vec A vector of smart pointers (unique) to plans from the
 *  previous step
 *  @param new_plans_ptr_vec A vector of smart pointers (unique) to plans which
 *  will be filled with plans that had a multidistrict split to make them
 *  @param tree_splitters_ptr_vec A vector of smart pointers (unique) to TreeSplitter
 *  objects which will be used to split a plan from `old_plans_ptr_vec` to make a 
 *  plan in `new_plans_ptr_vec`. 
 *  @param parent_index_vec A vector used to track the index of the previous plan
 *  sampled that was successfully split. The value of `parent_index_vec[i]` is the
 *  index of the old plan from which the new plan `new_plans_ptr_vec[i]` was
 *  successfully split from. In other words `new_plans_ptr_vec[i]` is equal to
 *  `attempt_region_split(old_plans_ptr_vec[parent_index_vec[i]], ...)`
 *  @param normalized_cumulative_weights A vector of weights used to sample indices
 *  of the `old_plans_ptr_vec`. The value of `normalized_cumulative_weights[i]` is
 *  the normalized cumulative probability of the weights up to index i 
 *  i.e. the probability that index `i` is selected is 
 *  normalized_cumulative_weights[i]-normalized_cumulative_weights[i-1]
 *  @param draw_tries_vec A vector used to keep track of how many plan split
 *  attempts were made for index i. The value `draw_tries_vec[i]` represents how
 *  many split attempts were made for the i-th new plan (including the successful
 *  split). For example, `draw_tries_vec[i] = 1` means that the first split
 *  attempt was successful.
 *  @param parent_unsuccessful_tries_vec A vector used to keep track of how many times the
 *  previous rounds plans were sampled and unsuccessfully split. The value
 *  `parent_unsuccessful_tries_vec[i]` represents how many times `old_plans_ptr_vec[i]` was sampled
 *  and then unsuccessfully split while creating all `M` of the new plans.
 *  THIS MAY NOT BE THREAD SAFE
 *  @param accept_rate The number of accepted splits over the total number of
 *  attempted splits. This is equal to `sum(draw_tries_vec)/M`
 *  @param n_unique_parent_indices The number of unique parent indices, ie the
 *  number of previous plans that had at least one descendant amongst the new
 *  plans. This is equal to `unique(parent_index_vec)`
 *  @param ancestors Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
 *  WHAT IT IS DOING
 *  @param lags Parameter from older `smc.cpp` code. I DON'T UNDERSTAND
 *  WHAT IT IS DOING
 *  @param split_district_only Whether or not to only allow for single district
 *  splits. If set to `true` will only attempt to split off one district at a
 *  time
 *  @param pool A threadpool for multithreading
 *  @param verbosity A parameter controlling the amount of detail printed out
 *  during the algorithms running
 *  @param diagnostic_level What type of diagnostics to track. Not implemented
 *  yet.
 * 
 *  @details Modifications
 *     - The `new_plans_ptr_vec` is updated with all the newly split plans
 *     - The `old_plans_ptr_vec` is updated with all the newly split plans as well.
 *     Note that the reason both this and `new_plans_ptr_vec` are updated is because
 *     of the nature of the code you need both vectors and so both are passed by
 *     reference to save memory.
 *     - The `original_ancestor_vec` is updated to contain the indices of the
 *     original ancestors of the new plans
 *     - The `parent_index_vec` is updated to contain the indices of the parents of the
 *     new plans
 *     - If two new valid regions are split then the new_region_ids is updated so the
 *     first entry is the first new region and the second entry is the second new region
 *     - The `draw_tries_vec` is updated to contain the number of tries for each
 *     of the new plans
 *     - The `parent_unsuccessful_tries_vec` is updated to contain the number of unsuccessful
 *     samples of the old plans
 *     - The `accept_rate` is updated to contain the average acceptance rate for
 *     this iteration
 *     - `n_unique_parent_indices` and `n_unique_original_ancestors` are updated
 *     with the unique number of parents and original ancestors for all the new
 *     plans respectively
 *     - `ancestors` is updated to something. THIS IS FROM ORIGINAL SMC CODE,
 *     I DO NOT KNOW WHAT IT MEANS
 * 
 */ 
void run_smc_step(
        const MapParams &map_params, const SplittingSchedule &splitting_schedule,
        ScoringFunction const &scoring_function,
        std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
        std::unique_ptr<PlanEnsemble> &old_plan_ensemble,
        std::unique_ptr<PlanEnsemble> &new_plan_ensemble,
        TreeSplitter const &tree_splitters,
        const arma::vec &normalized_cumulative_weights,
        SMCDiagnostics &smc_diagnostics,
        int const smc_step_num, int const step_num,
        umat &ancestors, const std::vector<int> &lags,
        RcppThread::ThreadPool &pool,
        int verbosity, int diagnostic_level
) {
    // important constants
    const int M = old_plan_ensemble->nsims;
    bool const smd_split_district_only = splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnlySMD;

    // PREVIOUS SMC CODE I DONT KNOW WHAT IT DOES
    const int dist_ctr = old_plan_ensemble->plan_ptr_vec.at(0)->num_regions;
    const int n_lags = lags.size();
    umat ancestors_new(M, n_lags); // lags/ancestor thing


    // Because of multithreading we have to add specific checks for if the user
    // wants to quit the program
    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50; // check for interrupts every _ iterations

    // The new region in the split plans is the number of regions in a split plan minus
    // one so the number of regions in a presplit plan
    int new_region_id = old_plan_ensemble->plan_ptr_vec.at(0)->num_regions;

    // we only save for linking edge
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;

    // These are Rcpp::IntegerMatrix::Column type
    Rcpp::IntegerMatrix::Column parent_index_vec = smc_diagnostics.parent_index_mat.column(
        smc_step_num
    );
    Rcpp::IntegerMatrix::Column draw_tries_vec = smc_diagnostics.draw_tries_mat.column(
        step_num
    );
    Rcpp::IntegerMatrix::Column parent_unsuccessful_tries_vec = smc_diagnostics.parent_unsuccessful_tries_mat.column(
        smc_step_num
    );
    
    // count the sizes we draw trees on
    std::vector<std::vector<int>> thread_tree_sizes(rng_states.size(), 
        std::vector<int>(map_params.total_seats, 0)
    );
    // count the sizes of regions successful trees drawn on
    std::vector<std::vector<int>> thread_successful_tree_sizes(rng_states.size(), 
        std::vector<int>(map_params.total_seats, 0)
    );

    // thread safe id counter for seeding RNG generator 
    std::atomic<int> thread_id_counter{0};
    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("About to start SMC Step for %d plans\n", M);
    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        static thread_local int thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
        static thread_local USTSampler ust_sampler(map_params, splitting_schedule);

        // REprintf("Thread ID %d\n", thread_id);
        // REprintf("Plan %d\n\n", i);
        bool ok = false;
        int idx;
        RcppThread::checkUserInterrupt(i % check_int == 0);

        while (!ok) {
            // increase the number of tries for particle i by 1
            draw_tries_vec[i]++;
            // sample previous plan
            idx = rng_states[thread_id].r_int_wgt(normalized_cumulative_weights);
            
            // Get region id the split
            int region_id_to_split;
            if(smd_split_district_only){
                // if just doing district splits just use remainder region
                // which is always the highest id
                region_id_to_split = old_plan_ensemble->plan_ptr_vec[idx]->num_regions - 1;
            }else{ 
                // if generalized split pick a region to try to split
                region_id_to_split = old_plan_ensemble->plan_ptr_vec[idx]->choose_multidistrict_to_split(
                    splitting_schedule.valid_region_sizes_to_split, rng_states[thread_id]
                );
            }
            int region_to_split_size = old_plan_ensemble->plan_ptr_vec[idx]->region_sizes[region_id_to_split];
            
            // int size_of_region_to_split = old_plan_ensemble->plan_ptr_vec[idx]->region_sizes[region_id_to_split];
            // Rprintf("Picked idx %d, Splitting Region size %d. Sizes to try:\n", idx, size_of_region_to_split);
            // for(auto i: splitting_schedule.all_regions_smaller_cut_sizes_to_try[size_of_region_to_split]){
            //     Rprintf("%d, ", i);
            // }
            // Rprintf("\nThe min possible =%d, max possible=%d\n", 
            // splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split].first,
            //     splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split].second);

            //increase the count 
            ++thread_tree_sizes[thread_id][region_to_split_size-1];
            // Rprintf("Count for %d is now %d\n", region_to_split_size, thread_tree_sizes[thread_id][region_to_split_size-1]);


            // Try to split the region 
            std::pair<bool, EdgeCut> edge_search_result = ust_sampler.attempt_to_find_valid_tree_split(
                rng_states[thread_id], tree_splitters,
                *old_plan_ensemble->plan_ptr_vec[idx], region_id_to_split,
                save_edge_selection_prob
            );
        
            // if successful update the new plan
            if(std::get<0>(edge_search_result)){
                if(DEBUG_GSMC_PLANS_VERBOSE){
                    Rprintf("Success, updating Plan %d\n", i);
                } 
                // make the new plan a copy of the old one 
                new_plan_ensemble->plan_ptr_vec[i]->shallow_copy(*old_plan_ensemble->plan_ptr_vec[idx]);
                // now split that region we found on the old one
                new_plan_ensemble->plan_ptr_vec[i]->update_from_successful_split(
                    tree_splitters,
                    ust_sampler, std::get<1>(edge_search_result),
                    region_id_to_split, new_region_id, 
                    true
                );
                // record index of new plan's parent
                parent_index_vec[i] = idx;
                // add as successful tree size 
                ++thread_successful_tree_sizes[thread_id][region_to_split_size-1];
                // make ok true to break out of loop
                ok = true;
            }else{ // else bad sample so try again
                 // check for user interrupt
                 RcppThread::checkUserInterrupt(draw_tries_vec[i] % reject_check_int == 0);
                 // if diagnostic level 2 or higher get unsuccessful count 
                 if(diagnostic_level >= 0){
                     // not atomic so technically not thread safe but doesn't seem to differ in practice
                     parent_unsuccessful_tries_vec[idx]++;
                 }
            }
        }

        // ORIGINAL SMC CODE I DONT KNOW WHAT THIS DOES
        // save ancestors/lags
        for (int j = 0; j < n_lags; j++) {
            if (dist_ctr <= lags[j]) {
                ancestors_new(i, j) = i;
            } else {
                ancestors_new(i, j) = ancestors(idx, j);
            }
        }
        if (verbosity >= 3) {
            ++bar;
        }
        
    });

    // Wait for all the threads to finish
    pool.wait();


    // now swap the old plans with the new ones. This avoids needing to actually copy
    std::swap(old_plan_ensemble, new_plan_ensemble);

    // update tree sizes counts 
    for (size_t region_size = 0; region_size < map_params.total_seats; region_size++)
    {
        for (size_t a_thread_id = 0; a_thread_id < rng_states.size(); a_thread_id++)
        {
            smc_diagnostics.tree_sizes_mat(region_size, step_num) += thread_tree_sizes[a_thread_id][region_size];
            smc_diagnostics.successful_tree_sizes_mat(region_size, step_num) += thread_successful_tree_sizes[a_thread_id][region_size];
        }
    }
    

    // now compute acceptance rate and unique parents and original ancestors
    double accept_rate = M / static_cast<double>(sum(draw_tries_vec));
    smc_diagnostics.acceptance_rates.at(step_num) = accept_rate;

    // Get number of unique parents
    std::set<int> unique_parents(parent_index_vec.begin(), parent_index_vec.end());
    smc_diagnostics.nunique_parents.at(smc_step_num) = unique_parents.size();
    if (verbosity >= 3) {
        Rcout << "  " << std::setprecision(2) << 100.0 * accept_rate << "% acceptance rate, " <<
       100.0 * smc_diagnostics.nunique_parents.at(smc_step_num) / M << "% of previous step's plans survived.\n";
    }

    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;
}


void run_merge_split_step_on_all_plans( 
    RcppThread::ThreadPool &pool,
    MapParams const &map_params, const SplittingSchedule &splitting_schedule,
    ScoringFunction const &scoring_function,
    std::vector<RNGState> &rng_states, SamplingSpace const sampling_space,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plan_ptrs_vec, 
    TreeSplitter const &tree_splitter,
    std::string const merge_prob_type, 
    double const rho, bool const is_final, 
    int const nsteps_to_run,
    int const merge_split_step_num, int const step_num,
    SMCDiagnostics &smc_diagnostics,
    int verbosity
){
    int const num_regions = plan_ptrs_vec[0]->num_regions;
    const int check_int = 50; // check for interrupts every _ iterations
    int nsims = (int) plan_ptrs_vec.size();
    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Going to run %d steps!\n", nsteps_to_run);

    // Diagnostics
    Rcpp::IntegerMatrix::Column success_count_vec = smc_diagnostics.merge_split_successes_mat.column(
        merge_split_step_num
    );

    // count the sizes we draw trees on
    std::vector<std::vector<int>> thread_tree_sizes(rng_states.size(), 
        std::vector<int>(map_params.total_seats, 0)
    );
    std::vector<std::vector<int>> thread_successful_tree_sizes(rng_states.size(), 
        std::vector<int>(map_params.total_seats, 0)
    );


    
    // thread safe id counter for seeding RNG generator 
    std::atomic<int> thread_id_counter{0};
    // create a progress bar
    RcppThread::ProgressBar bar(nsims, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, nsims, [&] (int i) {
        static thread_local int thread_id = thread_id_counter.fetch_add(1, std::memory_order_relaxed);
        static thread_local USTSampler ust_sampler(map_params, splitting_schedule);
        // Create variables needed for each 
        PlanMultigraph current_plan_multigraph(map_params);
        PlanMultigraph proposed_plan_multigraph(map_params);

        // store the number of succesful runs
        success_count_vec[i] = run_merge_split_steps(
            map_params, splitting_schedule, scoring_function,
            rng_states[thread_id], sampling_space,
            *plan_ptrs_vec[i], *new_plan_ptrs_vec[i], 
            ust_sampler, tree_splitter,
            current_plan_multigraph, 
            proposed_plan_multigraph,
            merge_prob_type,
            rho, is_final,
            nsteps_to_run,
            thread_tree_sizes[thread_id], thread_successful_tree_sizes[thread_id]
        );

        if (verbosity >= 3) {
            ++bar;
        }

        RcppThread::checkUserInterrupt(i % check_int == 0);

    });

    // Wait for all the threads to finish
    pool.wait();

    // update tree sizes counts 
    for (size_t region_size = 0; region_size < map_params.total_seats; region_size++)
    {
        for (size_t a_thread_id = 0; a_thread_id < rng_states.size(); a_thread_id++)
        {
            smc_diagnostics.tree_sizes_mat(region_size, step_num) += thread_tree_sizes[a_thread_id][region_size];
            smc_diagnostics.successful_tree_sizes_mat(region_size, step_num) += thread_successful_tree_sizes[a_thread_id][region_size];
        }
    }
    

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
//' @param ndists The number of districts the final plans will have
//' @param adj_list A 0-indexed adjacency list representing the undirected graph
//' which represents the underlying map the plans are to be drawn on
//' @param counties Vector of county labels of each vertex in `g`
//' @param pop A vector of the population associated with each vertex in `g`
//' @param target Ideal population of a valid district. This is what deviance is calculated
//' relative to
//' @param lower Acceptable lower bounds on a valid district's population
//' @param upper Acceptable upper bounds on a valid district's population
//' @param control Named list of additional parameters.
//' @param num_threads The number of threads the threadpool should use
//' @param verbosity What level of detail to print out while the algorithm is
//' running <ADD OPTIONS>
//' @export
List run_redist_gsmc(
    int const nsims, 
    int const total_seats, int const ndists, Rcpp::IntegerVector const district_seat_sizes,
    int const initial_num_regions, 
    List const &adj_list,
    arma::uvec const &counties, const arma::uvec &pop,
    Rcpp::CharacterVector const &step_types,
    double const target, double const lower, double const upper,
    double const rho, // compactness 
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    List const &control, // control has pop temper, and k parameter value, and splitting method are allowed
    List const &constraints, // constraints 
    int const verbosity, int const diagnostic_level,
    Rcpp::IntegerMatrix const &region_id_mat, 
    Rcpp::IntegerMatrix const &region_sizes_mat
){
    if (DEBUG_GSMC_PLANS_VERBOSE) REprintf("Inside c++ code!\n");
    bool diagnostic_mode = diagnostic_level == 1;
    // set the number of threads
    int num_threads = (int) control["num_threads"];
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    // get seq_alpha 
    double weights_alpha = as<double>(control["seq_alpha"]);
    bool const apply_weights_alpha = weights_alpha != 1;

    // re-seed MT so that `set.seed()` works in R
    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    int num_rng_states = num_threads > 0 ? num_threads : 1;
    std::vector<RNGState> rng_states;rng_states.reserve(num_rng_states);
    for (size_t i = 1; i <= num_rng_states; i++)
    {
        // same seed with i*3 long_jumps for state
        rng_states.emplace_back(global_rng_seed, i*3);
    }

    // Set the sampling space 
    SamplingSpace sampling_space = get_sampling_space(sampling_space_str);
    
    // Legacy, in future remove
    RNGState rng_state((int) Rcpp::sample(INT_MAX, 1)[0]);
    global_seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);
    if (DEBUG_GSMC_PLANS_VERBOSE) REprintf("RNG States created!\n");


    // unpack control params
    // lags thing (copied from original smc code, don't understand what its doing)
    std::vector<int> lags = as<std::vector<int>>(control["lags"]); arma::umat ancestors(nsims, lags.size(), fill::zeros);
    // weight type
    std::string wgt_type = as<std::string>(control["weight_type"]);
    // population tempering parameter 
    double pop_temper = as<double>(control["pop_temper"]);


    // total number of steps to run 
    int total_steps = static_cast<int>(step_types.size());
    int total_ms_steps = 0; int total_smc_steps = 0;
    std::vector<bool> merge_split_step_vec(step_types.size());
    for (size_t i = 0; i < step_types.size(); i++)
    {
        if(static_cast<std::string>(step_types.at(i)) == "smc"){
            merge_split_step_vec.at(i) = false;
            total_smc_steps++;
        }else if(static_cast<std::string>(step_types.at(i)) == "ms"){
            merge_split_step_vec.at(i) = true;
            total_ms_steps++;
        }else{
            REprintf("Invalid step type: %s\n", 
            static_cast<std::string>(step_types.at(i)).c_str());
            throw Rcpp::exception("Invalid step type passed!");
        }
    }
    // sanity check we're not splitting more than ndists districts
    if(initial_num_regions + total_smc_steps > ndists){
        REprintf("Trying to do %d splits with %d initial regions will "
        "create more than ndists=%d districts!\n", 
        total_smc_steps, initial_num_regions,ndists);
        throw Rcpp::exception(
            "Desired number of splits will produce more than ndist districts!"
            );
    }
    if (DEBUG_GSMC_PLANS_VERBOSE) REprintf("Step types vec created!\n");

    // see if we are splitting plans all the way or just creating partial plans
    bool splitting_all_the_way = ndists == initial_num_regions + total_smc_steps;

    // multipler for number of merge split steps 
    double ms_steps_multiplier;
    std::string merge_prob_type;
    if(total_ms_steps > 0){
        ms_steps_multiplier = as<double>(control["ms_moves_multiplier"]);
        merge_prob_type = as<std::string>(control["merge_prob_type"]);
    }

    double tol = std::max(target - lower, upper - target) / target;


    // get splitting type 
    SplittingMethodType splitting_method = get_splitting_type(
        static_cast<std::string>(control["splitting_method"])
        );

    // get the splitting size regime
    SplittingSizeScheduleType splitting_size_regime = get_splitting_size_regime(
        static_cast<std::string>(control["splitting_size_regime"])
    );
    auto splitting_schedule_ptr = get_splitting_schedule(
        total_smc_steps, ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        splitting_size_regime, control
    );
    if (DEBUG_GSMC_PLANS_VERBOSE) REprintf("Splitting Schedule Obj created!\n");

    // Whether or not to only do district splits only 
    bool const multi_member_districting = (
        splitting_size_regime == SplittingSizeScheduleType::DistrictOnlyMMD || 
        splitting_size_regime == SplittingSizeScheduleType::AnyValidSizeMMD
    );
    bool split_district_only = splitting_size_regime == SplittingSizeScheduleType::DistrictOnlySMD;
    bool use_graph_plan_space = sampling_space == SamplingSpace::GraphSpace;

    // Do some input checking 
    // Make sure first merge split argument isn't true 
    if(merge_split_step_vec.at(0)){
        throw Rcpp::exception("The first entry of merge_split_step_vec cannot be true.");
    };

    
    // Create map level graph and county level multigraph
    MapParams const map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper);
    int V = map_params.g.size();

    // Now create diagnostic information 
    SMCDiagnostics smc_diagnostics(
        sampling_space, splitting_method,
        splitting_size_regime, 
        merge_split_step_vec,
        V, nsims, ndists, total_seats, initial_num_regions,
        total_smc_steps, total_ms_steps,
        diagnostic_level, splitting_all_the_way, split_district_only
    );

    // Add scoring function (constraints)
    ScoringFunction const scoring_function(
        map_params, constraints, 
        pop_temper);

    // Create a threadpool
    RcppThread::ThreadPool pool(num_threads);

    // Now we add everything here to a scope since it won't be needed for the end
    // create the ensemble 
    std::unique_ptr<PlanEnsemble> plan_ensemble_ptr = get_plan_ensemble_ptr(
        map_params,
        initial_num_regions,
        nsims, sampling_space,
        region_id_mat, region_sizes_mat,
        rng_states, pool, verbosity 
    );


    arma::vec log_weights(nsims, arma::fill::zeros);
    {
    // Ensemble of dummy plans for copying 
    std::unique_ptr<PlanEnsemble> dummy_plan_ensemble_ptr = get_plan_ensemble_ptr(
        map_params,
        initial_num_regions,
        nsims, sampling_space,
        region_id_mat, region_sizes_mat,
        rng_states, pool, verbosity 
    );

    // Get the tree splitter
    std::unique_ptr<TreeSplitter> tree_splitter_ptr = get_tree_splitters(
        map_params,splitting_method, control, nsims
    );

    bool use_naive_k_splitter = splitting_method == SplittingMethodType::NaiveTopK;
    // adaptive k estimation threshold
    bool try_to_estimate_cut_k;
    double thresh;
    // k param values to potentially use. If set to 0 or lower then estimate
    std::vector<int> k_params;
    if(use_naive_k_splitter){
        try_to_estimate_cut_k = as<bool>(control["estimate_cut_k"]);
        if(try_to_estimate_cut_k){
            thresh = (double) control["adapt_k_thresh"];
            k_params.resize(total_smc_steps);
        }else{
            k_params = as<std::vector<int>>(control["manual_k_params"]);
        }
    }
    // Define output variables that must always be created

    // Start off all the unnormalized weights at 1
    arma::vec unnormalized_sampling_weights(nsims, arma::fill::ones);
    arma::vec normalized_cumulative_weights(nsims, arma::fill::value(1.0 / nsims));
    normalized_cumulative_weights = arma::cumsum(normalized_cumulative_weights);
    


    // Loading Info
    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << std::fixed << std::setprecision(4);
        if(!split_district_only){
            Rcout << "GENERALIZED SEQUENTIAL MONTE CARLO";
        }else{
            Rcout << "SEQUENTIAL MONTE CARLO";
        }
        if(total_ms_steps > 0){
            Rcout << " WITH MERGE SPLIT";
        }
        Rcout << "\n";
        Rcout << "Using " << sampling_space_to_str(sampling_space);
        Rcout << " Sampling space to sample " << nsims << " " << V << "-unit ";
        Rcout << "maps with " << ndists << " districts and population between "
              << lower << " and " << upper << " using " 
              << pool.getNumThreads() << " threads, "
              << total_ms_steps << " merge split steps, ";
        if(splitting_size_regime == SplittingSizeScheduleType::DistrictOnlySMD || splitting_size_regime == SplittingSizeScheduleType::DistrictOnlyMMD){
            Rcout << "and only performing 1-district splits.";
        }else if(splitting_size_regime == SplittingSizeScheduleType::AnyValidSizeSMD){
            Rcout << "and generalized region splits.";
        }else if(splitting_size_regime == SplittingSizeScheduleType::OneCustomSize){
            Rcout << "and custom size region splits.";
        }
        Rcout << " Using " << splitting_method_to_str(splitting_method) << " with " <<
        (wgt_type == "optimal" ? "Optimal" : "Simple") << " Weights!\n";
        if (map_params.cg.size() > 1){
            Rcout << "Ensuring no more than " << ndists - 1 << " splits of the "
                  << map_params.cg.size() << " administrative units.\n";
        }
        if(scoring_function.total_soft_constraints > 0){
            Rcout << "Applying " << scoring_function.total_soft_constraints << " constraints.\n";
        }
    }

    // keep track of if we need to swap at the end.
    // counts the number of smc steps
    int smc_step_num = 0;
    int merge_split_step_num = 0;

    std::string bar_fmt = "Split [{cli::pb_current}/{cli::pb_total}] {cli::pb_bar} | ETA{cli::pb_eta}";
    RObject bar = cli_progress_bar(total_steps, cli_config(false, bar_fmt.c_str()));
    // Now for each run through split the map
    try {
    for(int step_num=0; step_num < total_steps; step_num++){
        if(verbosity > 1){
            if(merge_split_step_vec[step_num]){
                Rprintf("Iteration %d: Merge Split Step %d \n", step_num+1, merge_split_step_num + 1);
            }else{
                Rprintf("Iteration %d: SMC Step %d of %d \n",  step_num+1, smc_step_num + 1, total_smc_steps);
            }
        }
        // its the final splitting step if step_num + 1 == total_smc steps
        bool const is_final_splitting_step = step_num + 1 == total_smc_steps;
        //  using std::chrono::high_resolution_clock;
        // using std::chrono::duration_cast;
        // using std::chrono::duration;
        // using std::chrono::milliseconds;

        // Check what step type
        if(!merge_split_step_vec[step_num]){

            // set the splitting schedule 
            splitting_schedule_ptr->set_potential_cut_sizes_for_each_valid_size(
                smc_step_num, plan_ensemble_ptr->plan_ptr_vec[0]->num_regions
                );

            // splitting_schedule_ptr->print_current_step_splitting_info();
            
            // Print if needed
            if(verbosity >= 3 && splitting_size_regime == SplittingSizeScheduleType::OneCustomSize){
                splitting_schedule_ptr->print_current_step_splitting_info();
            }

            // check if k is passed in or estimate 
            if(use_naive_k_splitter){
                if(try_to_estimate_cut_k){
                    // est k
                    int est_cut_k;
                    int last_k = smc_step_num == 0 ? std::max(1, V - 5) : k_params.at(smc_step_num-1);
                    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("About to try to estimate cut k!\n");
                    estimate_cut_k(
                        map_params, *splitting_schedule_ptr, rng_state, 
                        est_cut_k, last_k, unnormalized_sampling_weights, thresh,
                        tol, plan_ensemble_ptr->plan_ptr_vec, 
                        split_district_only, verbosity);
                    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Estimated cut k!\n");
                    k_params.at(smc_step_num) = est_cut_k;

                    if (verbosity >= 3) {
                        Rcout << " (using estimated k = " << k_params.at(smc_step_num) << ")\n";
                    }
                }else{
                    if (verbosity >= 3) {
                        Rcout << " (using input k = " << k_params.at(smc_step_num) << ")\n";
                    }
                }
                tree_splitter_ptr->update_single_int_param(k_params.at(smc_step_num));
                smc_diagnostics.cut_k_values.at(step_num) = tree_splitter_ptr->get_single_int_param();
            }

            
            // auto t1f = high_resolution_clock::now();

            if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("About to run smc step %d!\n", smc_step_num);
            // split the map
            run_smc_step(map_params, *splitting_schedule_ptr, scoring_function,
                rng_states, sampling_space,
                plan_ensemble_ptr, dummy_plan_ensemble_ptr, 
                *tree_splitter_ptr,
                normalized_cumulative_weights,
                smc_diagnostics,
                smc_step_num, step_num,
                ancestors, lags,
                pool,
                verbosity, diagnostic_mode ? 3 : 0
            );
            if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Ran smc step %d!\n", smc_step_num);
            if(smc_step_num == 0 && initial_num_regions == 1){
                // For the first ancestor one make every ancestor themselves
                std::iota(
                    smc_diagnostics.parent_index_mat.column(0).begin(), 
                    smc_diagnostics.parent_index_mat.column(0).end(), 
                    0);
            }
            // plan_ensemble_ptr->plan_ptr_vec[0]->Rprint(true);

            
            // auto t2f = std::chrono::high_resolution_clock::now();
            // /* Getting number of milliseconds as a double. */
            // std::chrono::duration<double, std::milli> ms_doublef = t2f - t1f;
            // Rcout << "Running SMC " << ms_doublef.count() << " ms\n";
            auto t1 = std::chrono::high_resolution_clock::now();

            // compute splitting probability if MMD or if Anysplits SMD and num regions isn't number of districts
            bool compute_log_splitting_prob = (
                multi_member_districting ||
                (
                splitting_schedule_ptr->schedule_type == SplittingSizeScheduleType::AnyValidSizeSMD &&
                plan_ensemble_ptr->plan_ptr_vec[0]->num_regions != ndists
                )
            );
            
            if(wgt_type == "optimal"){
                // TODO make more princicpal in the future 
                // for now its just if not district only and not final round 
                if (verbosity >= 3) Rprintf("Computing Optimal Weights:\n");
                compute_all_plans_log_optimal_incremental_weights(
                    pool,
                    map_params, *splitting_schedule_ptr, sampling_space,
                    scoring_function, rho,
                    plan_ensemble_ptr->plan_ptr_vec, *tree_splitter_ptr,
                    compute_log_splitting_prob, is_final_splitting_step,
                    smc_diagnostics.log_incremental_weights_mat.col(smc_step_num),
                    verbosity
                );
            }else if(wgt_type == "simple"){
                if (verbosity >= 3) Rprintf("Computing Simple Backwards Kernel Weights:\n");
                compute_all_plans_log_simple_incremental_weights(
                    pool,
                    map_params, *splitting_schedule_ptr,
                    sampling_space,
                    scoring_function, rho,
                    plan_ensemble_ptr->plan_ptr_vec, *tree_splitter_ptr,
                    compute_log_splitting_prob, is_final_splitting_step,
                    smc_diagnostics.log_incremental_weights_mat.col(smc_step_num),
                    verbosity
                );
            }else{
                throw Rcpp::exception("invalid weight type!");
            }
            if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Done computing weights!\n");
            
            auto t2 = std::chrono::high_resolution_clock::now();
            /* Getting number of milliseconds as a double. */
            std::chrono::duration<double, std::milli> ms_double = t2 - t1; 
            if(DEBUG_VERBOSE){
                Rcout << "Calculating log weights took" << ms_double.count() << " ms, ";
            }

            // do seq_alpha if needed
            if(apply_weights_alpha){
                // reindex first by making the log weight at index i the index of the parent
                // we use unnormalized_sampling_weights to hold new values then swap later
                for (size_t i = 0; i < nsims; i++)
                {
                    // REprintf("Sending %u to %d \n", i, smc_diagnostics.parent_index_mat(i, smc_step_num));
                    unnormalized_sampling_weights[i] = log_weights[smc_diagnostics.parent_index_mat(i, smc_step_num)];
                }
                std::swap(unnormalized_sampling_weights, log_weights);
                // add incremental weights to the current log weights 
                log_weights = log_weights + smc_diagnostics.log_incremental_weights_mat.col(smc_step_num);

                // if using seq_alpha then our sampling weights for next round are 
                // proportional to exp(alpha* (prev_log_weights + incremental_weights))
                unnormalized_sampling_weights = arma::exp(weights_alpha * log_weights);
                if(!is_final_splitting_step){
                    // if not the end then multiply by 1-alpha
                    log_weights = (1-weights_alpha) * log_weights;
                }
            }else{
                // if no seq alpha then log weights are just the incremental weights
                // and sampling weights are just exp of exponential weights 
                log_weights = smc_diagnostics.log_incremental_weights_mat.col(smc_step_num);
                unnormalized_sampling_weights = arma::exp(log_weights);
            }
            normalized_cumulative_weights =  arma::cumsum(unnormalized_sampling_weights);

            // compute log weight sd
            smc_diagnostics.log_wgt_stddevs.at(smc_step_num) = arma::stddev(log_weights);
            // compute effective sample size
            smc_diagnostics.n_eff.at(smc_step_num) = normalized_cumulative_weights[nsims-1] * normalized_cumulative_weights[nsims-1]  / arma::sum(arma::square(unnormalized_sampling_weights));
            // Now normalize the weights 
            normalized_cumulative_weights = normalized_cumulative_weights / normalized_cumulative_weights(
                nsims-1
            );

            if (verbosity >= 3) {
                Rcout << "  " << std::setprecision(2) 
                      << 100*smc_diagnostics.n_eff.at(smc_step_num)/nsims <<  "% efficiency."
                      << std::setprecision(4)
                      << " Log Weight Standard Deviation: " 
                      << smc_diagnostics.log_wgt_stddevs.at(smc_step_num) << std::endl;
            }

            // only increase if we have smc steps left else it will cause index issues
            // with merge split
            if(smc_step_num < total_smc_steps-1){
                smc_step_num++;
            }
        }else if(merge_split_step_vec[step_num]){ // check if its a merge split step
            // run merge split 
            // Set the number of steps to run at 1 over previous stage acceptance rate if not 0
            int prev_acceptance_index = merge_split_step_num == 0 ? step_num-1 : step_num - 2;
            double prev_acceptance_rate = smc_diagnostics.acceptance_rates.at(prev_acceptance_index);
            // if the acceptance is zero just default to 5
            prev_acceptance_rate = prev_acceptance_rate > 0 ? prev_acceptance_rate : .1;

            int nsteps_to_run = std::ceil(ms_steps_multiplier * std::ceil((1/prev_acceptance_rate))); 
            smc_diagnostics.num_merge_split_attempts_vec.at(merge_split_step_num) = nsteps_to_run;

            if (verbosity >= 3){
                Rprintf("  Running %d Merge Split Steps per plan, %d in total!\n", 
                    nsteps_to_run, nsteps_to_run*nsims);
            }


            splitting_schedule_ptr->update_cut_sizes_for_mergesplit_step(
                smc_step_num, plan_ensemble_ptr->plan_ptr_vec[0]->num_regions
            );


            // auto t1fm = high_resolution_clock::now();
            run_merge_split_step_on_all_plans(
                pool,
                map_params, *splitting_schedule_ptr, 
                scoring_function,
                rng_states, sampling_space, 
                plan_ensemble_ptr->plan_ptr_vec, dummy_plan_ensemble_ptr->plan_ptr_vec,
                *tree_splitter_ptr,
                merge_prob_type, 
                rho, is_final_splitting_step,
                nsteps_to_run,
                merge_split_step_num, step_num,
                smc_diagnostics,
                verbosity
            );

            if(use_naive_k_splitter){
                smc_diagnostics.cut_k_values.at(step_num) = tree_splitter_ptr->get_single_int_param();
            }
            

            // auto t2fm = high_resolution_clock::now();
            // /* Getting number of milliseconds as a double. */
            // duration<double, std::milli> ms_doublefm = t2fm - t1fm;
            // Rcout << "Running Merge split " << ms_doublefm.count() << " ms\n";

            // set the acceptance rate 
            int total_ms_successes = Rcpp::sum(smc_diagnostics.merge_split_successes_mat.column(merge_split_step_num));
            int total_ms_attempts = nsims * nsteps_to_run;

            smc_diagnostics.acceptance_rates.at(step_num) = total_ms_successes / static_cast<double>(total_ms_attempts);

            if (verbosity >= 3){
                Rprintf("  Acceptance Rate: %.2f\n",   
                    100.0*smc_diagnostics.acceptance_rates.at(step_num));
            }

            // Access the column
            IntegerMatrix::Column col = smc_diagnostics.draw_tries_mat(_, step_num);
            // Set all elements in the column to the value nsteps_to_run
            std::fill(col.begin(), col.end(), nsteps_to_run);

            merge_split_step_num++;
        }

        // Add details diagnostics if needed
        // Now update the diagnostic info if needed, region labels, dval column of the matrix
        if(diagnostic_mode){
            smc_diagnostics.add_full_step_diagnostics(
                total_steps, splitting_all_the_way,
                step_num, merge_split_step_num, smc_step_num,
                !merge_split_step_vec[step_num],
                sampling_space,
                pool,
                *plan_ensemble_ptr, *dummy_plan_ensemble_ptr,
                *splitting_schedule_ptr
            );
        }

        if (verbosity == 1 && CLI_SHOULD_TICK){
            cli_progress_set(bar, step_num);
        }
        Rcpp::checkUserInterrupt();

    }
    } catch (Rcpp::internal::InterruptedException e) {
        cli_progress_done(bar);
        return R_NilValue;
    }
    cli_progress_done(bar);
    

    // if not diagnostic reorder the plans 
    if(!diagnostic_mode){
        reorder_all_plans(pool, plan_ensemble_ptr->plan_ptr_vec, dummy_plan_ensemble_ptr->plan_ptr_vec);
    }
    // end of scope
    }

    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Exiting main loop and going to do diagnostics!\n");

    Rcpp::IntegerMatrix plan_mat = plan_ensemble_ptr->get_R_plans_matrix(); // integer matrix to store final plans
    // to try to save memory kill the vector 
    plan_ensemble_ptr->flattened_all_plans.clear(); plan_ensemble_ptr->flattened_all_plans.shrink_to_fit();
    
    std::vector<Rcpp::IntegerMatrix> plan_sizes_mat; // hacky way of potentially passing the plan sizes 
    // mat as output if we are not splitting all the way 
    plan_sizes_mat.reserve(1);
    bool plan_sizes_saved;

    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Plans saved!\n");

    // if only sampling partial plans or MMD plans then return the size matrix
    if(
        !splitting_all_the_way || 
        splitting_schedule_ptr->schedule_type == SplittingSizeScheduleType::AnyValidSizeMMD ||
        splitting_schedule_ptr->schedule_type == SplittingSizeScheduleType::DistrictOnlyMMD
    ){
        if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Getting ready to save region sizes!\n");
        plan_sizes_mat.push_back(plan_ensemble_ptr->get_R_sizes_matrix(pool));
        plan_sizes_saved = true;
    }else{
        // else create a dummy matrix 
        plan_sizes_mat.emplace_back(1,1);
        plan_sizes_saved = false;
    }

    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Plan matrix (and sizes potentially) saved!\n");


    // Return results
    List out = List::create(
        _["plans_mat"] = plan_mat,
        _["region_sizes_mat"] = plan_sizes_mat.at(0),
        _["plan_sizes_saved"] = plan_sizes_saved,
        _["log_weights"] = log_weights,
        _["ancestors"] = ancestors,
        _["step_types"] = step_types,
        _["merge_split_steps"] = merge_split_step_vec
    );

    // add all the diagnostics 
    smc_diagnostics.add_diagnostics_to_out_list(out);

    if(DEBUG_GSMC_PLANS_VERBOSE) Rprintf("Returning to R!\n");
    return out;

}
