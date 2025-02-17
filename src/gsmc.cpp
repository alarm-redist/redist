/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Run SMC with MCMC merge-split steps mixed in
********************************************************/

#include "gsmc.h"

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
        std::vector<std::unique_ptr<Plan>> &old_plans_ptr_vec, 
        std::vector<std::unique_ptr<Plan>> &new_plans_ptr_vec,
        std::vector<std::unique_ptr<TreeSplitter>> &tree_splitters_ptr_vec,
        Rcpp::IntegerMatrix::Column parent_index_vec,
        const arma::vec &normalized_cumulative_weights,
        Rcpp::IntegerMatrix::Column draw_tries_vec,
        Rcpp::IntegerMatrix::Column parent_unsuccessful_tries_vec,
        double &accept_rate,
        int &n_unique_parent_indices,
        umat &ancestors, const std::vector<int> &lags,
        bool const split_district_only,
        RcppThread::ThreadPool &pool,
        int verbosity, int diagnostic_level
                ) {
    // important constants
    const int V = map_params.V;
    const int M = static_cast<int>(old_plans_ptr_vec.size());

    // PREVIOUS SMC CODE I DONT KNOW WHAT IT DOES
    const int dist_ctr = old_plans_ptr_vec.at(0)->num_regions;
    const int n_lags = lags.size();
    umat ancestors_new(M, n_lags); // lags/ancestor thing


    // Because of multithreading we have to add specific checks for if the user
    // wants to quit the program
    const int reject_check_int = 200; // check for interrupts every _ rejections
    const int check_int = 50; // check for interrupts every _ iterations

    // The new region in the split plans is the number of regions in a split plan minus
    // one so the number of regions in a presplit plan
    int new_region_id = old_plans_ptr_vec.at(0)->num_regions;

    // create a progress bar
    RcppThread::ProgressBar bar(M, 1);
    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // REprintf("Plan %d\n\n", i);
        int reject_ct = 0;
        bool ok = false;
        int idx;
        int region_id_to_split, size_of_region_to_split;

        Tree ust = init_tree(V);
        std::vector<bool> visited(V);
        std::vector<bool> ignore(V);
        while (!ok) {
            // increase the number of tries for particle i by 1
            draw_tries_vec[i]++;
            // sample previous plan
            idx = r_int_wgt(M, normalized_cumulative_weights);
            
            *new_plans_ptr_vec.at(i) = *old_plans_ptr_vec.at(idx);

            if(split_district_only){
                // if just doing district splits just use remainder region
                region_id_to_split = new_plans_ptr_vec.at(i)->remainder_region;
            }else{                
                // if generalized split pick a region to try to split
                new_plans_ptr_vec.at(i)->choose_multidistrict_to_split(
                    region_id_to_split, splitting_schedule.valid_region_sizes_to_split
                );
            }
            size_of_region_to_split = new_plans_ptr_vec.at(i)->region_sizes(region_id_to_split);

            // Rprintf("Splitting Region size %d. Sizes to try:\n", size_of_region_to_split);
            // for(auto i: splitting_schedule.all_regions_smaller_cut_sizes_to_try[size_of_region_to_split]){
            //     Rprintf("%d, ", i);
            // }
            // Rprintf("\nThe min possible =%d, max possible=%d\n", 
            // splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split][0],
            //     splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split][1]);
            

            // Now try to split that region
            ok = new_plans_ptr_vec.at(i)->attempt_split(
                map_params, splitting_schedule, ust, *tree_splitters_ptr_vec.at(i), 
                visited, ignore, 
                splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split][0],
                splitting_schedule.all_regions_min_and_max_possible_cut_sizes[size_of_region_to_split][1],
                splitting_schedule.all_regions_smaller_cut_sizes_to_try[size_of_region_to_split],
                split_district_only,
                region_id_to_split, new_region_id
            );

            if(draw_tries_vec[i] == 40000 && draw_tries_vec[i] % 250 == 0){
                    Rprintf("Iter %d\n\n", i);
                    print_tree(ust);
            }

            
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
        if (verbosity >= 3) {
            ++bar;
        }
        RcppThread::checkUserInterrupt(i % check_int == 0);
    });

    // Wait for all the threads to finish
    pool.wait();


    // now replace the old plans with the new ones
    for(int i=0; i < M; i++){
        *old_plans_ptr_vec.at(i) = *new_plans_ptr_vec.at(i);
    }

    // now compute acceptance rate and unique parents and original ancestors
    accept_rate = M / static_cast<double>(sum(draw_tries_vec));

    // Get number of unique parents
    std::set<int> unique_parents(parent_index_vec.begin(), parent_index_vec.end());
    n_unique_parent_indices = unique_parents.size();
    if (verbosity >= 3) {
        Rcout << "  " << std::setprecision(2) << 100.0 * accept_rate << "% acceptance rate, " <<
       100.0 * n_unique_parent_indices / M << "% of previous step's plans survived.\n";
    }


    // ORIGINAL SMC CODE I DONT KNOW WHAT IT DOES
    ancestors = ancestors_new;

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
        int ndists, List adj_list,
        const arma::uvec &counties, const arma::uvec &pop,
        Rcpp::CharacterVector step_types,
        double target, double lower, double upper,
        arma::umat region_id_mat, arma::umat region_sizes_mat,
        std::string sampling_space,
        List control, // control has pop temper, and k parameter value, and whether only district splits are allowed
        int verbosity, bool diagnostic_mode){
    // re-seed MT so that `set.seed()` works in R
    seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    int nsims = static_cast<int>(region_id_mat.n_cols);
    // get the number of implied regions in the plans
    auto a_col = region_id_mat.col(0);
    std::unordered_set<arma::uword> unique_values;
    for (size_t i = 0; i < a_col.n_elem; ++i) {
        unique_values.insert(a_col(i));
    }
    int initial_num_regions = static_cast<int>(unique_values.size());

    
    // unpack control params
    // lags thing (copied from original smc code, don't understand what its doing)
    std::vector<int> lags = as<std::vector<int>>(control["lags"]); arma::umat ancestors(nsims, lags.size(), fill::zeros);


    // weight type
    std::string wgt_type = as<std::string>(control["weight_type"]);
    // population tempering parameter 
    double pop_temper = as<double>(control["pop_temper"]);

    // This is a total_ms_steps by nsims vector where [s][i] is the number of 
    // successful merge splits performed for plan i on merge split round s

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
    
    // see if we are splitting plans all the way or just creating partial plans
    bool splitting_all_the_way = ndists == initial_num_regions + total_smc_steps;

    // multipler for number of merge split steps 
    int ms_steps_multiplier = as<int>(control["ms_steps_multiplier"]);
    std::string merge_prob_type = as<std::string>(control["merge_prob_type"]);
    // set the number of threads
    int num_threads = (int) control["num_threads"];
    if (num_threads <= 0) num_threads = std::thread::hardware_concurrency();
    if (num_threads == 1) num_threads = 0;

    double tol = std::max(target - lower, upper - target) / target;

    // get the k used 
    Rcpp::IntegerVector cut_k_values(total_steps);

    // get splitting type 
    SplittingMethodType splitting_method = get_splitting_type(
        static_cast<std::string>(control["splitting_method"])
        );

    // get the splitting size regime
    SplittingSizeScheduleType splitting_size_regime = get_splitting_size_regime(
        static_cast<std::string>(control["splitting_size_regime"])
    );

    SplittingSchedule splitting_schedule(total_smc_steps, ndists, initial_num_regions, splitting_size_regime, control);

    std::vector<int> min_region_cut_sizes(total_smc_steps,0);
    std::vector<int> max_region_cut_sizes(total_smc_steps,0);
    // Whether or not to only do district splits only 
    bool split_district_only = splitting_size_regime == SplittingSizeScheduleType::DistrictOnly;

    bool use_graph_plan_space = sampling_space == "graph_plan_space";


    // Do some input checking 
    
    // Make sure first merge split argument isn't true 
    if(merge_split_step_vec.at(0)){
        throw Rcpp::exception("The first entry of merge_split_step_vec cannot be true.");
    };
    // TODO Check k params 

    // Create map level graph and county level multigraph
    MapParams map_params(adj_list, counties, pop, ndists, lower, target, upper);
    int V = map_params.g.size();


    // Now create diagnostic information 

    // Level 0
    std::vector<double> log_wgt_stddevs(total_smc_steps); // log weight std devs
    std::vector<double> acceptance_rates(total_steps); // Tracks the acceptance rate - total number of tries over nsims - for each round
    Rcpp::IntegerVector nunique_parents(total_smc_steps); // number of unique parents
    Rcpp::IntegerVector k_vals_used(total_steps); // k value used at each step
    std::vector<std::string> k_val_type(total_steps); // whether k val was estimated or passed in
    std::vector<double> n_eff(total_smc_steps, -1.0); // Tracks the effective sample size for the weights of each round
    Rcpp::IntegerMatrix plan_mat(V, nsims); // integer matrix to store final plans
    std::vector<Rcpp::IntegerMatrix> plan_sizes_mat; // hacky way of potentially passing the plan sizes 
    // mat as output if we are not splitting all the way 
    plan_sizes_mat.reserve(1);

    // Level 1
    // These are all nsims by number of smc steps 
    arma::dmat log_incremental_weights_mat(nsims, total_smc_steps, arma::fill::none); // entry [i][s] is the normalized weight of particle i AFTER split s
    Rcpp::IntegerMatrix draw_tries_mat(nsims, total_steps); // Entry [i][s] is the number of tries it took to form particle i on split s
    draw_tries_mat.fill(0); // fill it with zero
    Rcpp::IntegerMatrix parent_index_mat(nsims, total_smc_steps); // Entry [i][s] is the index of the parent of particle i at split s
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of 
    // successful merge splits performed for plan i on merge split round s
    Rcpp::IntegerMatrix merge_split_successes_mat(nsims, total_ms_steps);
    
    // Level 2
    Rcpp::IntegerMatrix parent_unsuccessful_tries_mat(nsims, total_smc_steps);
    parent_unsuccessful_tries_mat.fill(0);

    // Level 3 
    std::vector<Rcpp::IntegerMatrix> all_steps_plan_region_ids_list;
    all_steps_plan_region_ids_list.reserve(diagnostic_mode ? total_steps : 0);
    std::vector<std::vector<Graph>> all_steps_forests_adj_list;
    all_steps_forests_adj_list.resize((diagnostic_mode && !use_graph_plan_space) ? total_steps : 0);
    std::vector<std::vector<int>> all_steps_valid_region_sizes_to_split;
    all_steps_valid_region_sizes_to_split.resize(diagnostic_mode ? total_smc_steps : 0);
    std::vector<std::vector<int>> all_steps_valid_split_region_sizes;
    all_steps_valid_split_region_sizes.resize(diagnostic_mode ? total_smc_steps : 0);

    
    // Store size at every step but last one if needed
    int plan_dval_list_size = (diagnostic_mode & !split_district_only) ? total_steps-1 : 0;
    if(!splitting_all_the_way) plan_dval_list_size++;
    std::vector<Rcpp::IntegerMatrix> region_sizes_mat_list;
    region_sizes_mat_list.reserve(plan_dval_list_size);

    // If diagnostic mode track vertex region ids from every round
    if(diagnostic_mode){
        // The number of regions starts at 1
        int curr_num_regions = initial_num_regions;
        for (size_t i = 0; i < total_steps; i++)
        {
            all_steps_plan_region_ids_list.emplace_back(V, nsims);
            // Create V by nsims matrix for the plan
            // This is a vector where every entry is a V by nsims Rcpp::IntegerMatrix

            // increase number of regions by 1 if that step is an smc one
            if(!merge_split_step_vec.at(i)) curr_num_regions++;

            // If not doing district only splits, and its not the final one or 
            // we're only doing partial plans then make size matrix
            if(!split_district_only && (i < total_steps-1 || !splitting_all_the_way) ){
                // This is number of regions by nsims
                region_sizes_mat_list.emplace_back(curr_num_regions, nsims);
            }
        }
    }


    // Merge split related diagnostics

    // For each merge split step this counts the number of attempts that were made
    std::vector<int> num_merge_split_attempts_vec(total_ms_steps, -1);


    // Create copies of the matrices
    arma::umat dummy_region_id_mat = region_id_mat;
    arma::umat dummy_region_sizes_mat = region_sizes_mat;


    // Create a vector of pointers to Plan objects
    // b/c we're using abstract classes we must use pointers to the base class
    std::vector<std::unique_ptr<Plan>> plans_ptr_vec; plans_ptr_vec.reserve(nsims);
    std::vector<std::unique_ptr<Plan>> new_plans_ptr_vec; new_plans_ptr_vec.reserve(nsims);

    // Vector of splitters
    std::vector<std::unique_ptr<TreeSplitter>> tree_splitters_ptr_vec = get_tree_splitters(
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

    // Loop over each column of region_id_mat
    for (size_t i = 0; i < region_id_mat.n_cols; ++i) {
        // Create the underlying object for each unique pointer
        if(sampling_space == "graph_plan_space"){
        plans_ptr_vec.emplace_back(
            std::make_unique<GraphPlan>(
                region_id_mat.col(i), region_sizes_mat.col(i), 
                ndists, initial_num_regions,
                map_params.pop, split_district_only
            ));
        new_plans_ptr_vec.emplace_back(
            std::make_unique<GraphPlan>(
                dummy_region_id_mat.col(i), dummy_region_sizes_mat.col(i), 
                ndists, initial_num_regions,
                map_params.pop, split_district_only
            ));
        use_graph_plan_space = true;
        }else if(sampling_space == "spanning_forest_space"){
        plans_ptr_vec.emplace_back(
            std::make_unique<ForestPlan>(
                region_id_mat.col(i), region_sizes_mat.col(i), 
                ndists, initial_num_regions,
                map_params.pop, split_district_only
            ));
        new_plans_ptr_vec.emplace_back(
            std::make_unique<ForestPlan>(
                dummy_region_id_mat.col(i), dummy_region_sizes_mat.col(i), 
                ndists, initial_num_regions,
                map_params.pop, split_district_only
            ));
        }else{
            throw Rcpp::exception("Input is invalid\n");
        }
    }


    // Define output variables that must always be created

    // Start off all the unnormalized weights at 1
    std::vector<double> unnormalized_sampling_weights(nsims, 1.0);
    
    // Create a threadpool
    RcppThread::ThreadPool pool(num_threads);

    // Loading Info
    if (verbosity >= 1) {
        Rcout.imbue(std::locale(""));
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
        Rcout << "Sampling " << nsims << " " << V << "-unit ";
        Rcout << "maps with " << ndists << " districts and population between "
              << lower << " and " << upper << " using " 
              << (num_threads == 0 ? 1 : num_threads) << " threads, "
              << total_ms_steps << " merge split steps, ";
        if(splitting_size_regime == SplittingSizeScheduleType::DistrictOnly){
            Rcout << "and only performing 1-district splits.";
        }else if(splitting_size_regime == SplittingSizeScheduleType::AnyValidSize){
            Rcout << "and generalized region splits.";
        }else if(splitting_size_regime == SplittingSizeScheduleType::CustomSizes){
            Rcout << "and custom size region splits.";
        }
        Rcout << " Using " << splitting_method_to_str(splitting_method) << "!\n";
        if (map_params.cg.size() > 1){
            Rcout << "Ensuring no more than " << ndists - 1 << " splits of the "
                  << map_params.cg.size() << " administrative units.\n";
        }
    }

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
        //  using std::chrono::high_resolution_clock;
        // using std::chrono::duration_cast;
        // using std::chrono::duration;
        // using std::chrono::milliseconds;
        
        int prev_k;

        // Check what step type
        if(!merge_split_step_vec[step_num]){

            // set the splitting schedule 
            splitting_schedule.set_potential_cut_sizes_for_each_valid_size(
                smc_step_num, plans_ptr_vec[0]->num_regions
                );

            if(verbosity >= 3){
            Rprintf("The following region sizes can split:\n");
            for (size_t i = 1; i <= splitting_schedule.ndists; i++)
            {                
                if(splitting_schedule.valid_region_sizes_to_split[i]){
                Rprintf("\tRegion Size %d | ", (int) i);
                Rprintf("min/max (%d, %d)", 
                    splitting_schedule.all_regions_min_and_max_possible_cut_sizes[i][0],
                    splitting_schedule.all_regions_min_and_max_possible_cut_sizes[i][1]);
                Rprintf(" | possible split sizes: ");
                for(auto smaller_size: splitting_schedule.all_regions_smaller_cut_sizes_to_try[i]){
                    Rprintf("(%d, %d) ", smaller_size, i-smaller_size);
                }
                Rprintf("\n");
                }
                
            }
            Rprintf("All Possible Split Sizes are:(");
            for (size_t i = 1; i <= splitting_schedule.ndists; i++)
            {
                if(splitting_schedule.valid_split_region_sizes[i]){
                    Rprintf("%d, ", i);
                }
            }
            Rprintf(")\n");
            }

            // check if k is passed in or estimate 
            if(use_naive_k_splitter){
            if(try_to_estimate_cut_k){
                // est k
                int est_cut_k;
                int last_k = smc_step_num == 0 ? std::max(1, V - 5) : k_params.at(smc_step_num-1);
                // double thresh = (double) control["adapt_k_thresh"];

                estimate_cut_k(
                    map_params, splitting_schedule, 
                    est_cut_k, last_k, unnormalized_sampling_weights, thresh,
                    tol, plans_ptr_vec, 
                    split_district_only, verbosity);

                k_params.at(smc_step_num) = est_cut_k;

                if (verbosity >= 3) {
                    Rcout << " (using estimated k = " << k_params.at(smc_step_num) << ")\n";
                }
            }else{
                if (verbosity >= 3) {
                    Rcout << " (using input k = " << k_params.at(smc_step_num) << ")\n";
                }
            }
            }

            if(use_naive_k_splitter){
            for (size_t j = 0; j < tree_splitters_ptr_vec.size(); j++)
            {
                tree_splitters_ptr_vec.at(j)->update_single_int_param(k_params.at(smc_step_num));
            }
            }
            

            // auto t1f = high_resolution_clock::now();

            arma::vec cumulative_weights = arma::cumsum(
                arma::vec(unnormalized_sampling_weights)
            );

            arma::vec normalized_cumulative_weights = cumulative_weights / cumulative_weights(
                cumulative_weights.size()-1
            );
            
            // split the map
            run_smc_step(map_params, splitting_schedule,
                plans_ptr_vec, new_plans_ptr_vec, 
                tree_splitters_ptr_vec,
                parent_index_mat.column(smc_step_num),
                normalized_cumulative_weights,
                draw_tries_mat.column(step_num),
                parent_unsuccessful_tries_mat.column(smc_step_num),
                acceptance_rates.at(step_num),
                nunique_parents.at(smc_step_num),
                ancestors, lags,
                split_district_only,
                pool,
                verbosity, diagnostic_mode ? 3 : 0
            );


            if(use_naive_k_splitter){
                cut_k_values.at(step_num) = k_params.at(smc_step_num);
                prev_k = k_params.at(smc_step_num);
            }
            


            // for (size_t i = 0; i < plans_ptr_vec.size(); i++)
            // {
            //     plans_ptr_vec.at(i)->Rprint();
            // }
            


            // auto t2f = std::chrono::high_resolution_clock::now();
            // /* Getting number of milliseconds as a double. */
            // std::chrono::duration<double, std::milli> ms_doublef = t2f - t1f;
            // Rcout << "Running SMC " << ms_doublef.count() << " ms\n";
            auto t1 = std::chrono::high_resolution_clock::now();

            if(!use_graph_plan_space){
                if (verbosity >= 3) Rprintf("Computing Weights:\n");
                get_all_forest_plans_log_optimal_weights(
                    pool,
                    map_params, splitting_schedule, 
                    plans_ptr_vec, tree_splitters_ptr_vec,
                    split_district_only,
                    log_incremental_weights_mat.col(smc_step_num),
                    unnormalized_sampling_weights,
                    pop_temper, verbosity
            );
            }else if(wgt_type == "optimal"){
                // compute log incremental weights and sampling weights for next round
                get_all_plans_log_optimal_weights(
                    pool,
                    map_params, splitting_schedule,  
                    plans_ptr_vec,
                    split_district_only,
                    log_incremental_weights_mat.col(smc_step_num),
                    unnormalized_sampling_weights,
                    pop_temper
                );
            }else if(wgt_type == "adj_uniform"){
                get_all_plans_uniform_adj_weights(
                    pool,
                    splitting_schedule,
                    map_params.g,
                    plans_ptr_vec,
                    split_district_only,
                    log_incremental_weights_mat.col(smc_step_num),
                    unnormalized_sampling_weights,
                    target,
                    pop_temper
                );
            }else{
                throw Rcpp::exception("invalid weight type!");
            }
            
            
            auto t2 = std::chrono::high_resolution_clock::now();
            /* Getting number of milliseconds as a double. */
            std::chrono::duration<double, std::milli> ms_double = t2 - t1; 
            if(DEBUG_VERBOSE){
                Rcout << "Calculating log weights took" << ms_double.count() << " ms, ";
            }


            // compute log weight sd
            log_wgt_stddevs.at(smc_step_num) = arma::stddev(log_incremental_weights_mat.col(smc_step_num));
            // compute effective sample size
            n_eff.at(smc_step_num) = compute_n_eff(log_incremental_weights_mat.col(smc_step_num));

            if (verbosity >= 3) {
                Rcout << std::setprecision(1) << 100*n_eff.at(smc_step_num)/nsims <<  "% efficiency.\n";
            }

            if(smc_step_num == 0 && initial_num_regions == 1){
                // For the first ancestor one make every ancestor themselves
                std::iota(
                    parent_index_mat.column(0).begin(), 
                    parent_index_mat.column(0).end(), 
                    0);
            }

            if(diagnostic_mode){
                int current_num_regions = plans_ptr_vec[0]->num_regions;
                // save the acceptable split sizes 
                for (int region_size = 1; region_size <= ndists - current_num_regions + 2; region_size++)
                {
                    if(splitting_schedule.valid_split_region_sizes[region_size]){
                        all_steps_valid_split_region_sizes[smc_step_num].push_back(region_size);
                    }
                    if(splitting_schedule.valid_region_sizes_to_split[region_size]){
                        
                        all_steps_valid_region_sizes_to_split[smc_step_num].push_back(region_size);;
                    }
                }
                
            }

            // only increase if we have smc steps left else it will cause index issues
            // with merge split
            if(smc_step_num < total_smc_steps-1){
                smc_step_num++;
            }
        }else if(merge_split_step_vec[step_num]){ // check if its a merge split step
            // run merge split 
            // Set the number of steps to run at 1 over previous stage acceptance rate

            // TODO formalize this more 
            int prev_acceptance_index = merge_split_step_num == 0 ? step_num-1 : step_num - 2;

            // int nsteps_to_run = std::ceil(1/std::pow(acceptance_rates.at(step_num-1),1.5)); // * std::max(merge_split_step_num,1);
            int nsteps_to_run = ms_steps_multiplier * std::ceil(1/acceptance_rates.at(prev_acceptance_index)); // * std::max(merge_split_step_num,1);
            num_merge_split_attempts_vec.at(merge_split_step_num) = nsteps_to_run;

            if (verbosity >= 3){
                Rprintf("Ran %d Merge Split Attempts: ", 
                    nsteps_to_run);
            }

            // auto t1fm = high_resolution_clock::now();
            run_merge_split_step_on_all_plans(
                pool,
                map_params, splitting_schedule,
                plans_ptr_vec, new_plans_ptr_vec,
                tree_splitters_ptr_vec,
                split_district_only, merge_prob_type,
                nsteps_to_run,
                merge_split_successes_mat.column(merge_split_step_num)
            );


            cut_k_values.at(step_num) = prev_k;

            // auto t2fm = high_resolution_clock::now();
            // /* Getting number of milliseconds as a double. */
            // duration<double, std::milli> ms_doublefm = t2fm - t1fm;
            // Rcout << "Running Merge split " << ms_doublefm.count() << " ms\n";

            // set the acceptance rate 
            int total_ms_successes = Rcpp::sum(merge_split_successes_mat.column(merge_split_step_num));
            int total_ms_attempts = nsims * nsteps_to_run;

            acceptance_rates.at(step_num) = total_ms_successes / static_cast<double>(total_ms_attempts);

            if (verbosity >= 3){
                Rprintf("%d Successes out of %d attempts. Acceptance Rate: %.2f\n",   
                    total_ms_successes,total_ms_attempts, 100.0*acceptance_rates.at(step_num));
            }

            // Access the column
            IntegerMatrix::Column col = draw_tries_mat(_, step_num);
            
            // Set all elements in the column to the value nsteps_to_run
            std::fill(col.begin(), col.end(), nsteps_to_run);

            merge_split_step_num++;
        }


        
        // Now update the diagnostic info if needed, region labels, dval column of the matrix
        if(diagnostic_mode){ // record if in diagnostic mode and generalized splits
            // reorder the plans by oldest split if either we've done any merge split or
            // its generalized region splits 

            if(merge_split_step_num > 0 || !split_district_only){
                reorder_all_plans(pool, plans_ptr_vec, new_plans_ptr_vec);
            }

            // Copy the vertex plan matrix 
            copy_arma_to_rcpp_mat(pool, 
                region_id_mat.submat(0,0, region_id_mat.n_rows-1, region_id_mat.n_cols-1), 
                all_steps_plan_region_ids_list.at(step_num));
            

            // store the 
            if(!use_graph_plan_space){
                all_steps_forests_adj_list.at(step_num).reserve(nsims);
                for (size_t i = 0; i < nsims; i++)
                {
                    // add the forests from each plan at this step
                    all_steps_forests_adj_list.at(step_num).push_back(
                        plans_ptr_vec.at(i)->get_forest_adj()
                    );
                }
            }

            
            // Copy the sizes if neccesary 
            if(!split_district_only && (step_num < total_steps-1 || !splitting_all_the_way)){
                copy_arma_to_rcpp_mat(
                    pool, 
                    region_sizes_mat.submat(
                    0, 0,
                    plans_ptr_vec.at(0)->num_regions-1,
                    region_sizes_mat.n_cols-1
                    ), 
                    region_sizes_mat_list.at(step_num));
            }
            

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
        reorder_all_plans(pool, plans_ptr_vec, new_plans_ptr_vec);
    }

    // copy the plans into the plan_mat for returning
    copy_arma_to_rcpp_mat(pool, 
        region_id_mat.submat(0,0, region_id_mat.n_rows-1, region_id_mat.n_cols-1), 
        plan_mat);

    // if only sampling partial plans then return the size matrix
    if(!splitting_all_the_way){
        int num_final_regions = plans_ptr_vec.at(0)->num_regions;
        plan_sizes_mat.emplace_back(num_final_regions, nsims);
        copy_arma_to_rcpp_mat(
            pool, 
            region_sizes_mat.submat(
            0, 0,
            num_final_regions-1,
            region_sizes_mat.n_cols-1
            ), 
            plan_sizes_mat.at(0));
        // Rcpp::IntegerMatrix plan_mat(V, nsims); // integer matrix to store final plans
    }else{
        // else create a dummy matrix 
        plan_sizes_mat.emplace_back(1,1);
    }


    // Return results
    List out = List::create(
        _["plans_mat"] = plan_mat,
        _["plan_sizes_mat"] = plan_sizes_mat.at(0),
        _["parent_index"] = parent_index_mat,
        _["region_ids_mat_list"] = all_steps_plan_region_ids_list,
        _["region_sizes_mat_list"] = region_sizes_mat_list, 
        _["forest_adjs_list"] = all_steps_forests_adj_list,
        _["log_incremental_weights_mat"] = log_incremental_weights_mat,
        _["draw_tries_mat"] = draw_tries_mat,
        _["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat,
        _["acceptance_rates"] = acceptance_rates,
        _["nunique_parent_indices"] = nunique_parents,
        _["ancestors"] = ancestors,
        _["step_n_eff"] = n_eff,
        _["log_weight_stddev"] = log_wgt_stddevs,
        _["est_k"] = cut_k_values,
        _["step_types"] = step_types,
        _["min_region_cut_sizes"] = min_region_cut_sizes,
        _["max_region_cut_sizes"] = max_region_cut_sizes,
        _["valid_split_region_sizes_list"] = all_steps_valid_split_region_sizes,
        _["valid_region_sizes_to_split_list"] = all_steps_valid_region_sizes_to_split,
        _["merge_split_steps"] = merge_split_step_vec,
        _["merge_split_attempt_counts"] = num_merge_split_attempts_vec,
        _["merge_split_success_mat"] = merge_split_successes_mat
    );

    return out;

}
