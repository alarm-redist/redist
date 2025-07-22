/********************************************************
 * Author: Cory McCartan
 * Institution: Harvard University
 * Date Created: 2021/02
 * Purpose: Merge-split MCMC redistricting sampler
 * (like Carter et al. 2019 but using the SMC proposal)
 ********************************************************/

#include "merge_split.h"

constexpr bool DEBUG_PURE_MS_VERBOSE = false; // Compile-time constant


/*
 * Main entry point.
 *
 * USING MCMC
 * Sample `N` redistricting plans on map `g`, ensuring that the maximum
 * population deviation is between `lower` and `upper` (and ideally `target`)
 */
Rcpp::List ms_plans(
    int const nsims, int const warmup, int const thin, 
    int const ndists, int const total_seats, Rcpp::IntegerVector const &district_seat_sizes, 
    List const &adj_list, const arma::uvec &counties, const arma::uvec &pop,
    double const target, double const lower, double const upper,
    double const rho, // compactness 
    Rcpp::IntegerMatrix const &init_plan, Rcpp::IntegerMatrix const &init_seats,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &merge_prob_type, // method for setting probability of picking a pair to merge
    List const &control, // control has pop temper, and k parameter value, and whether only district splits are allowed
    List const &constraints, // constraints 
    int const verbosity, bool const diagnostic_mode
) {
    // whether or not to perform MH step
    bool do_mh = (bool) control["do_mh"];

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 1!\n");
    Rcpp::List out; // return 
    // re-seed MT so that `set.seed()` works in R
    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    RNGState rng_state(global_rng_seed);
    // Set the sampling space 
    SamplingSpace sampling_space = get_sampling_space(sampling_space_str);
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // TODO: Legacy, in future remove
    global_seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // make sure both are 1 column matrix
    if(init_plan.ncol() > 1 || init_seats.ncol() > 1){
        throw Rcpp::exception("Error!\n");
    }

    // Create map level graph and county level multigraph
    MapParams const map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper
    );
    int V = map_params.g.size();


    int initial_num_regions = static_cast<int>(ndists);
    if(init_seats.nrow() != initial_num_regions){
        REprintf("Inferred %d Initial Regions but Region Sizes Matrix has %u columns!\n",
            initial_num_regions, init_seats.nrow());
        throw Rcpp::exception("Initial plan region sizes should match number of regions");
    }
    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 2!\n");
    
    // get splitting type 
    SplittingMethodType splitting_method = get_splitting_type(
        static_cast<std::string>(control["splitting_method"])
        );


    // get the splitting size regime
    SplittingSizeScheduleType splitting_size_regime = SplittingSizeScheduleType::PureMergeSplitSize;

    auto splitting_schedule_ptr = std::make_unique<PureMSSplittingSchedule>(ndists, total_seats, as<std::vector<int>>(district_seat_sizes));
    // splitting_schedule_ptr->print_current_step_splitting_info();
    bool const mmd_plans = map_params.is_mmd;


    // Add scoring function (constraints)
    // No population tempering
    ScoringFunction const scoring_function(
        map_params, constraints, 0, false
    );



    // Now create diagnostic information 
    Rcpp::NumericVector log_mh_ratios(nsims); // stores log mh ratio 
    Rcpp::IntegerMatrix saved_plans_mat(V, nsims);
    Rcpp::IntegerMatrix saved_district_pops_mat(ndists, nsims);
    Rcpp::IntegerMatrix saved_plan_sizes(
        mmd_plans ? ndists : 1,
        mmd_plans ? nsims : 1
    );
    int current_plan_mat_col = 0;
    std::vector<int> tree_sizes(total_seats, 0);
    std::vector<int> successful_tree_sizes(total_seats, 0);

    // Level 3 
    // Saves proposal plans 
    Rcpp::IntegerMatrix proposed_plans_mat(
        diagnostic_mode ? V : 1, 
        diagnostic_mode ? nsims: 1);

    std::vector<std::vector<Graph>> all_steps_forests_adj_list;
    all_steps_forests_adj_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::ForestSpace) ? nsims : 0
    );

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 3!\n");

    // Create current plan and proposal

    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);


    

    Rcpp::IntegerVector mh_decisions(nsims);
    double mha;

    int total_post_warmup_steps = nsims * thin;
    int total_steps = total_post_warmup_steps + warmup;
    int start = 1 - warmup;

    // Track the total number of successes during the warmup 
    int warmup_acceptances = 0;
    // Track total number of successes after warmup
    int post_warump_acceptances = 0;

    {
    USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
    PlanMultigraph current_plan_multigraph(map_params);
    PlanMultigraph proposed_plan_multigraph(map_params);



    RcppThread::ThreadPool pool(1);
    // underlying vector from plan    
    PlanEnsemble plan_ensemble = get_plan_ensemble(
        map_params, *splitting_schedule_ptr,
        initial_num_regions, 
        1, sampling_space,
        init_plan, init_seats,
        rng_states, pool, verbosity
    );
    // plan_ensemble.plan_ptr_vec[0]->Rprint(true);
    // now get for proposal plan
    PlanEnsemble proposal_plan_ensemble = get_plan_ensemble(
        map_params, *splitting_schedule_ptr,
        initial_num_regions, 
        1, sampling_space,
        init_plan, init_seats,
        rng_states, pool, verbosity
    );

    


    // splitter
    std::unique_ptr<TreeSplitter> tree_splitter_ptr = get_tree_splitters(
        map_params, splitting_method, control, nsims
    );
    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 4!\n");
    // sanity check make sure the plan is ok 
    std::vector<bool> county_component_lookup(ndists * map_params.num_counties, false);
    bool hierarchically_valid = current_plan_multigraph.is_hierarchically_valid(
            *plan_ensemble.plan_ptr_vec[0], county_component_lookup
    );

    if(!hierarchically_valid){
        plan_ensemble.plan_ptr_vec[0]->Rprint(true);
        throw Rcpp::exception("Initial Plan is not hierarchically valid. Either turn off counties or pass in a hierarchically valid plan\n");
    }
    // build multigraph on current plan and get pairs of adj districts to merge
    auto build_result = plan_ensemble.plan_ptr_vec[0]->attempt_to_get_valid_mergesplit_pairs(
        current_plan_multigraph, *splitting_schedule_ptr, scoring_function
    );
    // shouldn't be possible but just a sanity check
    if (!build_result.first){
        throw Rcpp::exception("BIG ERROR: Plan registered as hierarchically valid but we failed to build hierarchical multigraph!\n");
    }
    auto current_plan_adj_region_pairs = plan_ensemble.plan_ptr_vec[0]->attempt_to_get_valid_mergesplit_pairs(
        current_plan_multigraph, *splitting_schedule_ptr, scoring_function
    ).second;

    // get weights 
    arma::vec current_plan_pair_unnoramalized_wgts = get_adj_pair_unnormalized_weights(
        *plan_ensemble.plan_ptr_vec[0],
        current_plan_adj_region_pairs,
        merge_prob_type
    );


    // Set or estimate k if doing graph space sampling
    if(sampling_space == SamplingSpace::GraphSpace){
        int cut_k;
        bool try_to_estimate_cut_k = as<bool>(control["estimate_cut_k"]);
        if(try_to_estimate_cut_k){
            double thresh = (double) control["adapt_k_thresh"];
            double tol = std::max(target - lower, upper - target) / target;

            cut_k = estimate_mergesplit_cut_k(
                    *plan_ensemble.plan_ptr_vec[0], current_plan_multigraph,  
                    *splitting_schedule_ptr,
                    thresh, tol, rng_state
                );

            if(verbosity >= 3){
                Rcout << " Using estimated k = " << cut_k << "\n";
            }
        }else{
            cut_k = as<int>(control["manual_k"]);
            if (verbosity >= 3){
                Rcout << "Using k = " << cut_k << "\n";
            }
        }
        // update the tree splitter
        tree_splitter_ptr->update_single_int_param(cut_k);
        out["cut_k"] = cut_k;
    }

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 5!\n");


    if(DEBUG_PURE_MS_VERBOSE){
        Rprintf("Adj regions are:");
        for (auto const &a_pair: current_plan_adj_region_pairs)
        {
            Rprintf("(%d, %d), ", a_pair.first, a_pair.second);
        }
        Rprintf("\n");
    }

    

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 6!\n");
    // Loading Info
    if (verbosity >= 1) {
        Rcout.imbue(std::locale::classic());
        Rcout << std::fixed << std::setprecision(0);
        Rcout << "MERGE SPLIT MONTE CARLO" << std::endl;
        Rcout << "Using " << sampling_space_to_str(sampling_space);
        Rcout << " Sampling space to sample " << nsims << " " << V << "-unit ";
        Rcout << "maps with " << ndists << " districts and population between "
              << lower << " and " << upper;
        Rcout << " Using " << splitting_method_to_str(splitting_method) << "!"
              << std::endl;
        if (map_params.cg.size() > 1){
            Rcout << "Sampling hierarchically with respect to the "
                << map_params.cg.size() << " administrative units." << std::endl;
        }
        if(scoring_function.total_soft_constraints > 0){
            Rcout << "Applying " << scoring_function.total_soft_constraints << " constraints"
            << std::endl;
        }
    }
    

    RObject bar = cli_progress_bar(total_steps, cli_config(false));
    for (int i = start, step_num = 0 ; i <= total_post_warmup_steps; i++, step_num++) {
        // Index 0 or less is warmup
        bool in_warmup = i <= 0;
        if(DEBUG_PURE_MS_VERBOSE){
            Rprintf("Iter %d and idx %d \n", i, current_plan_mat_col);
        }

        // attempt to mergesplit 
        std::tuple<bool, bool, double, int> mergesplit_result = attempt_mergesplit_step(
            map_params, *splitting_schedule_ptr, scoring_function,
            rng_state, sampling_space,
            *plan_ensemble.plan_ptr_vec[0], *proposal_plan_ensemble.plan_ptr_vec[0], 
            ust_sampler, *tree_splitter_ptr,
            current_plan_multigraph, proposed_plan_multigraph,
            merge_prob_type, save_edge_selection_prob,
            current_plan_adj_region_pairs,
            current_plan_pair_unnoramalized_wgts,
            rho, true, do_mh
        );
        // count size
        ++tree_sizes[std::get<3>(mergesplit_result)-1];

        // copy if needed
        if(!in_warmup && i % thin == 0){
            // Copy the plan into the matrix 
            // since a Plan's region IDs are associated with 
            std::copy(
                plan_ensemble.flattened_all_plans.begin(), 
                plan_ensemble.flattened_all_plans.end(),
                saved_plans_mat.column(current_plan_mat_col).begin() // Start of column in Rcpp::IntegerMatrix
            );
            // copy region populations
            std::copy(
                plan_ensemble.flattened_all_region_pops.begin(), 
                plan_ensemble.flattened_all_region_pops.end(),
                saved_district_pops_mat.column(current_plan_mat_col).begin() // Start of column in Rcpp::IntegerMatrix
            );
            // if mmd copy sizes
            if(mmd_plans){
                std::copy(
                    plan_ensemble.flattened_all_region_sizes.begin(), 
                    plan_ensemble.flattened_all_region_sizes.end(),
                    saved_plan_sizes.column(current_plan_mat_col).begin() // Start of column in Rcpp::IntegerMatrix
                );
            }
            if(diagnostic_mode){
                std::copy(
                    proposal_plan_ensemble.plan_ptr_vec[0]->region_ids.begin(), 
                    proposal_plan_ensemble.plan_ptr_vec[0]->region_ids.end(),
                    proposed_plans_mat.column(current_plan_mat_col).begin() // Start of column in Rcpp::IntegerMatrix
                );
            }
            
            // save ratio
            log_mh_ratios(current_plan_mat_col) = std::get<2>(mergesplit_result);
            mh_decisions(current_plan_mat_col) = std::get<1>(mergesplit_result);

            if(DEBUG_PURE_MS_VERBOSE) Rprintf("log_mh_ratio = %f\n", std::get<2>(mergesplit_result));

            // increase the count of columne we're on
            current_plan_mat_col++;
        }

        // if successful then update acceptance count
        if(std::get<1>(mergesplit_result)){
            ++successful_tree_sizes[std::get<3>(mergesplit_result)-1];
            if(in_warmup){
                ++warmup_acceptances;
            }else{
                ++post_warump_acceptances;
            }
        }

        if (verbosity >= 1 && CLI_SHOULD_TICK) {
            cli_progress_set(bar, step_num);
            mha = static_cast<double>(warmup_acceptances + post_warump_acceptances) / (step_num+1);
            cli_progress_set_format(bar, "{cli::pb_bar} {cli::pb_percent} | ETA: {cli::pb_eta} | MH Acceptance: %.2f", mha);
        }
        Rcpp::checkUserInterrupt();
 
    }
    cli_progress_done(bar);
    if(DEBUG_PURE_MS_VERBOSE) REprintf("Done with main MCMC loop\n");
    }
    
    if (verbosity >= 1) {
        Rcout << "Acceptance rate: " << std::setprecision(2) << 
        (100.0 * (warmup_acceptances + post_warump_acceptances)) / (total_steps) << 
        "%" << std::endl;
    }

    // now add 1 to plans 
    std::transform(
        saved_plans_mat.begin(), saved_plans_mat.end(), 
        saved_plans_mat.begin(), 
        [](int x) { return x + 1; }
    );

    if(DEBUG_PURE_MS_VERBOSE) REprintf("Added one to plans, now creating diagnostic list.\n");
    
    out["plans"] = saved_plans_mat;
    out["region_pops"] = saved_district_pops_mat;
    out["seats"] = saved_plan_sizes;
    out["mhdecisions"] = mh_decisions;
    out["total_steps"] = total_steps;
    out["warmup_acceptances"] = warmup_acceptances;
    out["post_warump_acceptances"] = post_warump_acceptances;
    out["log_mh_ratio"] = log_mh_ratios;
    out["tree_sizes"] = tree_sizes;

    if(diagnostic_mode){
        // now add 1 to plans 
        std::transform(
            proposed_plans_mat.begin(), proposed_plans_mat.end(), 
            proposed_plans_mat.begin(), 
            [](int x) { return x + 1; }
        );
        out["proposed_plans"] = proposed_plans_mat;
    }

    if(DEBUG_PURE_MS_VERBOSE) REprintf("Done now returning!\n");

    return out;
}
