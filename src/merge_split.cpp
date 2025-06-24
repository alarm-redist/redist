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
 * Select a pair of neighboring districts i, j
 */
void select_pair(int n_distr, const Graph &dist_g, int &i, int &j, RNGState &rng_state) {
    i = rng_state.r_int(n_distr);
    std::vector<int> nbors = dist_g[i];
    j = nbors[rng_state.r_int(nbors.size())];
}

/*
 * Choose k and multiplier for efficient, accurate sampling
 */
int adapt_ms_parameters(const Graph &g, int n_distr, double thresh,
                         double tol, PlanVector const &region_ids, const uvec &counties,
                         Multigraph const &cg, const uvec &pop, double target,
                         RNGState &rng_state) {
    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Estimation checkpint 1!");
    int k;
    // sample some spanning trees and compute deviances
    int V = g.size();
    // IN FUTURE USE MY OWN FUNCTION THIS HAS INEXING ERRORS
    Graph dist_g = district_graph(g, region_ids, n_distr, true);
    int k_max = std::min(20 + ((int) std::sqrt(V)), V - 1); // heuristic
    int N_adapt = (int) std::floor(4000.0 / sqrt((double) V));
    

    double lower = target * (1 - tol);
    double upper = target * (1 + tol);
    
    std::vector<std::vector<double>> devs;
    vec distr_ok(k_max+1, fill::zeros);
    int root;
    int max_ok = 0;
    std::vector<bool> ignore(V);
    std::vector<bool> visited(V);
    int distr_1, distr_2;
    int max_V = 0;
    Tree ust = init_tree(V);
    for (int i = 0; i < N_adapt; i++) {
        double joint_pop = 0;
        select_pair(n_distr, dist_g, distr_1, distr_2, rng_state);
        int n_vtx = 0;
        for (int j = 0; j < V; j++) {
            if (region_ids[j] == distr_1 || region_ids[j] == distr_2) {
                joint_pop += pop(j);
                ignore[j] = false;
                n_vtx++;
            } else {
                ignore[j] = true;
            }
        }
        if (n_vtx > max_V) max_V = n_vtx;

        clear_tree(ust);
        int result = sample_sub_ust(g, ust, V, root, visited, ignore,
                                    pop, lower, upper, counties, cg,
                                    rng_state);
        if (result != 0) {
            i--;
            continue;
        }

        devs.push_back(tree_dev(ust, root, pop, joint_pop, target));
        int n_ok = 0;
        for (int j = 0; j < V-1; j++) {
            if (ignore[j]) devs.at(i).at(j) = 2; // force not to work
            n_ok += devs.at(i).at(j) <= tol;
        }

        if (n_ok <= k_max)
            distr_ok(n_ok) += 1.0 / N_adapt;
        if (n_ok > max_ok && n_ok < k_max)
            max_ok = n_ok;
    }

    // For each k, compute pr(selected edge within top k),
    // among maps where valid edge was selected
    uvec idxs(N_adapt);
    for (k = 1; k <= k_max; k++) {
        idxs = as<uvec>(Rcpp::sample(k, N_adapt, true, R_NilValue, false));
        double sum_within = 0;
        int n_ok = 0;
        for (int i = 0; i < N_adapt; i++) {
            double dev = devs.at(i).at(idxs[i]);
            if (dev > tol) continue;
            else n_ok++;
            for (int j = 0; j < N_adapt; j++) {
                sum_within += ((double) (dev <= devs.at(j).at(k-1))) / N_adapt;
            }
        }
        if (sum_within / n_ok >= thresh) break;
    }

    if (k == k_max + 1) {
        Rcerr << "Warning: maximum hit; falling back to naive k estimator.\n";
        k = max_ok + 1;
    }

    k = std::min(k, max_V - 1);
    return(k);
}




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
    Rcpp::IntegerMatrix const &initial_plan, Rcpp::IntegerMatrix const &initial_region_sizes,
    std::string const &sampling_space_str, // sampling space (graphs, forest, etc)
    std::string const &merge_prob_type, // method for setting probability of picking a pair to merge
    List const &control, // control has pop temper, and k parameter value, and whether only district splits are allowed
    List const &constraints, // constraints 
    int const verbosity, bool const diagnostic_mode
) {
    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 1!\n");
    Rcpp::List out; // return 
    // re-seed MT so that `set.seed()` works in R
    int global_rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    RNGState rng_state(global_rng_seed, 42);
    // Set the sampling space 
    SamplingSpace sampling_space = get_sampling_space(sampling_space_str);
    bool save_edge_selection_prob = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // TODO: Legacy, in future remove
    global_seed_rng((int) Rcpp::sample(INT_MAX, 1)[0]);

    // make sure both are 1 column matrix
    if(initial_plan.ncol() > 1 || initial_region_sizes.ncol() > 1){
        throw Rcpp::exception("Error!\n");
    }


    int initial_num_regions = static_cast<int>(ndists);
    if(initial_region_sizes.nrow() != initial_num_regions){
        REprintf("Inferred %d Initial Regions but Region Sizes Matrix has %u columns!\n",
            initial_num_regions, initial_region_sizes.nrow());
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


    // Do some input checking 

    // Create map level graph and county level multigraph
    MapParams const map_params(
        adj_list, counties, pop, 
        ndists, total_seats, as<std::vector<int>>(district_seat_sizes),
        lower, target, upper
    );
    int V = map_params.g.size();

    // Add scoring function (constraints)
    // No population tempering
    ScoringFunction const scoring_function(
        map_params, constraints, 0
    );



    // Now create diagnostic information 
    Rcpp::NumericVector log_mh_ratios(nsims); // stores log mh ratio 
    Rcpp::IntegerMatrix saved_plans_mat(V, nsims);
    int current_plan_mat_col = 0;
    std::vector<int> tree_sizes(ndists, 0);
    std::vector<int> successful_tree_sizes(ndists, 0);

    // Level 3 
    // Saves proposal plans 
    Rcpp::IntegerMatrix proposed_plans_mat(
        diagnostic_mode ? V : 0, 
        diagnostic_mode ? nsims: 0);

    std::vector<std::vector<Graph>> all_steps_forests_adj_list;
    all_steps_forests_adj_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::ForestSpace) ? nsims : 0
    );

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 3!\n");


    bool split_district_only;

    // Create current plan and proposal

    int global_rng_seed2 = (int) Rcpp::sample(INT_MAX, 1)[0];
    std::vector<RNGState> rng_states;rng_states.reserve(1);
    rng_states.emplace_back(global_rng_seed2, 6);


    RcppThread::ThreadPool pool(1);
    // underlying vector from plan    
    PlanEnsemble plan_ensemble = get_plan_ensemble(
        map_params, initial_num_regions,
        1, sampling_space,
        initial_plan, initial_region_sizes,
        rng_states, pool, verbosity
    );
    // now get for proposal plan
    PlanEnsemble proposal_plan_ensemble = get_plan_ensemble(
        map_params, initial_num_regions,
        1, sampling_space,
        initial_plan, initial_region_sizes,
        rng_states, pool, verbosity
    );


    // splitter
    std::unique_ptr<TreeSplitter> tree_splitter_ptr = get_tree_splitters(
        map_params, splitting_method, control, nsims
    );
    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 4!\n");

    // Set or estimate k if doing graph space sampling
    if(sampling_space == SamplingSpace::GraphSpace){
        int cut_k;
        bool try_to_estimate_cut_k = as<bool>(control["estimate_cut_k"]);
        if(try_to_estimate_cut_k){
            double thresh = (double) control["adapt_k_thresh"];
            double tol = std::max(target - lower, upper - target) / target;
            cut_k = adapt_ms_parameters(
                map_params.g, ndists, thresh, tol, 
                plan_ensemble.plan_ptr_vec[0]->region_ids, 
                counties, map_params.cg, pop, target, rng_state);
            if(verbosity >= 3){
                Rcout << " Using estimated k = " << cut_k << ")\n";
            }
        }else{
            cut_k = as<int>(control["manual_k"]);
            if (verbosity >= 3){
                Rcout << "Using k = " << cut_k << "\n";
            }
        }
        // update the tree splitter
        tree_splitter_ptr->update_single_int_param(cut_k);
        out["est_k"] = cut_k;
    }

    if(DEBUG_PURE_MS_VERBOSE) Rprintf("Checkpoint 5!\n");


    Rcpp::IntegerVector mh_decisions(nsims);
    double mha;

    int total_post_warmup_steps = nsims * thin;
    int total_steps = total_post_warmup_steps + warmup;
    int start = 1 - warmup;

    // Track the total number of successes during the warmup 
    int warmup_acceptances = 0;
    // Track total number of successes after warmup
    int post_warump_acceptances = 0;

    USTSampler ust_sampler(map_params, *splitting_schedule_ptr);
    PlanMultigraph current_plan_multigraph(map_params);
    PlanMultigraph proposed_plan_multigraph(map_params);

    // build multigraph on current plan and get pairs of adj districts to merge
    auto current_plan_adj_region_pairs = plan_ensemble.plan_ptr_vec[0]->attempt_to_get_valid_mergesplit_pairs(
        current_plan_multigraph, *splitting_schedule_ptr
    ).second;
    // get weights 
    arma::vec current_plan_pair_unnoramalized_wgts = get_adj_pair_unnormalized_weights(
        *plan_ensemble.plan_ptr_vec[0],
        current_plan_adj_region_pairs,
        merge_prob_type
    );

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
        Rcout << std::fixed << std::setprecision(4);
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
            rho, true
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
    
    out["plans"] = saved_plans_mat;
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

    return out;
}
