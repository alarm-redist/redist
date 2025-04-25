/********************************************************
* Author: Philip O'Sullivan'
* Institution: Harvard University
* Date Created: 2024/10
* Purpose: Helper functions for all redist algorithm types
********************************************************/

#include "redist_alg_helpers.h"


Rcpp::List maximum_input_sizes(){
    // Return results
    List out = List::create(
        _["max_V"] = MAX_SUPPORTED_NUM_VERTICES,
        _["max_districts"] = MAX_SUPPORTED_NUM_DISTRICTS,
        _["max_counties"] = MAX_SUPPORTED_NUM_COUNTIES
    );

    return out;
}


// creates plan ensemble of blank plans
PlanEnsemble::PlanEnsemble(
    int const V, int const ndists, 
    int const total_pop, int const nsims, 
    SamplingSpace const sampling_space,
    RcppThread::ThreadPool &pool
):
    nsims(nsims), 
    V(V),
    flattened_all_plans(V*nsims, 0),
    flattened_all_region_sizes(ndists*nsims, 0),
    flattened_all_region_pops(ndists*nsims, 0),
    flattened_all_region_order_added(ndists*nsims, -1),
    plan_ptr_vec(nsims)
{
    bool const use_graph_space = sampling_space == SamplingSpace::GraphSpace;
    bool const use_forest_space = sampling_space == SamplingSpace::ForestSpace;
    bool const use_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;
    // create the plans 
    Rcpp::Rcout << "Creating Blank Plans!" << std::endl;
    RcppThread::ProgressBar bar(nsims, 1);
    pool.parallelFor(0, nsims, [&] (int i) {
        // create the plan attributes for this specific plan
        PlanVector plan_region_ids(flattened_all_plans, V * i, V * (i+1));
        RegionSizes plan_sizes(flattened_all_region_sizes, ndists * i, ndists * (i+1));
        IntPlanAttribute plan_pops(flattened_all_region_pops, ndists * i, ndists * (i+1));
        IntPlanAttribute plan_region_order_added(flattened_all_region_order_added, ndists * i, ndists * (i+1));
        // create the plans 
        if(use_graph_space){
            plan_ptr_vec[i] = std::make_unique<GraphPlan>(
                ndists, total_pop, 
                plan_region_ids, plan_sizes,
                plan_pops, plan_region_order_added
            );
        }else if(use_forest_space){
            plan_ptr_vec[i] = std::make_unique<ForestPlan>(
                ndists, total_pop, 
                plan_region_ids, plan_sizes,
                plan_pops, plan_region_order_added
            );
        }else if(use_linking_edge_space){
            plan_ptr_vec[i] = std::make_unique<LinkingEdgePlan>(
                ndists, total_pop, 
                plan_region_ids, plan_sizes,
                plan_pops, plan_region_order_added
            );
        }else{
            throw Rcpp::exception("Input is invalid\n");
        }
        ++bar;
    });

    pool.wait();
    
    
}

// creates plan ensemble of partial plans
PlanEnsemble::PlanEnsemble(
    int const V, int const ndists, int const num_regions,
    arma::uvec const &pop, int const nsims,
    SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, 
    Rcpp::IntegerMatrix const &region_sizes_mat,
    RcppThread::ThreadPool &pool 
):    
    nsims(nsims), 
    V(V),
    flattened_all_plans(plans_mat.begin(), plans_mat.end()),
    flattened_all_region_sizes(region_sizes_mat.begin(), region_sizes_mat.end()),
    flattened_all_region_pops(ndists*nsims, 0),
    flattened_all_region_order_added(ndists*nsims, -1),
    plan_ptr_vec(nsims)
{
    // check matrix dimensions 
    if(plans_mat.ncol() != nsims){
        REprintf("The number of columns (%u) in the initial plan matrix was not equal to nsims %d!\n",
            plans_mat.ncol(), nsims
        );
    }
    if(region_sizes_mat.ncol() != nsims){
        REprintf("The number of columns (%u) in the initial sizes matrix  was not equal to nsims %d!\n",
            region_sizes_mat.ncol(), nsims
        );
    }
    if(plans_mat.nrow() != V){
        REprintf("The number of rows (%u) in the initial plan matrix , was not equal to V %d!\n",
            plans_mat.nrow(), V
        );
    }
    if(region_sizes_mat.nrow() != ndists){
        REprintf("The number of rows (%u) in the initial sizes matrix, was not equal to ndists %d!\n",
            region_sizes_mat.nrow(), ndists
        );
    }

    // Now move the data in the matrix 


    bool const use_graph_space = sampling_space == SamplingSpace::GraphSpace;
    bool const use_forest_space = sampling_space == SamplingSpace::ForestSpace;
    bool const use_linking_edge_space = sampling_space == SamplingSpace::LinkingEdgeSpace;

    Rcpp::Rcout << "Loading Partial Plans!" << std::endl;
    RcppThread::ProgressBar bar(nsims, 1);
    pool.parallelFor(0, nsims, [&] (int i) {
        // create the plans 
        if(use_graph_space){
            // create the plan attributes for this specific plan
            PlanVector plan_region_ids(flattened_all_plans, V * i, V * (i+1));
            RegionSizes plan_sizes(flattened_all_region_sizes, ndists * i, ndists * (i+1));
            IntPlanAttribute plan_pops(flattened_all_region_pops, ndists * i, ndists * (i+1));
            IntPlanAttribute plan_region_order_added(flattened_all_region_order_added, ndists * i, ndists * (i+1));
            plan_ptr_vec[i] = std::make_unique<GraphPlan>(
                ndists, num_regions, pop, 
                plan_region_ids, plan_sizes,
                plan_pops, plan_region_order_added
            );

        }else{
            throw Rcpp::exception("This plan type not supported!\n");
        }
        ++bar;
    });

    pool.wait();
}


Rcpp::IntegerMatrix PlanEnsemble::get_R_plans_matrix(){
    // make the plans matrix
    Rcpp::IntegerMatrix plan_mat(V, nsims);
    // copy data over
    std::copy(
        flattened_all_plans.begin(),
        flattened_all_plans.end(),
        plan_mat.begin()
    );
    // now add 1 to everything 
    std::transform(
        plan_mat.begin(), plan_mat.end(), 
        plan_mat.begin(), 
        [](int x) { return x + 1; }
    );
    return plan_mat;
}

// PlanEnsemble::~PlanEnsemble(){
//     Rcpp::Rcout << "gone" << std::endl;
// }

PlanEnsemble get_plan_ensemble(
    int const V, int const ndists, int const num_regions,
    arma::uvec const &pop, int const nsims, 
    SamplingSpace const sampling_space,
    Rcpp::IntegerMatrix const &plans_mat, 
    Rcpp::IntegerMatrix const &region_sizes_mat,
    RcppThread::ThreadPool &pool 
){
    if(num_regions == 1){
        return PlanEnsemble(V, ndists, arma::sum(pop), nsims, sampling_space, pool);
    }else{
        return PlanEnsemble(
            V, ndists, num_regions, pop, nsims,
            sampling_space, plans_mat, region_sizes_mat, 
            pool);
    }
}

//' Copies data from an arma Matrix into an Rcpp Matrix
//'
//' Takes an arma matrix subview and copies all the data into an RcppMatrix
//' of the same size using the Rcpp Threadpool to copy in parallel. 
//'
//'
//' @title Copies data from an arma Matrix into an Rcpp Matrix
//'
//' @param pool A threadpool for multithreading
//' @param arma_mat Subview of an arma unsigned integer matrix 
//' @param rcpp_mat A matrix of integers with the same size as the arma_mat
//'
//' @details Modifications
//'    - The `rcpp_mat` is filled in with the data om the arma matrix subview
//'
//' @noRd
//' @keywords internal
void copy_plans_to_rcpp_mat(
    RcppThread::ThreadPool &pool,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec,
    Rcpp::IntegerMatrix &rcpp_mat,
    bool const copy_sizes_not_ids
){
    int const num_regions = plan_ptrs_vec[0]->num_regions;
    // check dimensions match 
    // number of columns should be size of vector
    // number of rows should match either num regions or V
    if(copy_sizes_not_ids){
        if(rcpp_mat.ncol() != plan_ptrs_vec.size() ||
           rcpp_mat.nrow() != num_regions){
            throw Rcpp::exception("Rcpp Matrix and Plans are not the same size");
        }
    }else{
        if(rcpp_mat.ncol() != plan_ptrs_vec.size() ||
           rcpp_mat.nrow() != plan_ptrs_vec[0]->region_ids.size()){
            throw Rcpp::exception("Rcpp Matrix and Plans are not the same size");
        }
    }

    // go by column because Rcpp matrices are column major
    int ncols = (int) rcpp_mat.ncol();
    

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, ncols, [&] (int j) {
        // Copy plan i into the rcpp matrix
        // either sizes or vertex ids
        if(copy_sizes_not_ids){
            std::copy(
                plan_ptrs_vec[j]->region_sizes.begin(),
                plan_ptrs_vec[j]->region_sizes.begin() + num_regions,
                rcpp_mat.column(j).begin() // Start of column in Rcpp::IntegerMatrix
            );
        }else{
            std::copy(
                plan_ptrs_vec[j]->region_ids.begin(),
                plan_ptrs_vec[j]->region_ids.end(),
                rcpp_mat.column(j).begin() // Start of column in Rcpp::IntegerMatrix
            );
        }
    });

    // Wait for all the threads to finish
    pool.wait();
    
}



//' Reorders all the plans in the vector by order a region was split
//'
//' Takes a vector of plans and uses the vector of dummy plans to reorder
//' each of the plans by the order a region was split.
//'
//'
//' @title Reorders all the plans in the vector by order a region was split
//'
//' @param pool A threadpool for multithreading
//' @param plans_vec A vector of plans
//' @param dummy_plans_vec A vector of dummy plans 
//'
//' @details Modifications
//'    - Each plan in the `plans_vec` object is reordered by when the region was split
//'    - Each plan is a shallow copy of the plans in `plans_vec`
//'
//' @noRd
//' @keywords internal
void reorder_all_plans(
    RcppThread::ThreadPool &pool,
    std::vector<std::unique_ptr<Plan>> &plan_ptrs_vec, 
    std::vector<std::unique_ptr<Plan>> &dummy_plan_ptrs_vec){

    int M = (int) plan_ptrs_vec.size();

    // Parallel thread pool where all objects in memory shared by default
    pool.parallelFor(0, M, [&] (int i) {
        // reorder every plan
        plan_ptrs_vec.at(i)->reorder_plan_by_oldest_split(*dummy_plan_ptrs_vec.at(i));
    });

    // Wait for all the threads to finish
    pool.wait();

    return;
    
}




std::unique_ptr<TreeSplitter> get_tree_splitters(
    MapParams const &map_params,
    SplittingMethodType const splitting_method,
    Rcpp::List const &control,
    int const nsims
){
    int V = map_params.V;
    double target = map_params.target;

    if(splitting_method == SplittingMethodType::NaiveTopK){
        // set splitting k to -1
        return std::make_unique<NaiveTopKSplitter>(V, -1);
    }else if(splitting_method == SplittingMethodType::UnifValid){
        return std::make_unique<UniformValidSplitter>(V);
    }else if(splitting_method == SplittingMethodType::ExpBiggerAbsDev){
        double alpha = as<double>(control["splitting_alpha"]);
        return std::make_unique<ExpoWeightedSplitter>(V, alpha, target);
    }else if(splitting_method == SplittingMethodType::ExpSmallerAbsDev){
        double alpha = as<double>(control["splitting_alpha"]);
        return std::make_unique<ExpoWeightedSmallerDevSplitter>(V, alpha, target);
    }else if(splitting_method == SplittingMethodType::Experimental){
        double epsilon = as<double>(control["splitting_epsilon"]);
        return std::make_unique<ExperimentalSplitter>(V, epsilon, target);
    }else{
        throw Rcpp::exception("Invalid Splitting Method!");
    }
}




SMCDiagnostics::SMCDiagnostics(
    SamplingSpace sampling_space, SplittingMethodType splitting_method_type,
    SplittingSizeScheduleType splitting_schedule_type, 
    std::vector<bool> &merge_split_step_vec,
    int V, int nsims,
    int ndists, int initial_num_regions,
    int total_smc_steps, int total_ms_steps,
    int diagnostic_level,
    bool splitting_all_the_way, bool split_district_only
): diagnostic_level(diagnostic_level), total_steps(total_smc_steps+total_ms_steps),
log_wgt_stddevs(total_smc_steps), acceptance_rates(total_steps),
nunique_parents(total_smc_steps), n_eff(total_smc_steps),
num_merge_split_attempts_vec(total_ms_steps),
cut_k_values(sampling_space == SamplingSpace::GraphSpace ? total_steps : 0)
{
    // Level 1 Diagnostics. Not too big relative to plan size
    log_incremental_weights_mat = arma::dmat(nsims, total_smc_steps, arma::fill::none); // entry [i][s] is the log unnormalized weight of particle i AFTER split s
    draw_tries_mat = Rcpp::IntegerMatrix(nsims, total_steps); // Entry [i][s] is the number of tries it took to form particle i on split s
    parent_index_mat = Rcpp::IntegerMatrix(nsims, total_smc_steps); // Entry [i][s] is the index of the parent of particle i at split s
    // This is a nsims by total_ms_steps matrix where [i][s] is the number of 
    // successful merge splits performed for plan i on merge split round s
    merge_split_successes_mat = Rcpp::IntegerMatrix(nsims, total_ms_steps);
    // counts the size of the trees
    tree_sizes_mat = Rcpp::IntegerMatrix(ndists, total_steps);
    successful_tree_sizes_mat = Rcpp::IntegerMatrix(ndists, total_steps);


    // Level 2
    parent_unsuccessful_tries_mat = Rcpp::IntegerMatrix(nsims, total_smc_steps);


    bool diagnostic_mode = diagnostic_level == 1;
    // level 3
    all_steps_plan_region_ids_list.reserve(diagnostic_mode ? total_steps : 0);
    all_steps_forests_adj_list.resize(
        (diagnostic_mode && sampling_space != SamplingSpace::GraphSpace) ? total_steps : 0
    );
    all_steps_linking_edge_list.resize(
        (diagnostic_mode && sampling_space == SamplingSpace::LinkingEdgeSpace) ? total_steps : 0
    );
    all_steps_valid_region_sizes_to_split.resize(diagnostic_mode ? total_smc_steps : 0);
    all_steps_valid_split_region_sizes.resize(diagnostic_mode ? total_smc_steps : 0);


    // Store size at every step but last one if needed
    int plan_dval_list_size = (diagnostic_mode & !split_district_only) ? total_steps-1 : 0;
    if(!splitting_all_the_way) plan_dval_list_size++;
    
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
}


void SMCDiagnostics::add_full_step_diagnostics(
    int const total_steps, bool const splitting_all_the_way,
    int const step_num, int const merge_split_step_num, int const smc_step_num,
    bool const is_smc_step,
    SamplingSpace const sampling_space,
    RcppThread::ThreadPool &pool,
    std::vector<std::unique_ptr<Plan>> &plans_ptr_vec, 
    std::vector<std::unique_ptr<Plan>> &new_plans_ptr_vec,
    SplittingSchedule const &splitting_schedule
){
    //if(diagnostic_mode){ // record if in diagnostic mode and generalized splits
    // reorder the plans by oldest split if either we'vxe done any merge split or
    // its generalized region splits 

    bool const split_district_only = splitting_schedule.schedule_type == SplittingSizeScheduleType::DistrictOnly;
    int const nsims = plans_ptr_vec.size();

    // if smc step update splitting step info
    if(is_smc_step){
        int current_num_regions = plans_ptr_vec[0]->num_regions;
        // save the acceptable split sizes 
        for (int region_size = 1; region_size <= splitting_schedule.ndists - current_num_regions + 2; region_size++)
        {
            if(splitting_schedule.valid_split_region_sizes[region_size]){
                all_steps_valid_split_region_sizes[smc_step_num].push_back(region_size);
            }
            if(splitting_schedule.valid_region_sizes_to_split[region_size]){
                
                all_steps_valid_region_sizes_to_split[smc_step_num].push_back(region_size);;
            }
        }
    }

    if(merge_split_step_num > 0 || !split_district_only){
        reorder_all_plans(pool, plans_ptr_vec, new_plans_ptr_vec);
    }

    // Copy the vertex plan matrix 
    copy_plans_to_rcpp_mat(pool, 
        plans_ptr_vec, 
        all_steps_plan_region_ids_list.at(step_num),
        false);

    // store the 
    if(!(sampling_space == SamplingSpace::GraphSpace)){
        all_steps_forests_adj_list.at(step_num).reserve(nsims);
        for (size_t i = 0; i < nsims; i++)
        {
            // add the forests from each plan at this step
            all_steps_forests_adj_list.at(step_num).push_back(
                plans_ptr_vec.at(i)->get_forest_adj()
            );
        }
        if(sampling_space == SamplingSpace::LinkingEdgeSpace){
            for (size_t i = 0; i < nsims; i++)
            {
                // add the forests from each plan at this step
                all_steps_linking_edge_list.at(step_num).push_back(
                    plans_ptr_vec.at(i)->get_linking_edges()
                );
            }
            
        }
    }
    
    // Copy the sizes if neccesary 
    if(!split_district_only && (step_num < total_steps-1 || !splitting_all_the_way)){
        copy_plans_to_rcpp_mat(
            pool, 
            plans_ptr_vec, 
            region_sizes_mat_list.at(step_num),
            true);
    }

    return;

}


void SMCDiagnostics::add_diagnostics_to_out_list(Rcpp::List &out){
    // make parent index 1 indexed in place
    std::transform(
        parent_index_mat.begin(), parent_index_mat.end(), 
        parent_index_mat.begin(), 
        [](int x) { return x + 1; }
    );

    out["acceptance_rates"] = acceptance_rates;
    out["draw_tries_mat"] = draw_tries_mat;
    out["parent_index"] = parent_index_mat;
    out["parent_unsuccessful_tries_mat"] = parent_unsuccessful_tries_mat;
    out["step_n_eff"] = n_eff;
    out["nunique_parent_indices"] = nunique_parents;
    out["tree_sizes"] = tree_sizes_mat;
    out["successful_tree_sizes"] = successful_tree_sizes_mat;
    out["log_weight_stddev"] = log_wgt_stddevs;
    out["est_k"] = cut_k_values;
    out["log_incremental_weights_mat"] = log_incremental_weights_mat;
    out["merge_split_attempt_counts"] = num_merge_split_attempts_vec;
    out["merge_split_success_mat"] = merge_split_successes_mat;
    out["region_ids_mat_list"] = all_steps_plan_region_ids_list;
    out["region_sizes_mat_list"] = region_sizes_mat_list;
    out["forest_adjs_list"] = all_steps_forests_adj_list;
    out["linking_edges_list"] = all_steps_linking_edge_list;
    out["valid_split_region_sizes_list"] = all_steps_valid_split_region_sizes;
    out["valid_region_sizes_to_split_list"] = all_steps_valid_region_sizes_to_split;

    return;
}


// 
Rcpp::IntegerVector resample_plans_lowvar(
    Rcpp::NumericVector const &normalized_weights,
    Rcpp::IntegerMatrix &plans_mat,
    Rcpp::IntegerMatrix &region_sizes_mat,
    bool const reorder_sizes_mat
){
    // generate resampling index
    int const N = normalized_weights.size();

    int rng_seed = (int) Rcpp::sample(INT_MAX, 1)[0];
    RNGState rng_state(rng_seed, 42);
    double r = rng_state.r_unif() / N;
    double cuml = normalized_weights[0];
    Rcpp::IntegerVector resample_index(N);
    std::vector<bool> index_unchanged(N);

    int i = 0;
    for (int n = 0; n < N; n++) {
        double u = r + n / (double) N;
        while (u > cuml) {
            cuml += normalized_weights[++i]; // increment then access
        }
        // resample_index maps entry i to its new value 
        // `resample_index[i] = k` means you should replace plan i with plan k
        resample_index[n] = i;
    }
    

    // Now we're going to reorder things one row at a time 
    // makes algout$plans[i] now equal to algout$plans[rs_idx[i]]
    // so we're mapping i -> rs_idx[i]
    std::vector<int> buffer(N);
    int const nrow = plans_mat.nrow();
    for (int row = 0; row < nrow; ++row) {
        // copy current row
        for (int col = 0; col < N; ++col) {
          buffer[col] = plans_mat(row, col);
        }
        // write back in new order
        for (int col = 0; col < N; ++col) {
            plans_mat(row, col) = buffer[resample_index[col]];
        }
    }

    // Reorder region sizes if needed 
    if(reorder_sizes_mat){
        std::vector<int> sizes_buffer(N);
        int const nrow = region_sizes_mat.nrow();
        for (int row = 0; row < nrow; ++row) {
            // copy current row
            for (int col = 0; col < N; ++col) {
              sizes_buffer[col] = region_sizes_mat(row, col);
            }
            // write back in new order
            for (int col = 0; col < N; ++col) {
                region_sizes_mat(row, col) = sizes_buffer[resample_index[col]];
            }
        }
    }

    // make resampling thing 1 indexed
    std::transform(
        resample_index.begin(), resample_index.end(), 
        resample_index.begin(), 
        [](int x) { return x + 1; }
    );

    return resample_index;
}