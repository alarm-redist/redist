#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/08/18
# Purpose: tidy R wrapper to run gSMC with merge split steps
# redistricting code
####################################################


#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `redist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams generic_redist_gsmc
#' @param manual_k_params Either a single value to use as the splitting parameter for
#' every round or a vector of length num_splitting_steps-1 where each value is the one to use for
#' a split.
#'
#' @return `redist_gsmc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @export
redist_gsmc <- function(
        map, nsims, counties = NULL,
        runs = 1L,
        num_splitting_steps = NULL,
        estimate_cut_k = TRUE,
       manual_k_params = NULL, adapt_k_thresh = .9999,
       split_district_only = FALSE, weight_type = "optimal",
       ms_freq = 0,
       ms_steps_multiplier = 1L,
       run_ms = 0 < ms_freq && ms_freq <= ndists,
       merge_prob_type = "uniform",
       resample = TRUE,
       num_processes=0L, num_threads_per_process=0L,
       multiprocess=FALSE,
       pop_temper = 0,
       init_region_ids_mat = NULL,
       init_region_sizes_mat = NULL,
       min_region_cut_sizes = NULL, max_region_cut_sizes = NULL,
       permitted_split_region_sizes_list = NULL,
       permitted_presplit_region_sizes_list = NULL,
       verbose = FALSE, silent = FALSE, diagnostic_mode = FALSE){

    splitting_params <- list()

    # check its a boolean
    if(!assertthat::is.flag(estimate_cut_k)){
        cli_abort("{.arg estimate_cut_k} must be a Boolean!")
    }
    splitting_params$estimate_cut_k <- estimate_cut_k


    if(estimate_cut_k){
        # check its a scalar
        if (!assertthat::is.number(adapt_k_thresh)){
            cli_abort("{.arg adapt_k_thresh} must be a number")
        }
        # now check its between 0 and 1
        if (adapt_k_thresh < 0 | adapt_k_thresh > 1){
            cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
        }
        splitting_params$adapt_k_thresh <- adapt_k_thresh
    }else{ # else check manual k parameter were passed in
        if(is.null(manual_k_params)){
            cli_abort("If not estimating k for Naive Top K splitting method then {.arg manual_k_params} must be non-null!")
        }
        # check all inputs are numbers
        if(!all(sapply(manual_k_params, assertthat::is.number))) {
            cli_abort("Manual splitting k parameter values all be numbers!")
        }
        # check k param input if not estimating
        if(any(manual_k_params < 1)) {
            cli_abort("Manual splitting k parameter values must be all at least 1.")
        }
        # if just a single number then repeat it
        if(length(manual_k_params) == 1 && floor(manual_k_params) == manual_k_params){
            manual_k_params <- rep(manual_k_params, total_smc_steps)
        }else if(length(manual_k_params) != total_smc_steps){
            cli_abort("K parameter input must be either 1 value or number of smc steps!")
        }else if(any(floor(manual_k_params) != manual_k_params)){
            # if either the length is not ndists-1 or its not all integers then throw
            # error
            cli_abort("K parameter values must be all integers")
        }
        splitting_params$manual_k_params <- manual_k_params
    }



    # figure out the alg type
    if(split_district_only && run_ms){
        alg_type <- "smc_ms"
    }else if(!split_district_only && run_ms){
        alg_type <- "gsmc_ms"
    }else if(split_district_only && !run_ms){
        alg_type <- "basic_smc"
    }else{
        alg_type <- "gsmc"
    }




    generic_redist_gsmc(
        map=map, nsims=nsims, counties = counties, runs = runs,
        alg_name=alg_type,
        split_district_only = split_district_only, weight_type = weight_type,
        sampling_space=GRAPH_PLAN_SPACE_SAMPLING,
        splitting_method=NAIVE_K_SPLITTING,
        splitting_params=splitting_params,
        ms_freq = ms_freq,
        ms_steps_multiplier = ms_steps_multiplier,
        run_ms = run_ms,
        merge_prob_type = merge_prob_type,
        resample = resample,
        num_processes=num_processes, num_threads_per_process=num_threads_per_process,
        multiprocess=multiprocess,
        pop_temper = pop_temper,
        init_region_ids_mat = init_region_ids_mat,
        init_region_sizes_mat = init_region_sizes_mat,
        permitted_split_region_sizes_list = permitted_split_region_sizes_list,
        permitted_presplit_region_sizes_list = permitted_presplit_region_sizes_list,
        num_splitting_steps = num_splitting_steps,
        verbose = verbose, silent = silent, diagnostic_mode = diagnostic_mode)

}

