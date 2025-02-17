#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2025/1/17
# Purpose: tidy R wrapper to run gSMC with merge split steps
# redistricting code on the space of spanning forests
####################################################



#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `treedist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints. It samples plans on the space of spanning
#' forests.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams redist_generic_gsmc
#'
#' @return `treedist_gsmc` returns a [redist_plans] object containing the simulated
#' plans.
#'
#' @export
treedist_gsmc <- function(
        map, nsims, counties = NULL, runs = 1L, alg_name,
        split_district_only = FALSE, weight_type = "optimal",
        splitting_method = UNIF_VALID_EDGE_SPLITTING,
        splitting_method_params = NULL,
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
        permitted_split_region_sizes_list = NULL,
        permitted_presplit_region_sizes_list = NULL,
        num_splitting_steps = NULL,
        verbose = FALSE, silent = FALSE, diagnostic_mode = FALSE
){

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
        sampling_space=FOREST_SPACE_SAMPLING,
        splitting_method=splitting_method,
        splitting_params=splitting_method_params,
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
