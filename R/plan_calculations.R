#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2025/4/2
# Purpose: Functions for computing unnormalized target density of plans
# and optimal weight functions
####################################################


#' Computes the unnormalized log target density for a ONE INDEXED! plans matrix
#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `generic_redist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints.
#'
#' Sampling can be performed either on the space of plans (graph partitions) or
#' on the space of spanning forests. This is an internal function designed to
#' support either sampling space and any arbitrary tree splitting method. The
#' functions `redist_gsmc` and `treedist_gsmc` are largely wrappers to calling
#' this function.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams generic_redist_gsmc
#' @param num_threads The number of threads used for computing
#' the target densities.
#' @param plan_matrix A matrix of 1-indexed plans
#' @param init_region_sizes_mat A matrix of region sizes. The matrix must have
#' the same number of columns as `plan_matrix`, `ndists` rows, and each column
#' must sum to `ndists`.
#'
#' @return A vector of the unnormalized target density of each plan
#'
#' @md
#' @order 1
#' @export
compute_log_target_density <- function(
        map, plan_matrix, sizes_matrix = NULL,
        counties = NULL, compactness = 1L,
        constraints = list(),
        num_threads = 0, pop_temper = 0L
        ){
    ndists <- attr(map, "ndists")
    # if its just a single plan make it a matrix
    if(is.vector(plan_matrix)){
        plan_matrix <- as.matrix(plan_matrix, cols = 1)
    }
    # if sizes matrix is null then assume all regions are the same size
    if(is.null(sizes_matrix)){
        sizes_matrix <- matrix(1L, nrow = ndists, ncol = ncol(plan_matrix))
    }


    counties_q <- rlang::enquo(counties)
    # get validated inputs
    map_params <- get_map_parameters(map, counties_q=counties_q, use_counties_q = TRUE)
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    # need to pass in quosure
    constraints_q <- rlang::enquo(constraints)
    constraints <- validate_constraints(map, constraints_q=constraints_q, use_constraints_q=TRUE)


    num_regions <- dplyr::n_distinct(plan_matrix[,1])
    unnormalized_log_density <- compute_log_unnormalized_plan_target_density(
        adj_list, counties, pop,
        constraints, pop_temper, rho=compactness,
        ndists=ndists, num_regions=num_regions,
        lower=pop_bounds[1],
        target=pop_bounds[2],
        upper=pop_bounds[3],
        region_ids=plan_matrix-1L,
        region_sizes=sizes_matrix,
        num_threads=num_threads
    )

    # check if anything was -Inf so probability 0
    num_prob_zero <- sum(!is.finite(unnormalized_log_density))

    if(num_prob_zero > 0){
        cli::cli_warn("{num_prob_zero} of the {length(unnormalized_log_density)} plans have probability 0!")
    }

    return(unnormalized_log_density)
}




#' Computes the unnormalized log target density for a ONE INDEXED! plans matrix
#' gSMC Redistricting Sampler (O'Sullivan, McCartan and Imai ???)
#'
#' `generic_redist_gsmc` uses a Sequential Monte Carlo algorithm (O'Sullivan, McCartan and Imai ???)
#' to generate representative samples of congressional or legislative
#' redistricting plans according to contiguity, population, compactness, and
#' administrative boundary constraints.
#'
#' Sampling can be performed either on the space of plans (graph partitions) or
#' on the space of spanning forests. This is an internal function designed to
#' support either sampling space and any arbitrary tree splitting method. The
#' functions `redist_gsmc` and `treedist_gsmc` are largely wrappers to calling
#' this function.
#'
#' This function draws samples from a specific target measure controlled by
#' the `map` parameters.
#'
#' @inheritParams generic_redist_gsmc
#' @param num_threads The number of threads used for computing
#' the target densities.
#' @param plan_matrix A matrix of 1-indexed plans
#' @param init_region_sizes_mat A matrix of region sizes. The matrix must have
#' the same number of columns as `plan_matrix`, `ndists` rows, and each column
#' must sum to `ndists`.
#'
#' @return A vector of the unnormalized target density of each plan
#'
#' @md
#' @order 1
#' @export
compute_log_optimal_weights <- function(
        map, plan_matrix, sizes_matrix = NULL,
        counties = NULL, constraints = list(),
        splitting_schedule = "any_valid_sizes",
        num_threads = 0, compactness = 1L, pop_temper = 0L
){
    ndists <- attr(map, "ndists")
    # if its just a single plan make it a matrix
    if(is.vector(plan_matrix)){
        plan_matrix <- as.matrix(plan_matrix, cols = 1)
    }
    # if sizes matrix is null then assume all regions are the same size
    if(is.null(sizes_matrix)){
        sizes_matrix <- matrix(1L, nrow = ndists, ncol = ncol(plan_matrix))
    }


    counties_q <- rlang::enquo(counties)
    # get validated inputs
    map_params <- get_map_parameters(map, counties_q=counties_q, use_counties_q = TRUE)
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    # need to pass in quosure
    constraints_q <- rlang::enquo(constraints)
    constraints <- validate_constraints(map, constraints_q=constraints_q, use_constraints_q=TRUE)


    num_regions <- dplyr::n_distinct(plan_matrix[,1])

    unnormalized_log_density <- compute_plans_log_optimal_weights(
        adj_list, counties, pop,
        constraints, pop_temper, rho=compactness,
        splitting_schedule,
        ndists, num_regions=num_regions,
        lower=pop_bounds[1],
        target=pop_bounds[2],
        upper=pop_bounds[3],
        plan_matrix-1L, sizes_matrix,
        num_threads
    )

    # check if anything was -Inf so probability 0
    num_prob_zero <- sum(!is.finite(unnormalized_log_density))

    if(num_prob_zero > 0){
        cli::cli_warn("{num_prob_zero} of the plans have probability 0!")
    }

    return(unnormalized_log_density)
}
