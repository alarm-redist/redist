#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2025/4/2
# Purpose: Functions for computing unnormalized target density of plans
# and optimal weight functions
####################################################


#' Computes the unnormalized log target density for plans
#'
#' Computes the unnormalized log target density for the standard redist
#' target measure. The density for a plan  can be decomposed into a plan-term
#' (if relevant) and the procuct over regions of the regions compactness and
#' its score
#'
#' @inheritParams redist_gsmc
#' @param plans Either a `redist_plans` object or a matrix or vector of plan labels
#' @param num_threads The number of threads used for computing
#' the target densities.
#' @param sizes_matrix A matrix of region sizes. The matrix must have
#' the same number of columns as `plans`, `ndists` rows, and each column
#' must sum to `ndists`.
#'
#' @return A vector of the unnormalized target density of each plan
#'
#' @md
#' @order 1
#' @export
compute_log_target_density <- function(
        map, plans, sizes_matrix = NULL,
        counties = NULL, compactness = 1L,
        constraints = list(),
        num_threads = 0, pop_temper = 0L
        ){

    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(enum_plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }
    num_regions <- dplyr::n_distinct(plan_matrix[,1])

    counties_q <- rlang::enquo(counties)
    # get validated inputs
    map_params <- get_map_parameters(map, counties_q=counties_q, use_counties_q = TRUE)
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    ndists <- attr(map, "ndists")
    total_seats <- attr(map, "total_seats")
    district_seat_sizes <- attr(map, "district_seat_sizes")
    storage.mode(district_seat_sizes) <- "integer"

    if (inherits(plans, "redist_plans")){
        sizes_matrix <- get_nseats_matrix(plans)
    }else if(is.null(sizes_matrix)){
        # infer
        prec_pop <- map[[attr(map, "pop_col")]]
        distr_pop <- pop_tally(plan_matrix, prec_pop, num_regions)
        sizes_matrix <- infer_region_sizes(
            distr_pop,
            attr(map, "pop_bounds")[1], attr(map, "pop_bounds")[3],
            total_seats
        )
    }


    # need to pass in quosure
    constraints_q <- rlang::enquo(constraints)
    constraints <- validate_constraints(map, constraints_q=constraints_q, use_constraints_q=TRUE)

    unnormalized_log_density <- compute_log_unnormalized_plan_target_density(
        adj_list, counties, pop,
        constraints, pop_temper, rho=compactness,
        ndists=ndists, total_seats = total_seats, num_regions=num_regions,
        district_seat_sizes=district_seat_sizes,
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



#' Computes the unnormalized log target density contribution for each region in a plan
#'
#' Computes the unnormalized log target density contribution for each region in
#' a plan. This amounts to the log number of spanning trees minus the score. Note
#' that if the target measure has a plan-wide component then the sum of the log
#' densities over regions is not neccesarily equal to the log density of an entire plan.
#'
#' @inheritParams redist_gsmc
#' @param num_threads The number of threads used for computing
#' the target densities.
#' @param plan_matrix A matrix of 1-indexed plans
#' @param init_region_sizes_mat A matrix of region sizes. The matrix must have
#' the same number of columns as `plan_matrix`, `ndists` rows, and each column
#' must sum to `ndists`.
#'
#' @return A matrix of number of regions by plans of the unnormalized target density of each region in a plan
#'
#' @md
#' @order 1
#' @export
compute_log_target_density_by_region <- function(
        map, plans, sizes_matrix = NULL,
        counties = NULL, compactness = 1L,
        constraints = list(),
        num_threads = 0, pop_temper = 0L
){
    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(enum_plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }
    num_regions <- dplyr::n_distinct(plan_matrix[,1])


    counties_q <- rlang::enquo(counties)
    # get validated inputs
    map_params <- get_map_parameters(map, counties_q=counties_q, use_counties_q = TRUE)
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    ndists <- attr(map, "ndists")
    total_seats <- attr(map, "total_seats")
    district_seat_sizes <- attr(map, "district_seat_sizes")
    storage.mode(district_seat_sizes) <- "integer"

    if (inherits(plans, "redist_plans")){
        sizes_matrix <- get_nseats_matrix(plans)
    }else if(is.null(sizes_matrix)){
        # infer
        prec_pop <- map[[attr(map, "pop_col")]]
        distr_pop <- pop_tally(plan_matrix, prec_pop, num_regions)
        sizes_matrix <- infer_region_sizes(
            distr_pop,
            attr(map, "pop_bounds")[1], attr(map, "pop_bounds")[3],
            total_seats
        )
    }


    # need to pass in quosure
    constraints_q <- rlang::enquo(constraints)
    constraints <- validate_constraints(map, constraints_q=constraints_q, use_constraints_q=TRUE)

    unnormalized_log_region_densities <- compute_log_unnormalized_region_target_density(
        adj_list, counties, pop,
        constraints, pop_temper, rho=compactness,
        ndists=ndists, total_seats = total_seats, num_regions=num_regions,
        district_seat_sizes=district_seat_sizes,
        lower=pop_bounds[1],
        target=pop_bounds[2],
        upper=pop_bounds[3],
        region_ids=plan_matrix-1L,
        region_sizes=sizes_matrix,
        num_threads=num_threads
    )

    # check if anything was -Inf so probability 0
    num_prob_zero <- sum(!is.finite(unnormalized_log_region_densities))

    if(num_prob_zero > 0){
        cli::cli_warn("{num_prob_zero} of the {prod(dim(unnormalized_log_region_densities))} regions have probability 0!")
    }

    return(unnormalized_log_region_densities)
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
#' @inheritParams redist_gsmc
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
    }else if(is.vector(sizes_matrix)){
        sizes_matrix <- as.matrix(sizes_matrix, cols = 1)
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
