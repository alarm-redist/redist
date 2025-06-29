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
#' @inheritParams redist_smc
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
        map, plans, counties = NULL,
        sizes_matrix = NULL, compactness = 1L,
        constraints = list(),
        num_threads = 0, pop_temper = 0L
){

    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }
    num_regions <- dplyr::n_distinct(plan_matrix[,1])


    # get validated inputs
    map_params <- get_map_parameters(map, !!rlang::enquo(counties))
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


    # need to pass in the defused quosure
    constraints <- validate_constraints(map=map, constraints=rlang::enquo(constraints))

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
#' @inheritParams redist_smc
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
        map, plans, counties = NULL,
        sizes_matrix = NULL, compactness = 1L,
        constraints = list(),
        num_threads = 0, pop_temper = 0L
){
    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }
    num_regions <- dplyr::n_distinct(plan_matrix[,1])


    # get validated inputs
    map_params <- get_map_parameters(map, !!rlang::enquo(counties))
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


    # need to pass in the defused quosure
    constraints <- validate_constraints(map=map, constraints=rlang::enquo(constraints))

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


#' Computes the unnormalized log incremental weights for plan
#'
#'
#' @inheritParams redist_smc
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
        map, plans, counties = NULL,
        sizes_matrix = NULL, constraints = list(),
        splitting_schedule = "any_valid_sizes",
        num_threads = 0, compactness = 1L, pop_temper = 0L
){
    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }
    num_regions <- dplyr::n_distinct(plan_matrix[,1])


    # get validated inputs
    map_params <- get_map_parameters(map, !!rlang::enquo(counties))
    map <- map_params$map
    V <- map_params$V
    adj_list <- map_params$adj_list
    counties <- map_params$counties
    pop <- map_params$pop
    pop_bounds <- map_params$pop_bounds

    ndists <- attr(map, "ndists")
    total_seats <- attr(map, "total_seats")
    district_seat_sizes <- attr(map, "district_seat_sizes")
    districting_scheme <- attr(map, "districting_scheme")
    if(districting_scheme == "MMD"){
        splitting_schedule <- "split_district_only_mmd"
    }
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


    # need to pass in the defused quosure
    constraints <- validate_constraints(map=map, constraints=rlang::enquo(constraints))

    unnormalized_log_density <- compute_plans_log_optimal_weights(
        adj_list, counties, pop,
        constraints, pop_temper, rho=compactness,
        splitting_schedule,
        ndists, total_seats, num_regions=num_regions,
        district_seat_sizes,
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
