#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/12/26
# Purpose: Helper functions shared across all redist algorithm types
####################################################




# simple but common arg checking functions
is_scalar <- function(x){
    return(length(x) == 1)
}



#' Takes a map and possibly a county label and returns map parameters for redist algorithms
#'
#' Takes a `redist_map` and possibly a county label vector and after checking
#' the inputs are valid it returns a list with the adjacency map, county vector,
#' and population variables
#'
#' @param map A `redist_map`
#' @param counties Either a numeric or factor vector mapping each vertex in the
#' graph to a county.
#'
#' @returns A list with 0-indexed adjancency list, county label vector,
#' population bounds vector, and population vector.
#' @noRd
get_map_parameters <- function(map, counties=NULL){

    # get the map in adjacency form
    map <- validate_redist_map(map)
    V <- nrow(map)
    adj_list <- get_adj(map)
    # get seat related values
    ndists <- attr(map, "ndists")
    nseats <- attr(map, "nseats")
    seats_range <- attr(map, "seats_range")
    storage.mode(seats_range) <- "integer"
    districting_scheme <- attr(map, "districting_scheme")


    counties <- rlang::eval_tidy(rlang::enquo(counties), map)

    # validate the counties
    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        # check if just the column name was passed in
        if(assertthat::is.scalar(counties) && assertthat::is.string(counties)){
            if(!counties %in% names(map)){
                cli::cli_abort("{counties} is not in the map!")
            }else{
                counties <- map[[counties]]
            }
        }

        if(any(is.na(counties))){
            cli::cli_abort("County vector must not contain missing values.")
        }

        # handle discontinuous counties
        component <- contiguity(adj_list, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
                                   paste0(as.character(counties), "-", component),
                                   as.character(counties)) %>%
            as.factor() %>%
            as.integer()
        if (any(component > 1)) {
            cli::cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }



    # get population stuff
    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    smallest_district_size <- attr(map, "seats_range") |> min()
    if (any(pop >= pop_bounds[3] * smallest_district_size)) {
        too_big <- as.character(which(pop >= pop_bounds[3] * smallest_district_size))
        cli::cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                population larger than the district target.",
                    "x" = "Redistricting impossible."))
    }

    # now check we don't have too many regions or counties
    max_possible_sizes <- maximum_input_sizes()
    if(dplyr::n_distinct(counties) > max_possible_sizes$max_counties ){
        cli::cli_abort("The maximum number of supported counties is {max_possible_sizes$max_counties}!")
    }
    if(V > max_possible_sizes$max_V ){
        cli::cli_abort("The maximum number of supported vertices in a map is {max_possible_sizes$max_V}!")
    }
    if(attr(map, "ndists") > max_possible_sizes$max_districts){
        cli::cli_abort("The maximum number of supported districts in a map is {max_possible_sizes$max_districts}!")
    }


    return(list(
        map=map,
        adj_list=adj_list,
        V=V,
        ndists=ndists,
        nseats=nseats,
        seats_range=seats_range,
        districting_scheme=districting_scheme,
        counties=counties,
        pop=pop,
        pop_bounds=pop_bounds
    ))

}


#' Checks that a constraint input is valid and returns it as a list ready for c++
#'
#' Checks that the sampling space/splitting method combination given is valid
#' and that the splitting parameters are valid.
#'
#' @inheritParams redist_smc
#'
#'
#' @returns A list of constraints which is safe to pass into c++ code
#' @noRd
validate_constraints <- function(map, constraints, compactness = 1L){
    # Other constraints
    if (!inherits(constraints, "redist_constr")) {
        constraints <- new_redist_constr(rlang::eval_tidy(rlang::enquo(constraints), map))
    }
    if (any(c("edges_removed", "log_st") %in% names(constraints))) {
        cli::cli_warn(c("{.var edges_removed} or {.var log_st} constraint found in
           {.arg constraints} and will be ignored.",
                   ">" = "Adjust using {.arg compactness} instead."))
    }
    if (any(c("poslby", "fry_hold") %in% names(constraints)) && compactness == 1) {
        cli::cli_warn("{.var polsby} or {.var fry_hold} constraint found in {.arg constraints}
                 with {.arg compactness == 1). This may disrupt efficient sampling.")
    }
    constraints <- as.list(constraints) # drop data attribute
    return(constraints)
}


#' Checks the sampling space and splitting method are valid
#'
#' Checks that the sampling space/splitting method combination given is valid
#' and that the splitting parameters are valid.
#'
#' @inheritParams redist_smc
#'
#'
#' @returns A list of splitting parameters which is safe to pass into c++ code
#' @noRd
validate_sample_space_and_splitting_method <- function(
        sampling_space, splitting_method, splitting_params, num_splitting_steps
    ){

    # first check its a valid sampling space
    sampling_space <- rlang::arg_match(sampling_space, values = VALID_SAMPLING_SPACES)
    # check splitting method
    splitting_method <- rlang::arg_match(splitting_method, values = VALID_FORWARD_KERNEL_TYPES)

    # create the
    cleaned_splitting_params <- list()

    # next check if its graph space then must be naive k
    if(sampling_space == GRAPH_PLAN_SPACE_SAMPLING){
        if(splitting_method != NAIVE_K_SPLITTING){
            cli::cli_abort("{.arg splitting_method} must be {NAIVE_K_SPLITTING} when sampling on Graph Plan Space")
        }
        if(!"estimate_cut_k" %in% names(splitting_params)){
            # If no adapt k mentioned default is to estimate
            splitting_params$estimate_cut_k <- T
            #cli::cli_abort("For Naive Top K splitting method {.arg estimate_cut_k} must be included in {.arg splitting_params}")
        }
        # check its a boolean
        if(!assertthat::is.flag(splitting_params$estimate_cut_k)){
            cli::cli_abort("{.arg estimate_cut_k} must be a Boolean!")
        }
        cleaned_splitting_params$estimate_cut_k <- splitting_params$estimate_cut_k

        # now check if estimating cut k then adapt_k_thresh must be in there
        if(splitting_params$estimate_cut_k){
            # now check its included
            if(!"adapt_k_thresh" %in% names(splitting_params)){
                # cli::cli_abort("If estimating k for Naive Top K splitting method then {.arg estimate_cut_k} must be included in {.arg splitting_params}")
                cli::cli_warn("No adaptive threshold set for Naive Top K estimation. Defaulting to {.arg .99}")
                # default is .99
                splitting_params$adapt_k_thresh <- .99
            }
            # check its a scalar
            if (!assertthat::is.number(splitting_params$adapt_k_thresh)){
                cli::cli_abort("{.arg adapt_k_thresh} must be a number")
            }
            # now check its between 0 and 1
            if (splitting_params$adapt_k_thresh < 0 | splitting_params$adapt_k_thresh > 1){
                cli::cli_abort("{.arg adapt_k_thresh} must lie in [0, 1].")
            }
            cleaned_splitting_params$adapt_k_thresh <- splitting_params$adapt_k_thresh
        }else{ # else check manual k parameter were passed in
            if(!"manual_k_params" %in% names(splitting_params)){
                cli::cli_abort("If not estimating k for Naive Top K splitting method then {.arg manual_k_params} must be included in {.arg splitting_params}")
            }
            # check all inputs are numbers
            if(!all(sapply(splitting_params$manual_k_params, assertthat::is.number))) {
                cli::cli_abort("Manual splitting k parameter values all be numbers!")
            }
            # check k param input if not estimating
            if(any(splitting_params$manual_k_params < 1)) {
                cli::cli_abort("Manual splitting k parameter values must be all at least 1.")
            }
            # if just a single number then repeat it
            if(length(splitting_params$manual_k_params) == 1 && floor(splitting_params$manual_k_params) == splitting_params$manual_k_params){
                splitting_params$manual_k_params <- rep(splitting_params$manual_k_params, num_splitting_steps)
            }else if(length(splitting_params$manual_k_params) != num_splitting_steps){
                cli::cli_abort("K parameter input must be either 1 value or number of smc steps!")
            }else if(any(floor(splitting_params$manual_k_params) != splitting_params$manual_k_params)){
                # if either the length is not ndists-1 or its not all integers then throw
                # error
                cli::cli_abort("K parameter values must be all integers")
            }
            cleaned_splitting_params$manual_k_params <- splitting_params$manual_k_params
        }
    }else if(sampling_space == FOREST_SPACE_SAMPLING || sampling_space == LINKING_EDGE_SPACE_SAMPLING){
        # check splitting method is not naive k
        if(splitting_method == NAIVE_K_SPLITTING){
            cli::cli_abort("{.arg splitting_method} cannot be {NAIVE_K_SPLITTING} when sampling on Spanning Forest Space")
        }else if(!splitting_method %in% VALID_FOREST_SPLITTING_METHODS){
            cli::cli_abort("{.arg splitting_method} of {splitting_method} is not a valid splitting method. It must be one of[{VALID_FOREST_SPLITTING_METHODS}]")
        }
        if(splitting_method == EXP_BIGGER_ABS_DEV_SPLITTING){
            if(!"splitting_alpha" %in% names(splitting_params)){
                cli::cli_abort("For {.arg splitting_method} of {splitting_method} then {.arg splitting_alpha} must be included in {.arg splitting_params}")
            }
            cleaned_splitting_params$splitting_alpha <- splitting_params$splitting_alpha
        }
    }

    return(cleaned_splitting_params)
}



#' Checks that a region id matrix is valid
#'
#' Checks that every plan in a matrix of partial plans is valid. A partial
#' plans matrix with `init_num_regions` is valid if it is a `V` by `nsims`
#' matrix where for each column it only has values between
#' 0,...,`init_num_regions`. Will throw an error if anything is wrong.
#'
#' @param init_region_ids_mat A V by nsims matrix of zero-indexed partial plans
#' @param V The number of vertices in the plan graph.
#' @param nsims The number of simulations being run.
#' @param @init_num_regions The number of regions in the partial plans stored in the
#' `init_region_ids_mat`
#'
#' @noRd
validate_initial_region_id_mat <- function(init_region_ids_mat, V, nsims, init_num_regions){
    # check that matrix dimension is V by nsims
    if(!is.matrix(init_region_ids_mat)){
        cli::cli_abort("{.arg init_region_ids_mat} must be a matrix.")
    }
    if (nrow(init_region_ids_mat) != V)
        cli::cli_abort("{.arg init_region_ids_mat} must have as many rows as {.arg map} has precincts.")
    if (ncol(init_region_ids_mat) != nsims)
        cli::cli_abort("{.arg init_region_ids_mat} must have {.arg nsims} columns.")

    # check that every column only has values from 0 to num_regions-1
    if(any(colmin(init_region_ids_mat) > 0))
        cli::cli_abort("{.arg init_region_ids_mat} must have at least one region with id 0.")
    if(any(colmin(init_region_ids_mat) < 0))
        cli::cli_abort("{.arg init_region_ids_mat} can't have number less than 0.")
    if(any(colmax(init_region_ids_mat) != init_num_regions-1))
        cli::cli_abort("{.arg init_region_ids_mat} can't have number greater than {.arg init_num_regions}-1.")

    expected_region_ids <- seq_len(init_num_regions)-1
    cols_as_expected <- apply(init_region_ids_mat,2, function(a_col) base::setequal(expected_region_ids, a_col))

    if(!any(cols_as_expected))
        cli::cli_abort("{.arg init_region_ids_mat} can only have values between 0,...,{.arg init_num_regions}-1.")

}



#' Checks that a region seats matrix is valid
#'
#' Checks that a region seats matrix is valid. For plans with `init_num_regions`
#' the seats matrix should:
#'      - Have dimensions `init_num_regions` by `nsims`
#'      - Each value must be at least `min(seats_range)` and no bigger than `nseats`
#'      - Each column must sum to `nseats`
#'
#' If any of these conditions are not true then the function will throw an error.
#' Else returns the validated matrix
#'
#'
#' @param init_seats A ndists by nsims matrix of region sizes
#' @param nseats The total number of seats in the plan
#' @param seats_range The number of seats a district is allowed to have.For
#' single member districting schemes this is always 1.
#' @param nsims The number of simulations being run.
#' @param init_num_regions The number of regions in partial plans stored in the
#' `init_seats`
#' @param split_districts_only Whether the plans should be all districts except
#' for possible the remainder (which should have the largest index).
#' @param ncores The number of threads to use when checking. Generally `1` is
#' fine but for a lot of plans more threads can speed things up.
#'
#' @noRd
validate_init_seats <- function(
        init_seats, nseats, seats_range,
        nsims, init_num_regions,
        split_districts_only,
        ncores = 1L){
    # check that its a matrix
    if(!is.matrix(init_seats)){
        cli::cli_abort("{.arg init_seats} must be a matrix.")
    }
    # check its integerish
    if(!rlang::is_integerish(init_seats)){
        cli::cli_abort("{.arg init_seats} must be all integers")
    }
    # update the storage mode if needed
    if(storage.mode(init_seats) != "integer"){
        storage.mode(init_seats) <- "integer"
    }
    # now check dimensions
    if (nrow(init_seats) != init_num_regions){
        cli::cli_abort("{.arg init_seats} must have {.arg init_num_regions} rows.")
    }
    if (ncol(init_seats) != nsims){
        cli::cli_abort("{.arg init_seats} must have {.arg nsims} columns.")
    }
    # now check values in c++
    validate_init_seats_cpp(
        init_seats, init_num_regions,
        nseats, seats_range,
        split_districts_only,
        ncores
    )

    # now return
    return(init_seats)
}




#' Attempts to infer the number of seats for each district in a plan
#'
#' Attempts to guess the number of seats assigned to each region in a plans
#' matrix.
#'
#' @param plans Either a [redist_plans] object or a 1-indexed plans matrix
#' @param nseats The total number of seats in the plan
#' @param precint_pops A vector of precinct population assignments with same
#' number of rows as the plans matrix.
#' @param lower_pop_bound The minimum population bounds for a single seat
#' @param upper_pop_bound The maximum population bounds for a single seat
#' @param num_threads The number of threads to use while calculating. Defaults to
#' maximum availible on the machine.
#'
#'
#' @returns A matrix of the number of seats for each region in `plans`
#' @noRd
infer_plan_seats <- function(
        plans, nseats, precint_pops,
        lower_pop_bound, upper_pop_bound,
        num_threads = 0L){

    if (inherits(plans, "redist_plans")){
        plan_matrix <- get_plans_matrix(plans)
    }else if(is.matrix(plans)){
        plan_matrix <- plans
    }else if(is.vector(plans) && is.numeric(plans)){
        plan_matrix <- as.matrix(plans, cols = 1)
    }else{
        cli::cli_abort("{.arg plans} must be a matrix or {.cls redist_plans} type!")
    }

    num_regions <- length(unique(plan_matrix[,1]))
    distr_pop <- pop_tally(plan_matrix, precint_pops, num_regions)
    sizes_matrix <- infer_region_seats(
        distr_pop,
        lower_pop_bound, upper_pop_bound,
        nseats, num_threads
    )
    return(sizes_matrix)

}
