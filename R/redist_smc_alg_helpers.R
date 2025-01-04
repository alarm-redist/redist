#####################################################
# Author: Philip O'Sullivan
# Institution: Harvard University
# Date Created: 2024/12/26
# Purpose: Helper functions shared across all smc algorithm types
####################################################




#' Checks a county label input is valid and if so returns it.
#'
#' Checks that a county label input is valid and if so returns it.
#'
#' @param map A `redist_map`
#' @param V The number of vertices in the plan graph.
#' @param counties The county label input
#'
#'
#' @returns A counties label vector
validate_counties <- function(map, adj_list, V, counties){

    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    if (is.null(counties)) {
        counties <- rep(1, V)
    } else {
        if (any(is.na(counties)))
            cli_abort("County vector must not contain missing values.")

        # handle discontinuous counties
        component <- contiguity(adj_list, vctrs::vec_group_id(counties))
        counties <- dplyr::if_else(component > 1,
                                   paste0(as.character(counties), "-", component),
                                   as.character(counties)) %>%
            as.factor() %>%
            as.integer()
        if (any(component > 1)) {
            cli_warn("Counties were not contiguous; expect additional splits.")
        }
    }

    return(counties)
}



#' Takes a map and possibly a county label and returns graph parameters for SMC algorithms
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
get_map_parameters <- function(map, counties=NULL){

    # get the map in adjacency form
    map <- validate_redist_map(map)
    V <- nrow(map)
    adj <- get_adj(map)


    counties <- validate_counties(map, adj, V, counties)


    # get population stuff
    pop_bounds <- attr(map, "pop_bounds")
    pop <- map[[attr(map, "pop_col")]]
    if (any(pop >= pop_bounds[3])) {
        too_big <- as.character(which(pop >= pop_bounds[3]))
        cli_abort(c("Unit{?s} {too_big} ha{?ve/s/ve}
                population larger than the district target.",
                    "x" = "Redistricting impossible."))
    }

    return(list(
        map=map,
        adj_list=adj,
        V=V,
        counties=counties,
        pop=pop,
        pop_bounds=pop_bounds
    ))

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
validate_initial_region_id_mat <- function(init_region_ids_mat, V, nsims, init_num_regions){
    # check that matrix dimension is V by nsims
    if(!is.matrix(init_region_ids_mat)){
        cli_abort("{.arg init_region_ids_mat} must be a matrix.")
    }
    if (nrow(init_region_ids_mat) != V)
        cli_abort("{.arg init_region_ids_mat} must have as many rows as {.arg map} has precincts.")
    if (ncol(init_region_ids_mat) != nsims)
        cli_abort("{.arg init_region_ids_mat} must have {.arg nsims} columns.")

    # check that every column only has values from 0 to num_regions-1
    if(any(colmin(init_region_ids_mat) != 0))
        cli_abort("{.arg init_region_ids_mat} can't have number less than 0.")
    if(any(colmax(init_region_ids_mat) != init_num_regions-1))
        cli_abort("{.arg init_region_ids_mat} can't have number greater than {.arg init_num_regions}-1.")

    expected_region_ids <- seq_len(init_num_regions)-1
    cols_as_expected <- apply(init_region_ids_mat,2, function(a_col) base::setequal(expected_region_ids, a_col))

    if(!any(cols_as_expected))
        cli_abort("{.arg init_region_ids_mat} can only have values between 0,...,{.arg init_num_regions}-1.")

}



#' Checks that a region size matrix is valid
#'
#' Checks that every size matrix is valid. A partial plan with `init_num_regions`
#' has a valid size matrix if it is `ndists` by `nsims` where
#'      - For each column rows 1:`init_num_regions` sum to 1 and each entry is
#'      between 1 and `ndists` inclusive
#'      - For each column rows `init_num_regions`+1:`ndists` are all zero
#'
#' Will throw an error if anything is wrong.
#'
#' @param init_dvals_mat A ndists by nsims matrix of region sizes
#' @param ndists The number of final districts in the plan.
#' @param nsims The number of simulations being run.
#' @param init_num_regions The number of regions in the partial plans stored in the
#' `init_region_ids_mat`
#'
validate_initial_region_sizes_mat <- function(init_dvals_mat, ndists, nsims, init_num_regions){
    # check that matrix dimension is ndists by nsims
    if(!is.matrix(init_dvals_mat)){
        cli_abort("{.arg init_dvals_mat} must be a matrix.")
    }
    if (nrow(init_dvals_mat) != ndists)
        cli_abort("{.arg init_dvals_mat} must have {.arg ndists} rows.")
    if (ncol(init_dvals_mat) != nsims)
        cli_abort("{.arg init_dvals_mat} must have {.arg nsims} columns.")

    # check each value is an integer
    if(any(init_dvals_mat != round(init_dvals_mat))){
        cli_abort("Each entry of {.arg init_dvals_mat} must sum an integer!")
    }

    # check that each column from 1:init_num_regions sums to ndists
    if(any(colSums(init_dvals_mat[1:init_num_regions, ,drop=FALSE]) != ndists))
        cli_abort("Each column of {.arg init_dvals_mat} must sum to {.arg ndists}.")
    # check that init_num_regions+1 to the end of the matrix is 0
    if(any(init_dvals_mat[(init_num_regions+1):ndists,] != 0))
        cli_abort("Rows {.arg init_num_regions}+1 to {.arg ndists} of each column of {.arg init_dvals_mat} must be 0.")
    # check that 1:init_num_regions is between 1 and ndists
    if(any(init_dvals_mat[1:init_num_regions,] < 1))
        cli_abort("All values in rows 1 to {.arg init_num_regions} of each column of {.arg init_dvals_mat} must be greater than 0.")
    if(any(init_dvals_mat[1:init_num_regions,] > ndists))
        cli_abort("All values in rows 1 to {.arg init_num_regions} of each column of {.arg init_dvals_mat} must be less or equal to {.arg ndists}.")

}



# TODO: Make this work with partial splits. Right now only works with
# starting from the beginning

#' Crude check that cut sizes are valid
#'
#' Crude check that the minimum and maximum cut sizes are valid
#'
#' @param ndists The number of districts in the final plan.
#' @param num_splits The number of splits to be performed.
#' @param min_region_cut_sizes Vector of minimum region cut sizes
#' for each split.
#' @param max_region_cut_sizes Vector of maximum region cut sizes
#' for each split.
#'
#' @returns A counties label vector
validate_cut_sizes <- function(ndists, num_splits,
        min_region_cut_sizes, max_region_cut_sizes
){
    # check the cut sizes are between 1 and num_regions at that step minus 1
    if(any(min_region_cut_sizes < 0)){
        cli_abort("Minimum region cut sizes must be 1 or greater!")
    }
    if(any(min_region_cut_sizes > rev(seq_len(ndists-1)))){
        cli_abort("Minimum region cut sizes must be less than {.arg ndists}!")
    }
    if(any(max_region_cut_sizes < 0)){
        cli_abort("Maximum region cut sizes must be 1 or greater!")
    }
    if(any(max_region_cut_sizes > rev(seq_len(ndists-1)))){
        cli_abort("Maximum region cut sizes must be less than {.arg ndists}!")
    }

    # check min cut size is always smaller
    if(any(min_region_cut_sizes > max_region_cut_sizes)){
        cli_abort("Maximum region cut sizes must be less than {.arg ndists}!")
    }


    # tests by trying the minimum split size each time
    potentential_sizes <- ndists
    for (i in 1:num_splits) {
        # check there's at least one region that can be split
        if(all(potentential_sizes <= min_region_cut_sizes[i])){
            cli_abort("For split {i} all region sizes are less than {min_region_cut_sizes[i]}")
        }

        # find the index of the biggest region by size
        biggest_piece_to_cut_index <- which.max(potentential_sizes)


        # get the biggest region size
        biggest_piece_to_cut <- potentential_sizes[biggest_piece_to_cut_index]

        bigger_cut_size <- max(
            biggest_piece_to_cut - min_region_cut_sizes[i],
            min_region_cut_sizes[i]
            )

        smaller_cut_size <- min(
            biggest_piece_to_cut - min_region_cut_sizes[i],
            min_region_cut_sizes[i]
        )

        untouched_region_indices <- setdiff(
            seq_len(length(potentential_sizes)),
            biggest_piece_to_cut_index
            )


        # split it into smallest cut size and whatever is left
        potentential_sizes <- c(
            potentential_sizes[untouched_region_indices],
            smaller_cut_size,
            bigger_cut_size
        )

    }


}
