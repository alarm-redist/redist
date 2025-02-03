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



#' Takes a map and possibly a county label and returns map parameters for SMC algorithms
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
OLD_validate_cut_sizes <- function(ndists, num_splits,
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


# Starts a list
# check that for each split the allowable sizes make sense
# also check that sizes make sense
validate_cut_sizes <- function(
        ndists, num_splits, initial_num_regions,
        permitted_split_region_sizes_list){
    if(initial_num_regions > 1){
        cli_abort("Not supported yet!")
    }
    # check actually a list
    if(!is.list(permitted_split_region_sizes_list)){
        cli_abort(
            "{.arg permitted_split_region_sizes_list} must be a list of numeric vectors!"
        )
    }
    # check each element of the list is a numeric vector
    if(!all(sapply(permitted_split_region_sizes_list, is.numeric))){
        cli_abort(
            "Each element of {.arg permitted_split_region_sizes_list} must be a numeric vector!"
        )
    }
    # check each element is an integer
    if(!all(sapply(permitted_split_region_sizes_list, function(x) all.equal(x, as.integer(x))))){
        cli_abort(
            "Each element of each entry in {.arg permitted_split_region_sizes_list} must be an integer!"
        )
    }
    # check for each element the max is less than ndists - num_splits
    for (split_num_i in seq_len(num_splits)) {
        # the number of regions after the split
        n_regions_after_split <- initial_num_regions + split_num_i
        max_possible_region_size <- ndists - n_regions_after_split + 1
        # check that all sizes are between 1 (inclusive) and n_regions_after_split -1
        if(any(permitted_split_region_sizes_list[[split_num_i]] < 1)){
            cli_abort(
                "An element of entry {split_num_i} of {.arg permitted_split_region_sizes_list} was less than 1! All valid region sizes must be at least 1!"
            )
        }
        # check that everything is less than or eqaul to n_regions_after_split -1
        if(any(permitted_split_region_sizes_list[[split_num_i]] > max_possible_region_size)){
            cli_abort(
                "An element of entry {split_num_i} of {.arg permitted_split_region_sizes_list} was greater than {max_possible_region_size}! All valid region sizes at step {split_num_i} must be less than {max_possible_region_size}!"
            )
        }
    }

    # for now just return. No other validation performed
    return()


    # Now check that the sizes passed for each step are all
    # it is possible to obtain them from the set of valid sizes of the step before
    # there are no unaccounted for sizes
    prev_possible_sizes <- c(ndists)
    for (split_i in seq_len(num_splits)) {
        # vector of previous step region sizes which can be split in this step
        # we want to make sure by the end this has at least one element
        prev_possible_sizes_that_can_be_split <- vector(mode = "numeric", length=0L)
        # the user input of permitted sizes
        split_i_permitted_region_sizes <- permitted_split_region_sizes_list[[split_i]]
        # vector used to track the actual sizes permitted by this input
        # and the previous round
        split_i_possible_region_sizes <- vector(mode = "numeric", length=0L)
        # check for each value in the permitted size that we could get there
        for (a_region_size in split_i_permitted_region_sizes) {
            # for each region
            # 1. find all regions in previous step bigger than that
            # 2. check that the other pieces made by splitting are in here

            # Get the potential sizes made if we split `a_region_size` from each of the
            # previous possible sizes
            possible_other_split_sizes <- prev_possible_sizes - a_region_size
            # boolean for each previous possible size if it can be split
            # which means the other split size is a valid size
            split_possible_mask <- possible_other_split_sizes %in% split_i_permitted_region_sizes

            # allow for this region to left unsplit
            # NOTE: This might be missing some edge cases but idk rn
            split_i_possible_region_sizes <- append(
                split_i_possible_region_sizes,
                a_region_size
            )

            # if all false no split possible and this region will be left untouched
            # but if at least one is true it can be split
            if(any(split_possible_mask)){
                # else at least one split is possible so record that
                split_i_possible_region_sizes <- append(
                    split_i_possible_region_sizes,
                    possible_other_split_sizes[split_possible_mask]
                )
                # add those to vector of previous regions that we can split
                prev_possible_sizes_that_can_be_split <- append(
                    prev_possible_sizes_that_can_be_split,
                    prev_possible_sizes[split_possible_mask]
                )
            }
        }

        # now make them unique
        prev_possible_sizes_that_can_be_split <- unique(prev_possible_sizes_that_can_be_split)
        split_i_possible_region_sizes <- unique(split_i_possible_region_sizes)

        # Now check the following
        # 1. at least one of the previous regions can be split
        if(length(prev_possible_sizes_that_can_be_split) == 0){
            cli_abort(
            "It is not possible to split any of the valid sizes from before split {split_i} into a valid size for split {split_i}"
            )
        }
        # 2. There are no implied sizes that were not entered
        # e.g. ndists=30 and only allow size 9 (when 21 is also implied)
        if(!base::setequal(split_i_possible_region_sizes, split_i_permitted_region_sizes)){
            cli_abort(
            "Inputted permissable sizes for split {split_i} are {sort(split_i_permitted_region_sizes)} but given the previous sizes this would actually produce sizes {sort(split_i_possible_region_sizes)} which is not the same as the input!"
            )
        }
        # now set prev possible sizes to the ones from this round and repeat
        prev_possible_sizes <- split_i_possible_region_sizes
    }



    # need to check
}
