#' Pull back plans to unmerged units
#'
#' Merging map units through \code{\link{merge_by}} or \code{\link{summarize}}
#' changes the indexing of each unit.  Use this function to take a set of
#' redistricting plans from a \code{redist} algorithm and re-index them to
#' be compatible with the original set of units.
#'
#' @param plans a \code{redist_plans} object
#' @param map optionally, a \code{redist_map} object, which will be used to set the new population vector
#'
#' @returns a new, re-indexed, \code{redist_plans} object
#'
#' @concept analyze
#' @export
pullback <- function(plans, map = NULL) {
    if (!inherits(plans, "redist_plans")) cli_abort("{.arg plans} must be a {.cls redist_plans}")

    merge_idx <- attr(plans, "merge_idx")
    if (is.null(merge_idx)) {
        cli_warn("No merged indexing found.")
        return(plans)
    }

    attr(plans, "merge_idx") <- NULL
    if (inherits(map, "redist_map")) {
        attr(plans, "prec_pop") <- map[[attr(map, "pop_col")]]
    } else {
        attr(plans, "prec_pop") <- NULL
    }

    set_plan_matrix(plans, get_plans_matrix(plans)[merge_idx, ])
}



#' Helper function to check types for tidy wrappers
#' @noRd
check_tidy_types <- function(map, .data) {
    if (!is.null(map) && !inherits(map, "data.frame"))
        cli_abort("{.arg map} must be a data frame")
    if (is.null(.data))
        cli_abort("Must provide {.arg .data} if not called within a {.pkg dplyr} verb")
    if (!inherits(.data, "redist_plans"))
        cli_abort("{.arg data} must be a {.cls redist_plans}")
}


#' Tally a variable by district
#'
#' @param map a `redist_map` object
#' @param x a variable to tally. Tidy-evaluated.
#' @param .data a `redist_plans` object or matrix of plans
#' @param num_threads The number of threads to use for the calculation. Can
#' significantly speedup calculations for a large number of plans. 
#'
#' @return a vector containing the tallied values by district and plan (column-major)
#'
#' @concept analyze
#' @export
tally_var <- function(map, x, .data = pl(), num_threads = 1) {
    check_tidy_types(map, .data)
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    x <- rlang::eval_tidy(rlang::enquo(x), map)
    as.numeric(pop_tally(get_plans_matrix(.data), x, attr(.data, "ndists"), num_threads))
}

#' @rdname redist.group.percent
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object or matrix of plans
#' @param num_threads The number of threads to use for the calculation. Can
#' significantly speedup calculations for a large number of plans.
#'
#' @concept analyze
#' @export
group_frac <- function(map, group_pop, total_pop = map[[attr(map, "pop_col")]],
                       .data = pl(), num_threads = 1) {
    check_tidy_types(map, .data)
    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    group_pop <- rlang::eval_tidy(rlang::enquo(group_pop), map)
    total_pop <- rlang::eval_tidy(rlang::enquo(total_pop), map)

    if (!is.numeric(group_pop) || !is.numeric(total_pop))
        cli_abort("{.arg group_pop} and {.arg total_pop} must be numeric vectors.")

    plans = get_plans_matrix(.data)
    if (length(total_pop) != nrow(plans))
        cli_abort("{.arg .data} and {.arg total_pop} must have the same number of precincts.")
    if (length(group_pop) != nrow(plans))
        cli_abort("{.arg .data} and {.arg group_pop} must have the same number of precincts.")

    as.numeric(group_pct(plans, group_pop, total_pop, attr(.data, "ndists"), num_threads))
}


#' Average a variable by precinct (Deprecated)
#'
#' Deprecated in favor of [proj_avg()].
#' Takes a column of a `redist_plans` object and averages it across a set of
#' `draws` for each precinct.
#'
#' @param plans a `redist_plans` object
#' @param x an expression to average. Tidy-evaluated in `plans`.
#' @param draws which draws to average. `NULL` will average all draws, including
#' reference plans. The special value `NA` will average all sampled draws. An
#' integer, logical, or character vector indicating specific draws may also be
#' provided.
#'
#' @return a vector of length matching the number of precincts, containing the average.
#'
#' @export
avg_by_prec <- function(plans, x, draws = NA) {
    .Deprecated("proj_avg")
    proj_avg(plans, {{ x }}, draws)
}


#' @rdname redist.parity
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param num_threads The number of threads to use for the calculation. Can
#' significantly speedup calculations for a large number of plans.
#' @param rep_by_district Whether or not to return a duplicated value for each 
#' district (so size `nsims*districts`) or return a single value for each plan.
#'
#' @concept analyze
#' @export
plan_parity <- function(map, .data = pl(), num_threads = 1, rep_by_district = TRUE) {
    check_tidy_types(map, .data)
    ndists <- attr(map, "ndists")
    total_pop <- map[[attr(map, "pop_col")]]
    if (is.null(total_pop)) cli_abort("Population vector missing from {.arg map}")

    if(rep_by_district){
        return(rep(max_dev(get_plans_matrix(.data), total_pop, ndists, num_threads),
                each = ndists))
    }else{
        return(max_dev(get_plans_matrix(.data), total_pop, ndists, num_threads))
    }
}


#' Identify which counties are split by a plan
#'
#' @param plan A vector of precinct/unit assignments
#' @param counties A vector of county names or county ids.
#'
#' @return A logical vector which is \code{TRUE} for precincts belonging to
#' counties which are split
#'
#' @concept analyze
#' @export
is_county_split <- function(plan, counties) {
    counties <- vctrs::vec_group_id(counties)
    as.logical((tapply(plan, counties, FUN = function(y) length(unique(y))) > 1)[counties])
}




#' Extract the last plan from a set of plans
#'
#' @param plans A \code{\link{redist_plans}} object
#'
#' @returns An integer vector containing the final plan assignment.
#'
#' @concept analyze
#' @export
last_plan <- function(plans) {
    plan_m <- get_plans_matrix(plans)
    plan_m[, ncol(plan_m)]
}

#' Extract the district assignments for a precinct across all simulated plans
#'
#' @param prec the precinct number
#' @param .data a \code{\link{redist_plans}} object
#'
#' @return integer vector, a row from a plans matrix
#'
#' @concept analyze
#' @export
prec_assignment <- function(prec, .data = pl()) {
    check_tidy_types(NULL, .data)

    m <- get_plans_matrix(.data)
    if (is.integer(prec)) {
        if (prec <= 0 || prec > nrow(m))
            cli_abort(c("{.arg prec} out of bounds",
                "i" = "There are {nrow(m)} precincts in these plans."))
    } else {
        cli_abort("{.arg prec} must be an integer index")
    }

    assignment <- m[prec, , drop = FALSE]
    if ("district" %in% names(.data) && is.factor(.data$district)) {
        lev <- levels(.data$district)
        assignment <- factor(lev[assignment], lev, ordered = is.ordered(.data$district))
    }

    assignment
}

#' Compute a matrix of precinct co-occurrences
#'
#' For a map with `n` precincts Returns an `n`-by-`n` matrix, where each
#' entry measures the fraction of the plans in which the row and column
#' precincts were in the same district.
#'
#' @param plans a [redist_plans] object.
#' @param which [`<data-masking>`][dplyr::dplyr_data_masking] which plans to
#' compute the co-occurrence over.  Defaults to all.
#' @param sampled_only if `TRUE`, do not include reference plans.
#' @param ncores the number of parallel cores to use in the computation.
#'
#' @return a symmetric matrix the size of the number of precincts.
#'
#' @concept analyze
#' @md
#' @export
prec_cooccurrence <- function(plans, which = NULL, sampled_only = TRUE, ncores = 1) {
    if (sampled_only)
        plans <- subset_sampled(plans)
    which <- eval_tidy(enquo(which), plans)
    plan_m <- get_plans_matrix(plans)
    if (is.null(which))
        which <- seq_len(ncol(plan_m))
    prec_cooccur(plan_m, which, ncores)
}
