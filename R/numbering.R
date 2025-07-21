
# helper function for match_numbers
find_numbering <- function(plan, ref, pop, tot_pop) {
    joint <- plan_joint(ref, plan, pop)

    renumb <- solve_hungarian(1 - joint/tot_pop)[, 2]

    list(renumb = renumb,
        shared = sum(diag(joint[, renumb]))/tot_pop)
}

#' Renumber districts to match an existing plan
#'
#' District numbers in simulated plans are by and large random.  This
#' function attempts to renumber the districts across all simulated plans to
#' match the numbers in a provided plan, using the Hungarian algorithm.
#'
#' @param data a \code{redist_plans} object.
#' @param plan a character vector giving the name of the plan to match to (e.g.,
#' for a reference plan), or an integer vector containing the plan itself.
#' @param total_pop a vector of population counts. Should not be needed for most
#' \code{redist_plans} objects.
#' @param col the name of a new column to store the vector of population overlap
#' with the reference plan: the fraction of the total population who are in
#' the same district under each plan and the reference plan. Set to
#' \code{NULL} if no column should be created.
#' renumbering options in any plan.
#'
#' @returns a modified \code{redist_plans} object. New district numbers will be
#' stored as an ordered factor variable in the \code{district} column. The
#' district numbers in the plan matrix will match the levels of this factor.
#'
#' @examples
#' data(iowa)
#'
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)
#' plans <- redist_smc(iowa_map, 100, silent = TRUE)
#' match_numbers(plans, "cd_2010")
#'
#' @concept analyze
#' @export
match_numbers <- function(data, plan, total_pop = attr(data, "prec_pop"), col = "pop_overlap") {
    if (!inherits(data, "redist_plans")) cli::cli_abort("{.arg data} must be a {.cls redist_plans}")
    if (!"district" %in% names(data)) cli::cli_abort("Missing {.field district} column in {.arg data}")

    plan_mat <- get_plans_matrix(data)
    if (is.character(plan)) plan <- plan_mat[, plan]
    plan <- factor(plan, ordered = TRUE)
    ndists <- length(levels(plan))


    if (is.null(total_pop))
        cli::cli_abort("Must provide {.arg total_pop} for this {.cls redist_plans} object.")
    if (max(plan_mat[, 1]) != ndists)
        cli::cli_abort("Can't match numbers on a subset of a {.cls redist_plans}")
    if (length(plan) != nrow(plan_mat))
        cli::cli_abort(c("{.arg plan} doesn't have the right length.",
                    "i"="{.code length(plan)} should match the number of precincts,
                    i.e., {.code nrow(get_plans_matrix(data))}."))

    # compute renumbering and extract info
    best_renumb <- apply(plan_mat, 2, find_numbering,
        plan = as.integer(plan),
        pop = total_pop, tot_pop = sum(total_pop))
    renumb <- as.integer(vapply(best_renumb, function(x) x$renumb, integer(ndists)))

    if (!is.null(col))
        data[[col]] <- as.numeric(vapply(best_renumb, function(x) rep(x$shared, ndists),
            numeric(ndists)))

    renumb_mat <- renumber_matrix(plan_mat, renumb)
    colnames(renumb_mat) <- colnames(plan_mat)
    data <- set_plan_matrix(data, renumb_mat)
    data$district <- factor(levels(plan)[renumb], levels(plan), ordered = TRUE)

    orig_groups <- dplyr::group_vars(data)
    dplyr::group_by(data, .data$draw) %>%
        dplyr::arrange(.data$district, .by_group = TRUE) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(orig_groups)))
}

#' Renumber districts to match a quantity of interest
#'
#' District numbers in simulated plans are by and large random.  This
#' function will renumber the districts across all simulated plans in order
#' of a provided quantity of interest.
#'
#' @param data a \code{redist_plans} object
#' @param x \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the quantity of interest.
#' @param desc \code{TRUE} if district should be sorted in descending order.
#'
#' @returns a modified \code{redist_plans} object. New district numbers will be
#' stored as an ordered factor variable in the \code{district} column. The
#' district numbers in the plan matrix will match the levels of this factor.
#'
#' @concept analyze
#' @export
number_by <- function(data, x, desc = FALSE) {
    if (!inherits(data, "redist_plans")) cli::cli_abort("{.arg data} must be a {.cls redist_plans}")
    if (!"district" %in% names(data)) cli::cli_abort("Missing {.field district} column in {.arg data}")

    ord <- 1 - 2*desc
    m <- get_plans_matrix(data)
    orig_groups <- dplyr::group_vars(data)
    dplyr::group_by(data, .data$draw) %>%
        dplyr::mutate(district = rank(ord*{{ x }}, ties.method = "random")) %>%
        set_plan_matrix(`colnames<-`(renumber_matrix(m, .$district), colnames(m))) %>%
        dplyr::arrange(district, .by_group = TRUE) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(orig_groups)))
}

