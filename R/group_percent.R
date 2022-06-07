#' Calculate Group Proportion by District
#'
#' \code{redist.group.percent} computes the proportion that a group makes up in
#' each district across a matrix of maps.
#'
#' @param plans A matrix with one row
#' for each precinct and one column for each map. Required.
#' @param group_pop A numeric vector with the population of the group for every precinct.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#' @return matrix with percent for each district
#'
#' @export
#' @concept analyze
#'
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' cd <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#'
#' redist.group.percent(plans = cd,
#'     group_pop = fl25$BlackPop,
#'     total_pop = fl25$TotPop)
#'
redist.group.percent <- function(plans, group_pop, total_pop, ncores = 1) {
    if (!is.numeric(group_pop) || !is.numeric(total_pop))
        cli_abort("{.arg group_pop} and {.arg total_pop} must be numeric vectors.")
    if (!is.matrix(plans))
        cli_abort("{.arg plans} must be a matrix.")

    if (!is.matrix(plans)) {
        plans <- as.matrix(plans)
    }

    if (length(total_pop) != nrow(plans))
        cli_abort("{.arg plans} and {.total_pop} must have the same number of precincts.")
    if (length(group_pop) != nrow(plans))
        cli_abort("{.arg plans} and {.group_pop} must have the same number of precincts.")

    ndists <- max(plans[, 1])
    if (ndists ==  length(unique(plans[, 1])) - 1) {
        plans <- plans + 1
        ndists <- ndists + 1
    }
    group_pct(plans, group_pop, total_pop, ndists)
}
