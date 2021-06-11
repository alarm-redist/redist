#' Calculate Group Percent by District
#'
#' \code{redist.group.percent} computes the percentage that a group makes up in
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
#'                     group_pop = fl25$BlackPop,
#'                     total_pop = fl25$TotPop)
#'
redist.group.percent <- function(plans, group_pop, total_pop, ncores = 1) {

    if (!any(class(total_pop) %in% c('numeric', 'integer')))
        stop('Please provide "total_pop" as a numeric vector.')
    if (!any(class(group_pop) %in% c('numeric', 'integer')))
        stop('Please provide "groupop" as a numeric vector.')

    if (!any(class(plans) %in% c('integer', 'numeric', 'matrix')))
        stop('Please provide "plans" as a matrix.')

    if (!is.matrix(plans)) {
        plans <- as.matrix(plans)
    }

    if (length(total_pop) != nrow(plans))
        stop('Arguments "plans" and "total_pop" do not have same number of precincts.')
    if (length(group_pop) != nrow(plans))
        stop('Arguments "plans" and "group_pop" do not have same number of precincts.')

    ndists <- max(plans[, 1])
    if(ndists ==  length(unique(plans[,1])) - 1 ){
        plans <- plans + 1
        ndists <- ndists + 1
    }
    group_pct(plans, group_pop, total_pop, ndists)
}
