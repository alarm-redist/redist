#' Calculate Group Percent by District
#'
#' \code{redist.group.percent} computes the percentage that a group makes up in
#' each district across a matrix of maps.
#'
#' @param district_membership A matrix with one row
#' for each precinct and one column for each map. Required.
#' @param group_pop A numeric vector with the population of the group for every precinct.
#' @param full_pop A numeric vector with the population for every precinct.
#'
#' @return A matrix with a row for each district and a column for each plan,
#'  with the entries representing the percentage of the group in that district,
#'  relative to the total population of the district.
#'
#' @export
#'
#' @examples \dontrun{
#' data("fl25")
#' data(algdat)
#' cd <- algdat.p10$cdmat[,1:5]
#'
#' redist.group.percent(district_membership = cd,
#'                     group_pop = fl25$BlackPop,
#'                     full_pop = fl25$TotPop)
#' }
redist.group.percent <- function(district_membership, group_pop, full_pop){
    if (!any(class(full_pop) %in% c('numeric', 'integer')))
        stop('Please provide "full_pop" as a numeric vector.')
    if (!any(class(group_pop) %in% c('numeric', 'integer')))
        stop('Please provide "group_pop" as a numeric vector.')

    if (!any(class(district_membership) %in% c('integer', 'numeric', 'matrix')))
        stop('Please provide "district_membership" as a matrix.')
    if (!is.matrix(district_membership))
        district_membership <- as.matrix(district_membership)
    if (min(district_membership[,1]) == 0)
        district_membership = 1 + district_membership

    if (length(full_pop) != nrow(district_membership))
        stop('Arguments "district_membership" and "fullpop" do not have same number of precincts.')
    if (length(group_pop) != nrow(district_membership))
        stop('Arguments "district_membership" and "groupop" do not have same number of precincts.')

    n_distr = max(district_membership[,1])
    group_pct(district_membership, group_pop, full_pop, n_distr)
}
