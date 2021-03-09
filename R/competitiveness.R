#' Compute Competitiveness
#'
#' Currently only implements the competitiveness function in equation (5)
#' of Cho & Liu 2016.
#'
#'
#' @param district_membership A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param rvote A numeric vector with the Republican vote for each precinct.
#' @param dvote A numeric vector with the Democratic vote for each precinct.
#' @return Numeric vector with competitiveness scores
#' @export
#'
#' @concept analyze
#' @examples \dontrun{
#' data("algdat.p10")
#' comp <- redist.competitiveness(algdat.p10$cdmat,
#' algdat.p10$precinct.data$repvote,
#' algdat.p10$precinct.data$demvote)
#' }
redist.competitiveness <- function(district_membership, rvote, dvote){
  # dont just copy this
  dvs <- redist.metrics(district_membership = district_membership, "DVS", rvote, dvote)
  dvs <- matrix(data = c(dvs$DVS), ncol = ncol(district_membership))
  # when pasting %TODO% -- just use dvs as is or will break
  nd <- length(unique(district_membership[,1]))
  talisman(dvs = dvs, nd = nd)
}
