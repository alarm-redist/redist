#' Compute Competitiveness
#'
#' Currently only implements the competitiveness function in equation (5)
#' of Cho & Liu 2016.
#'
#'
#' @param plans A numeric vector (if only one map) or matrix with one row
#' for each precinct and one column for each map. Required.
#' @param rvote A numeric vector with the Republican vote for each precinct.
#' @param dvote A numeric vector with the Democratic vote for each precinct.
#' @return Numeric vector with competitiveness scores
#' @export
#'
#' @concept analyze
#' @examples \dontrun{
#' data(fl25)
#' data(fl25_enum)
#'
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' comp <- redist.competitiveness(plans_05, fl25$mccain, fl25$obama)
#' }
redist.competitiveness <- function(plans, rvote, dvote){
  # dont just copy this
  dvs <- redist.metrics(plans = plans, "DVS", rvote, dvote)
  dvs <- matrix(data = c(dvs$DVS), ncol = ncol(plans))
  # when pasting %TODO% -- just use dvs as is or will break
  nd <- length(unique(plans[,1]))
  talisman(dvs = dvs, nd = nd)
}
