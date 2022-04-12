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
#' @param alpha A numeric value for the alpha parameter for the talisman metric
#' @param beta A numeric value for the beta parameter for the talisman metric
#' @return Numeric vector with competitiveness scores
#' @export
#'
#' @concept analyze
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' comp <- redist.competitiveness(plans_05, fl25$mccain, fl25$obama)
#'
redist.competitiveness <- function(plans, rvote, dvote, alpha = 1, beta = 1){
    nd <- length(unique(plans[, 1]))
    redistmetrics::compet_talisman(plans = plans, shp = data.frame(),
                                   rvote = rvote, dvote = dvote,
                                   alpha = alpha, beta = beta)[seq(1, nd * ncol(plans), by = nd)]
}
