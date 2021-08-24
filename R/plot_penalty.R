#' Visualize VRA Penalty
#'
#' Plots the shape of the VRA Gibbs penalty.
#'
#' This function allows you to plot the un-exponentiated VRA Gibbs penalty implemented as
#' `vra` within MCMC, Merge-Split, and SMC. The function takes two key inputs,
#' `tgt_min` and `tgt_other` which center the minimum penalty spots. A higher y-value
#' indicates a higher penalty and incentivizes moving towards a spot with a lower y-value.
#' The x-axis indicates the group population proportion in a given district. The default
#' indicates the default settings within \code{redist_smc()}.
#'
#' @param tgt_min double, defaults to 0.55. The minority target percent.
#' @param tgt_other double, defaults to 0.25. The other group target percent.
#' @param strength_vra double, strength of the VRA constraint.
#' @param pow_vra double, exponent of the VRA constraint.
#' @param limits Whether to limit y axis to 0,500. Default is TRUE for
#' comparability across values.
#'
#' @return ggplot
#'
#' @concept prepare
#' @export
#'
#' @importFrom ggplot2 lims labs geom_path
#'
#' @examples
#' redist.plot.penalty()
#'
redist.plot.penalty <- function(tgt_min = 0.55, tgt_other = 0.25,
                                     strength_vra = 2500, pow_vra = 1.5, limits = TRUE){

  tb <- tibble(x = seq(0, 1, by = .001))
  tb <- tb %>% mutate(y = strength_vra*(abs(x - tgt_min)^pow_vra)*(abs(x-tgt_other)^pow_vra))


  if(limits) {
    ret <- tb %>% filter(y < 500) %>%
      ggplot(aes(x = x, y = y)) +
      geom_path() +
      theme_bw() +
      lims(x = c(0,1), y = c(0,500)) +
      labs(x = 'Group Population', y = 'Penalty Size')
  } else {
    ret <- tb %>%
      ggplot(aes(x = x, y = y)) +
      geom_path() +
      theme_bw() +
      labs(x = 'Group Population', y = 'Penalty Size')
  }


  return(ret)
}

globalVariables(c('x', 'y'))
