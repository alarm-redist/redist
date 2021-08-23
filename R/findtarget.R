#' Find Majority Minority Remainder
#'
#' Given a percent goal for majority minority districts, this computes the average
#' value of minority in non-majority minority districts. This value is "tgt_other"
#' in \code{redist.flip} and \code{redist_smc}.
#'
#' @param tgt_min target group population for majority minority district
#' @param group_pop A vector of populations for some subgroup of interest.
#' @param total_pop A vector containing the populations of each geographic unit.
#' @param ndists The number of congressional districts.
#' @param nmmd The number of majority minority districts.
#'
#' @return numeric value to target
#'
#' @concept prepare
#' @export
redist.find.target <- function(tgt_min, group_pop, total_pop, ndists, nmmd){
  totpop <- sum(total_pop)
  targetpop <- totpop/ndists
  tmm <- nmmd*tgt_min*targetpop
  totgroup <- sum(group_pop)
  tgt_other <- (totgroup - tmm)/((ndists-nmmd)*targetpop)
  #(sum(group_pop) - nmmd*tgt_min*targetpop)/((ndists-nmmd)*sum(total_pop)/ndists)
  return(c(tgt_other = tgt_other))
}

#' Create Constraints for SMC
#'
#' @param constraints Vector of constraints to include. Currently only 'vra' implemented.
#' @param tgt_min Defaults to 0.55. If 'vra' included, the minority percent to encourage in each district.
#' @param group_pop A vector of populations for some subgroup of interest.
#' @param total_pop A vector containing the populations of each geographic unit.
#' @param ndists The total number of districts.
#' @param nmmd The number of majority minority districts to target for 'vra' constraint
#' @param strength_vra The strength of the 'vra' constraint. Defaults to 2500.
#' @param pow_vra  The exponent for the 'vra' constraint. Defaults to 1.5.
#'
#' @return list of lists for each constraint selected
#'
#' @concept prepare
#' @export
redist.constraint.helper <- function(constraints = 'vra', tgt_min = 0.55,
                                     group_pop, total_pop, ndists, nmmd,
                                     strength_vra = 2500, pow_vra = 1.5){


  ret <- list()

  if('vra' %in% constraints){
    tgt_other <- redist.find.target(tgt_min, group_pop, total_pop, ndists, nmmd)

    ret['vra'] <- list(strength = strength_vra,
                       min_pop = group_pop,
                       tgt_vra_min = tgt_min,
                       tgt_vra_other = tgt_other,
                       pow_vra = 1.5)

  }
  return(ret)
}

globalVariables(c('vra'))
