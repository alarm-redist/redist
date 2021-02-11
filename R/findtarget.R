#' Find Majority Minority Remainder
#' 
#' Given a percent goal for majority minority districts, this computes the average 
#' value of minority in non-majority minority districts. This value is "tgt_other" 
#' in \code{redist.mcmc} and \code{redist.smc}.
#'
#' @param tgt_min target group population for majority minority district
#' @param grouppop A vector of populations for some subgroup of interest.
#' @param fullpop A vector containing the populations of each geographic unit.
#' @param ndists The number of congressional districts.
#' @param nmmd The number of majority minority districts.
#'
#' @return numeric value to target
#' @export
redist.find.target <- function(tgt_min, grouppop, fullpop, ndists, nmmd){
  totpop <- sum(fullpop)
  targetpop <- totpop/ndists
  tmm <- nmmd*tgt_min*targetpop
  totgroup <- sum(grouppop)
  tgt_other <- (totgroup - tmm)/((ndists-nmmd)*targetpop)
  #(sum(grouppop) - nmmd*tgt_min*targetpop)/((ndists-nmmd)*sum(fullpop)/ndists)
  return(c(tgt_other = tgt_other))
}

#' Create Constraints for SMC
#'
#' @param constraints Vector of constraints to include. Currently only 'vra' implemented.
#' @param tgt_min Defaults to 0.55. If 'vra' included, the minority percent to encourage in each district.
#' @param grouppop The minority group population by precinct, used by 'vra' constraint.
#' @param fullpop The full population by precinct.
#' @param ndists The total number of districts.
#' @param nmmd The number of majority minority districts to target for 'vra' constraint
#' @param strength_vra The strength of the 'vra' constraint. Defaults to 2500.
#' @param pow_vra  The exponent for the 'vra' constraint. Defaults to 1.5.
#'
#' @return list of lists for each constraint selected
#' @export
redist.constraint.helper <- function(constraints = 'vra', tgt_min = 0.55, grouppop, fullpop, ndists, nmmd, strength_vra = 2500, pow_vra = 1.5){
  ret <- list()
  
  if(vra %in% constraints){
    tgt_other <- redist.find.target(tgt_min, grouppop, fullpop, ndists, nmmd)
    
    ret['vra'] <- list(strength = strength, 
                       min_pop = grouppop, 
                       tgt_vra_min = tgt_min, 
                       tgt_vra_other = tgt_other,
                       pow_vra = 1.5)
    
  }
  return(ret)
}

globalVariables(c('vra'))