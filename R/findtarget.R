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

