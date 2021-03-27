#' Dangerous but Skinny Flip
#' 
#' Runs flip silently without checking input quality for use within other contexts
#' which already check things. It returns just a matrix of plans.
#'
#' @param adj zero indexed adjacency list
#' @param init_plan initial plan
#' @param total_pop total population
#' @param pop_tol maximum population deviance allowed
#' @param nsims number of steps to take
#' @param eprob edge cut probability
#' @param lambda number of components to swap
#' @param constraints constraint list
#'
#' @return matrix  with 1 indexed plans
#' 
#' @noRd
#' 
skinny_flips <- function(adj, init_plan, total_pop, pop_tol, nsims, eprob, lambda, constraints){
  
  
  algout <- swMH(aList = adj,
                 cdvec = init_plan,
                 cdorigvec = constraints$similarity$plan,
                 popvec = total_pop,
                 grouppopvec = constraints$group_pop,
                 areas_vec = constraints$compact$areas,
                 county_membership = constraints$counties,
                 borderlength_mat = constraints$compact$borderlength_mat,
                 nsims = nsims,
                 eprob = eprob,
                 pct_dist_parity = pop_tol,
                 beta_sequence = c(1,1,1,1),
                 beta_weights = c(1,1,1,1),
                 ssdmat = constraints$compact$ssdmat,
                 lambda = lambda,
                 beta = 0,
                 weight_population = constraints$population$weight,
                 weight_compact = constraints$compact$weight,
                 weight_segregation = constraints$segregation$weight,
                 weight_vra = constraints$vra$weight,
                 weight_similar = constraints$similarity$weight,
                 weight_countysplit = constraints$countysplit$weight,
                 weight_partisan = constraints$partisan$weight,
                 weight_minority = constraints$minority$weight,
                 weight_hinge = constraints$hinge$weight,
                 adapt_beta = 'none',
                 adjswap = TRUE,
                 exact_mh = FALSE,
                 adapt_lambda = FALSE,
                 adapt_eprob = FALSE,
                 compactness_measure = constraints$compact$metric,
                 partisan_measure = constraints$partisan$metric,
                 ssd_denom = constraints$compact$ssd_denom,
                 tgt_min = constraints$vra$target_min,
                 tgt_other = constraints$vra$target_other,
                 rvote = constraints$partisan$rvote,
                 dvote = constraints$partisan$dvote,
                 minorityprop = constraints$hinge$minorityprop,
                 verbose = FALSE)

  return(algout$plans + 1)
}