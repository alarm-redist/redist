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
                 popvec = total_pop,
                 constraints = as.list(constraints),
                 nsims = nsims,
                 eprob = eprob,
                 pct_dist_parity = pop_tol,
                 beta_sequence = c(1 ,1, 1, 1),
                 beta_weights = c(1, 1, 1, 1),
                 lambda = lambda,
                 beta = 0,
                 adapt_beta = 'none',
                 adjswap = TRUE,
                 exact_mh = FALSE,
                 adapt_lambda = FALSE,
                 adapt_eprob = FALSE,
                 verbose = FALSE)

  algout$plans + 1
}
