#' Compare the Population Overlap Across Plans at the District Level
#'
#' This implements Crespin's 2005 measure of district continuity, as applied to
#' the geographies represented by a plan, typically precincts or voting districts.
#' This implementation assumes none of the precincts in plan_old or plan_new are
#' split.
#'
#' @param plan_old The reference or original plan to compare against
#' @param plan_new The new plan to compare to the reference plan
#' @param total_pop The total population by precinct This can also take a redist_map
#' object and will use the population in that object. If nothing is provided, it weights
#' all entries in plan equally.
#' @param normalize_rows Default TRUE. Normalize populations by row. If FALSE,
#' normalizes by column. If NULL, does not normalize.
#'
#' @return matrix with length(unique(plan_old)) rows and length(unique(plan_new)) columns
#' @export
#'
#' @references
#' "Using Geographic Information Systems to Measure District Change, 2000-02",
#' Michael Crespin, Political Analysis (2005) 13(3): 253-260
#'
#' @examples
#' set.seed(5)
#' data(iowa)
#' iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 4)
#' plans <- redist_smc(iowa_map, 2)
#' plans_mat <- get_plans_matrix(plans)
#' ov <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map)
#' round(ov, 2)
#'
#' ov_col <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map, normalize_rows = FALSE)
#' round(ov_col, 2)
#'
#' ov_un_norm <- redist.dist.pop.overlap(plans_mat[, 1], plans_mat[, 2],
#' iowa_map, normalize_rows = NULL)
#' round(ov_un_norm, 2)
#'
#' iowa_map_5 <- iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 5)
#' plan_5 <- get_plans_matrix(redist_smc(iowa_map_5, 1))
#' ov4_5 <- redist.dist.pop.overlap(plans_mat[, 1], plan_5, iowa_map)
#' round(ov4_5, 2)
#'
redist.dist.pop.overlap <- function(plan_old, plan_new, total_pop, normalize_rows = TRUE) {
  if (missing(plan_old)) {
    stop('Please pass an argument to `plan_old`.')
  }
  if (missing(plan_new)) {
    stop('Please pass an argument to `plan_new`.')
  }

  if (length(plan_old) != length(plan_new)) {
    stop('plan_old and plan_new must have the same number of entries.')
  }

  if (missing(total_pop)) {
    total_pop <- rep(1L, length(plan_old))
  } else {
    if (inherits(total_pop, 'redist_map')) {
      total_pop <- total_pop[[attr(total_pop, 'pop_col')]]
    }
  }

  out <- tapply(total_pop, list(plan_old, plan_new), sum)
  out[is.na(out)] <- 0

  if (is.null(normalize_rows)) {
    return(out)
  }

  if (normalize_rows) {
    pops <- tapply(total_pop, plan_old, sum)
    out <- out / c(pops)
  } else {
    pops <- tapply(total_pop, plan_new, sum)
    out <- t(t(out) / c(pops))
  }

  out
}

#' Compare the Population Overlap Across Plans at the Precinct Level
#'
#' @param plan_old The reference plan to compare against
#' @param plan_new The new plan to compare to the reference plan
#' @param total_pop The total population by precinct This can also take a redist_map
#' object and will use the population in that object. If nothing is provided, it weights
#' all entries in plan equally.
#' @param weighting Should weighting be done by sum of populations `'s'`, mean of
#' populations `'m'`, geometric mean of populations `'g'`, or none `'n'`
#' @param normalize Should entries be normalized by the total population
#' @param index_only Default is FALSE. TRUE returns only one numeric index, the
#' mean of the upper triangle of the matrix, under the weighting and normalization
#' chosen.
#' @param return_mat Defaults to FALSE, where it returns the summary by row. If
#' TRUE returns matrix with length(plan_old) rows and columns. Ignored if
#' index_only = TRUE.
#'
#' @return numeric vector with length(plan_old) entries
#' @export
#'
#' @examples
#' set.seed(5)
#' data(iowa)
#' iowa_map <- redist_map(iowa, total_pop = pop, pop_tol = 0.01, ndists = 4)
#' plans <- redist_smc(iowa_map, 2, silent = TRUE)
#' plans_mat <- get_plans_matrix(plans)
#' ov_vec <- redist.prec.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map)
#' redist.prec.pop.overlap(plans_mat[, 1], plans_mat[, 2], iowa_map,  weighting = 's',
#' normalize = FALSE, index_only = TRUE)
#'
redist.prec.pop.overlap <- function(plan_old, plan_new, total_pop, weighting = 's',
                                    normalize = TRUE, index_only = FALSE, return_mat = FALSE) {

  weighting <- match.arg(weighting, choices = c('s','m','g','n'))

  if (missing(plan_old)) {
    stop('Please pass an argument to `plan_old`.')
  } else {
    plan_old <- as.matrix(plan_old)
  }
  if (missing(plan_new)) {
    stop('Please pass an argument to `plan_new`.')
  } else {
    plan_new <- as.matrix(plan_new)
  }

  if (length(plan_old) != length(plan_new)) {
    stop('plan_old and plan_new must have the same number of entries.')
  }

  nprec <- length(plan_old)

  if (missing(total_pop)) {
    total_pop <- rep(1L, nprec)
  } else {
    if (inherits(total_pop, 'redist_map')) {
      total_pop <- total_pop[[attr(total_pop, 'pop_col')]]
    }
  }

  total_pop <- matrix(total_pop, nrow = length(total_pop))

  if (weighting == 's') { # weight by sum of ith and jth pop
    cmat <- matrix(rep(total_pop, nprec), ncol = nprec)
    rmat <- t(cmat)
    wts <- (cmat + rmat)
  } else if (weighting == 'g') { # weight by geometric mean of ith and jth pop
    wts <- sqrt(total_pop %*% t(total_pop))
  } else if (weighting == 'm') {# weight by mean of ith and jth pop
    cmat <- matrix(rep(total_pop, nprec), ncol = nprec)
    rmat <- t(cmat)
    wts <- (cmat + rmat)/2
  } else { #if (weighting == 'n') { # no weights (matrix of 1s)
    wts <- matrix(1, nprec, nprec)
  }

  mat_old <- prec_cooccur(plan_old, 1)
  mat_new <- prec_cooccur(plan_new, 1)

  wted_mat <- wts * abs(mat_old - mat_new)

  if (normalize) {
    wted_mat <- wted_mat / sum(total_pop)
  }

  if (index_only) {
    return(mean(wted_mat[upper.tri(wted_mat)]))
  }

  if (return_mat) {
    return(wted_mat)
  }
  apply(wted_mat, MARGIN = 1, FUN = mean)
}
