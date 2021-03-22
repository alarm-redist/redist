#' Flip Constraints Helper
#'
#' Creates a constraints list for use with  \code{redist_flip}.
#'
#'
#' @param map A \code{\link{redist_map}} object.
#' @param constraint character vector of constraints to use. Current accepted are "compact",
#' "population", "countysplit", "hinge", "vra", "minority", "similarity", "partisan", and "segregation."
#' Defaults to compact.
#' @param constraintweight corresponding weights to use with constraint. Weights must be nonzero if provided. Defaults to
#' a weak compactness constraint
#' @param compactness_metric character with "edges-removed", "polsby-popper", or "fryer-holden". Default is edges-removed.
#' @param areas areas to use with compact:polsby-popper. Computed from map if not provided and needed.
#' @param borderlength_mat border lengths to use with compact:polsby-popper.
#' Computed from map if not provided and needed.
#' @param ssdmat squared distance matrix to use with compact:fryer-holden. Computed from map if not provided and needed.
#' @param ssd_denom Denominator for compact:fryer-holden. Defaults to 1.0
#' @param counties tidy evaluated column with counties in map.
#' @param partisan_metric character with "proportional-representation" or "efficiency-gap"
#' @param rvote numeric vector of votes for Republicans. used with partisan constraint.
#' @param dvote numeric vector of votes for Democrats used with partisan constraint.
#' @param group_pop tidy evaluated column with group populations in map. Used with hinge, segregation, vra, and minority.
#' @param target_min numeric with higher minority target for vra. Defaults to 0.55.
#' @param target_other numeric with higher minority target for vra. Defaults to 0.25.
#' @param minorityprop minority proportions to target. Used by minority and hinge.Defaults to c(0.55) if none specified.
#'
#' @return list with an entry for each constraint usable with redist_flip()
#' @export 
#'
#' @examples \dontrun{
#' data(iowa)
#' iowa_map <- redist_map(iowa, existing_plan = cd_2010)
#'
#' cons <- cons <- flip_constraints_helper(
#'   map = iowa_map, constraint = c('compact', 'vra'),
#'   constraintweight = c(0.05, 10), target_min = 0.05,
#'   target_other = 0.01, group_pop = black
#' )
#'
#' redist_flip(iowa_map, 10, constraints = cons)
#' }
flip_constraints_helper <- function(map,
                                    constraint = 'compact',
                                    constraintweight = 0.6,
                                    compactness_metric = 'edges-removed',
                                    areas,
                                    borderlength_mat,
                                    ssdmat,
                                    ssd_denom,
                                    counties = NULL,
                                    partisan_metric = 'efficiency-gap',
                                    rvote,
                                    dvote,
                                    group_pop,
                                    target_min = 0.55,
                                    target_other = 0.25,
                                    minorityprop) {
  if (missing(map)) {
    stop('Please provide a redist_map object to map.')
  }

  constraints_list <- process_flip_constr(
    constraints = list(),
    group_pop = rep(0, nrow(map)),
    counties = rep(1, nrow(map))
  )

  if (is.null(constraint)) {
    constraints_list$compact$weight <- 0

    return(constraints_list)
  }

  match.arg(constraint, choices = names(constraints_list), several.ok = TRUE)

  if (any(constraint %in% c('hinge', 'vra', 'segregation', 'minority'))) {
    group_pop <- eval_tidy(enquo(group_pop), map)
    if (is.null(group_pop)) {
      stop('column matching group_pop could not be found in object map.')
    }
  }

  if (length(constraint) != length(constraintweight)) {
    stop('constraint and constraintweight must have the same number of entries.')
  }

  constraints_list <- process_flip_constr(
    constraints = list(),
    group_pop = rep(0, nrow(map)),
    counties = rep(1, nrow(map))
  )

  if ('compact' %in% constraint) {
    constraints_list$compact$weight <- constraintweight[which(constraint == 'compact')]

    match.arg(arg = compactness_metric, c('edges-removed', 'fryer-holden', 'polsby-popper'))

    if (compactness_metric == 'fryer-holden') {
      if (missing(ssdmat)) {
        suppressWarnings(ssdmat <- st_distance(st_centroid(map), st_centroid(map)))
      }
      constraints_list$compact$ssdmat <- as.matrix(ssdmat)
      if (missing(ssd_denom)) {
        ssd_denom <- 1.0
      }
      constraints_list$compact$ssd_denom
    }

    if (compactness_metric == 'polsby-popper') {
      if (missing(areas)) {
        areas <- st_area(map)
      }
      constraints_list$compact$areas <- areas
      if (missing(borderlength_mat)) {
        borderlength_mat <- as.matrix(redist.prep.polsbypopper(map, get_adj(map)))
      }
      constraints_list$compact$borderlength_mat <- borderlength_mat
    }
  }

  if ('population' %in% constraint) {
    constraints_list$population$weight <- constraintweight[which(constraint == 'compact')]
  }

  if ('countysplit' %in% constraint) {
    constraints_list$countysplit$weight <- constraintweight[which(constraint == 'countysplit')]
    counties <- eval_tidy(enquo(counties), map)
    if (is.null(counties)) {
      stop('column matching counties could not be found in object map.')
    } else {
      counties <- redist.county.id(counties)
      components <- contiguity(adj, counties)
      if (any(components > 1)) {
        warning('counties are not contiguous. Expect additional splits.')
        counties <- redist.county.id(redist.county.relabel(adj, counties))
      }
    }
    constraints_list$countysplit$counties <- counties
  }

  if ('vra' %in% constraint) {
    constraints_list$vra$weight <- constraintweight[which(constraint == 'vra')]

    if (!class(target_min) %in% c('numeric', 'integer')) {
      stop('target_min must be a numeric or integer.')
    }
    if (!class(target_other) %in% c('numeric', 'integer')) {
      stop('target_other must be a numeric or integer.')
    }

    constraints_list$vra$target_min <- target_min
    constraints_list$vra$target_other <- target_other
    constraints_list$vra$group_pop <- group_pop
  }


  if ('hinge' %in% constraint) {
    constraints_list$hinge$weight <- constraintweight[which(constraint == 'hinge')]
    if (missing(minorityprop)) {
      minorityprop <- c(0.55)
    }
    constraints_list$minority$minorityprop <- minorityprop

    constraints_list$minority$group_pop <- group_pop
  }

  if ('minority' %in% constraint) {
    constraints_list$minority$weight <- constraintweight[which(constraint == 'minority')]

    if (missing(minorityprop)) {
      minorityprop <- c(0.55)
    }
    constraints_list$minority$minorityprop <- minorityprop

    constraints_list$minority$group_pop <- group_pop
  }
  if ('segregation' %in% constraint) {
    constraints_list$segregation$weight <- constraintweight[which(constraint == 'segregation')]
    constraints_list$segregation$group_pop <- group_pop
  }
  if ('similarity' %in% constraint) {
    constraints_list$similarity$weight <- constraintweight[which(constraint == 'similarity')]
  }
  if ('partisan' %in% constraint) {
    constraints_list$partisan$weight <- constraintweight[which(constraint == 'partisan')]
    if (missing(rvote)) {
      stop('Please provide an argument to rvote with constraint partisan.')
    }
    if (missing(dvote)) {
      stop('Please provide an argument to dvote with constraint partisan.')
    }

    match.args(partisan_metric, c('proportional-representation', 'efficiency-gap'))

    constraints_list$partisan$metric <- partisan_metric
  }
}
