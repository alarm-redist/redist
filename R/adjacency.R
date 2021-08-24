#' Adjacency List functionality for redist
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#' @param plan A numeric vector (if only one map) or matrix with one row
#'
#' @return Adjacency list
#' @description Creates an adjacency list that is zero indexed with no skips
#'
#'
#' @importFrom sf st_relate
#' @export
redist.adjacency <- function(shp, plan){
  # Check input
  if (!any(c('sf','SpatialPolygonsDataFrame') %in% class(shp))) {
    stop('Please provide "shp" as an sf or sp object.')
  }

  # Create the adjacency with spdep function
  # Get standard rooks contiguity
  adj <- suppressMessages(st_relate(shp, shp, pattern = "F***1****"))
  # items contained entirely within ~ even if validly 'rooks' adjacent ~ do not meet this, you need:
  withinadj <- suppressMessages(st_relate(x = shp, pattern = "2121**2*2"))
  adj <- lapply(1:nrow(shp), function(x) c(adj[[x]], withinadj[[x]]))

  # Check for zero indexing
  zero <- min(unlist(adj)) == 0

  # Make zero indexed if not
  min <- min(unlist(adj))
  if(!zero){
  adj <- lapply(adj, function(x){x-min})
  }

  # Check that no numbers are skipped
  # low resolution shp files may result in skips, this fixes most issues
  correct_n <- nrow(shp) == length(unique(unlist(adj)))

  if (!correct_n) {
    warning('At least one precinct had no adjacent precincts.')
  } else if (any(contiguity(adj = adj, group = rep(1, length(adj))) > 1)) {
    warning('All precincts have at least one neighbor, but the graph is disconnected.')
  }


  if (!missing(plan)) {
    cont <- contiguity(adj, plan)
    if (any(cont > 1)) {
      warning(paste0('District', unique(plan[cont>1]), ' was not contiguous.'))
    }
  }

  # return a checked adjacency list
  return(adj)
}

#' Reduce Adjacency List
#'
#' Tool to help reduce adjacency lists for analyzing subsets of maps.
#'
#' @param adj A zero-indexed adjacency list. Required.
#' @param keep_rows row numbers of precincts to keep
#'
#' @return zero indexed adjacency list with max value length(keep_rows) - 1
#'
#' @concept prepare
#' @export
#'
#' @examples
#' data(fl25_adj)
#' redist.reduce.adjacency(fl25_adj, c(2, 3, 4, 6, 21))
#' 
redist.reduce.adjacency <- function(adj, keep_rows) {
    # Check inputs:
    if (!(class(keep_rows) %in% c('numeric', 'integer'))) {
        stop('Please provide "keep_rows" as a numeric or integer vector.')
    }
    if (min(unlist(adj)) != 0) {
        stop('Please provide "adj" as a 0-indexed list.')
    }
    if (max(unlist(adj)) != (length(adj) - 1)) {
        warning('"adj" did not have typical values of 0:(length(adj)-1)')
    }

    # Prep objects for Rcpp
    prec_map = rep(-1L, length(adj))
    #prec_map[keep_rows] = order(keep_rows) - 1L
    prec_map <- dplyr::coalesce(match(1:length(adj), keep_rows) - 1L, -1L)

    # Reduce!
    reduce_adj(adj, prec_map, length(keep_rows))
}


#' Coarsen Adjacency List
#'
#' @param adj A zero-indexed adjacency list. Required.
#' @param groups integer vector of elements of adjacency to group
#'
#' @return adjacency list coarsened
#'
#' @concept prepare
#' @export
redist.coarsen.adjacency <- function(adj, groups) {
    if (min(unlist(adj)) != 0) {
        stop('Please provide "adj" as a 0-indexed list.')
    }
    if (max(unlist(adj)) != (length(adj) - 1)) {
        warning('"adj" did not have typical values of 0:(length(adj)-1)')
    }
    if (length(groups) != length(adj)) {
        stop('groups and adj have sizes which do not conform.')
    }
    if (min(groups) != 0) {
        groups <- groups - min(groups)
    }

    groups <- as.integer(groups)

    coarsen_adjacency(adj, groups)
}
