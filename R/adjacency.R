#' Adjacency List functionality for redist
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#' @param plan A numeric vector (if only one map) or matrix with one row 
#' @param district_membership Deprecated -- Use plan. A numeric vector (if only one map) or matrix with one row. 
#' for each precinct and one column for each map. Optional. Checks for contiguity within 
#' districts if provided.
#'
#' @return Adjacency list
#' @description Creates an adjacency list that is zero indexed with no skips
#' 
#' 
#' @importFrom sf st_relate
#' @export
redist.adjacency <- function(shp, plan, district_membership){
  if(!missing(district_membership)){
    plan <- district_membership
    .Deprecated('plan')
  }
  # Check input
  if(!any(c('sf','SpatialPolygonsDataFrame') %in% class(shp))){
    stop('Please provide "shp" as an sf or sp object.')
  }
  
  # Create the adjacency with spdep function
  # Get standard rooks contiguity
  adj <- sf::st_relate(shp, shp, pattern = "F***1****")
  # items contained entirely within ~ even if validly 'rooks' adjacent ~ do not meet this, you need:
  withinadj <- st_relate(x = shp, pattern = "2121**2*2")
  adj <- lapply(1:nrow(shp), function(x){c(adj[[x]], withinadj[[x]])})
  
  # Check for zero indexing
  zero <- min(unlist(adj)) == 0
  
  # Make zero indexed if not
  min <- min(unlist(adj))
  if(!zero){
  adj <- lapply(adj, function(x){x-min})
  }
  
  # Check that no numbers are skipped
  # low resolution shp files may result in skips, this fixes most issues
  skip <- !all(sort(unique(unlist(adj))) == 0:(nrow(shp)-1))
  correct_n <- nrow(shp) == length(unique(unlist(adj))) 
  
  if(skip){
    warning('At least one precinct had no adjacent precincts.')
  }
  
  
  if(!missing(plan)){
    cont <- contiguity(adj, plan)
    if(any(cont > 1 )){
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
#' @param adjacency  Deprecated. Use adj. A zero-indexed adjacency list. 
#'
#' @return zero indexed adjacency list with max value length(keep_rows) - 1
#' @export
#'
#' @examples \dontrun{
#' data("algdat.p10")
#' redist.reduce.adjacency(algdat.p10$adjlist, c(2, 3, 4, 6, 21))
#' }
redist.reduce.adjacency <- function(adj, keep_rows, adjacency){
  
  if(~missing(adjacency)){
    adj <- adjacency
    .Deprecated('adj', old = 'adjacency')
  }
  # Check inputs:
  if(!(class(keep_rows) %in% c('numeric', 'integer'))){
    stop('Please provide "keep_rows" as a numeric or integer vector.')
  }
  if(min(unlist(adj)) != 0){
    stop('Please provide "adj" as a 0-indexed list.')
  }
  if(max(unlist(adj))!= (length(adj)-1)){
    warning('"adj" did not have typical values of 0:(length(adj)-1)')
  }
  
  # Prep objects for Rcpp
  prec_keep <- rep(0L, length(adj))
  prec_keep[keep_rows] <- 1L
  keep_rows <- keep_rows - 1
  keep_rows <- as.integer(keep_rows)
  
  # Reduce!
  return(reduce_adj(adj_list = adj, prec_keep = prec_keep, 
                    prec_idx = keep_rows))
}

#' Coarsen Adjacency List
#'
#' @param adj A zero-indexed adjacency list. Required.
#' @param groups integer vector of elements of adjacency to group
#' @param adjacency Deperecated -- use adj. A zero-indexed adjacency list
#'
#' @return adjacency list coarsened
#' @export
redist.coarsen.adjacency <- function(adj, groups, adjacency){
  if(!missing(adjacency)){
    .Deprecated('adj',  old = 'adjacency')
      adj <- adjacency
  }
  
  if(min(unlist(adj)) != 0){
    stop('Please provide "adj" as a 0-indexed list.')
  }
  if(max(unlist(adj))!= (length(adj)-1)){
    warning('"adj" did not have typical values of 0:(length(adj)-1)')
  }
  if(length(groups) != length(adj)){
    stop('groups and adj have sizes which do not conform.')
  }
  if(min(groups) != 0){
    groups <- groups - min(groups)
  }
  
  groups <- as.integer(groups)
  
  return(coarsen_adjacency(adj, groups))

}
