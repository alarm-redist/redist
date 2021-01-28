#' Adjacency List functionality for redist
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#'
#' @return Adjacency list
#' @description Creates an adjacency list that is zero indexed with no skips
#'
#'
#' @importFrom sf st_relate
#' @export
redist.adjacency <- function(shp){

  # Check input
  if(!any(c('sf','SpatialPolygonsDataFrame') %in% class(shp))){
    stop('Please provide "shp" as an sf or sp object.')
  }

  # Create the adjacency graph
  adj <- sf::st_relate(shp, shp, pattern = "F***1****")

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

  # if there are skips, sink
  while(skip){
    arr <- sort(unique(unlist(adj)))
    index <- min(which(arr != 0:(nrow(shp)-1)))-1
    adj <- lapply(adj, function(x){
      replace(x, list = which(x == arr[index+1]), index)
    })
    skip <- !all(sort(unique(unlist(adj))) == 0:(nrow(shp)-1))
  }

  # return a checked adjacency list
  return(adj)
}

#' Reduce Adjacency List
#'
#' Tool to help reduce adjacency lists for analyzing subsets of maps.
#'
#' @param adjacency  A zero-indexed adjacency list. Required.
#' @param keep_rows row numbers of precincts to keep
#'
#' @return zero indexed adjacency list with max value length(keep_rows) - 1
#' @export
#'
#' @examples \dontrun{
#' data("algdat.p10")
#' redist.reduce.adjacency(algdat.p10$adjlist, c(2, 3, 4, 6, 21))
#' }
redist.reduce.adjacency <- function(adjacency, keep_rows){
  # Check inputs:
  if(!(class(keep_rows) %in% c('numeric', 'integer'))){
    stop('Please provide "keep_rows" as a numeric or integer vector.')
  }
  if(min(unlist(adjacency)) != 0){
    stop('Please provide "adjacency" as a 0-indexed list.')
  }
  if(max(unlist(adjacency))!= (length(adjacency)-1)){
    warning('"adjacency" did not have typical values of 0:(length(adjacency)-1)')
  }

  # Prep objects for Rcpp
  prec_keep <- rep(0L, length(adjacency))
  prec_keep[keep_rows] <- 1L
  keep_rows <- keep_rows - 1
  keep_rows <- as.integer(keep_rows)

  # Reduce!
  return(reduce_adj(adj_list = adjacency, prec_keep = prec_keep,
                    prec_idx = keep_rows))
}
