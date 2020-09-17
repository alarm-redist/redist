#' Adjacency List functionality for redist
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#'
#' @return Adjacency list
#' @description Creates an adjacency list that is zero indexed with no skips
#' 
#' @export
redist.adjacency <- function(shp){
  
  # Check input
  if(!any(c('sf','SpatialPolygonsDataFrame') %in% class(shp))){
    stop('Please provide "shp" as an sf or sp object.')
  }
  
  # Create the adjacency with spdep function
  adj <- spdep::poly2nb(shp, queen = FALSE)
  
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