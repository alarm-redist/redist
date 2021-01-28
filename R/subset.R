#' Subset a shp
#' 
#' Subsets a shp object along with its adjacency. Useful for running smaller analyses
#' on pairs of districts. Provide population, ndist, popcons, and subndist to get proper 
#' population parity constraints on subsets.
#'
#' @param shp  An sf object
#' @param adjacency A zero-indexed adjacency list. Created with 
#' \code{redist.adjacency} if not supplied.
#' @param keep_rows row numbers of precincts to keep. Random submap selected if not supplied.
#' @param population numeric vector with one entry for the population of each precinct.
#' @param ndist integer, number of districts in whole map
#' @param popcons The strength of the hard population constraint.
#' @param subndist integer, number of districts in subset map
#'
#' @return a list containing the following components:
#' \item{shp}{The subsetted shp object}
#' \item{adjacency}{The subsetted adjacency list for shp}
#' \item{keep_rows}{The indices of the rows kept.}
#' \item{subndist}{The number of districts in the subset.}
#' \item{subpopcons}{The new parity constraint for a subset.}
#' @export
#'
redist.subset <- function(shp, adjacency, keep_rows, population, ndist, popcons, subndist){
  if(missing(shp)){
    stop('Please provide an argument to "shp". Use redist.reduce.adjacency to subset adjacency lists.')
  }
  if(!('sf' %in% class(shp))){
    stop('Please provide "shp" as an sf object.')
  }
  
  if(missing(adjacency)){
    adjacency <- redist.adjacency(shp)
  }
  
  if(missing(keep_rows)){
    n <- sample(1:nrow(shp),1)
    keep_rows <- redist.random.subgraph(shp, n, adjacency)$keep_rows
  }
  
  if(!missing(population)&!missing(ndist)&!missing(popcons)&!missing(subndist)){
    parpop <- sum(population)/ndist
    subparpop <- sum(population[keep_rows])/subndist
    subdev <- min(abs(subparpop-parpop*(1-popcons)), abs(subparpop-parpop*(1+popcons)))
    subpopcons <- subdev/subparpop
  } else{
    subndist <- NA_real_
    subpopcons <- NA_real_
  }
  
  
  rlist <- list(shp = shp %>% dplyr::slice(keep_rows),
                adjacency = redist.reduce.adjacency(adjacency, keep_rows = keep_rows),
                keep_rows = keep_rows,
                subndist = subndist,
                subpopcons = subpopcons)
  
  return(rlist)
  
}
