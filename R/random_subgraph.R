#' Preprocess shapefile arguments
#'
#' @param shp sf or SpatialPolygonsDataFrame
#'
#' @return shp as sf object
#'
#' @examples \dontrun{
#' data("fl25")
#' preproc.shp(fl25)
#' }
#' @noRd
preproc.shp <- function(shp){
  if(!is.null(shp)){
    if('SpatialPolygonsDataFrame' %in% class(shp)){
      shp <- shp %>%  st_as_sf()
    } else if(!('sf' %in% class(shp))){
      stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
    }
  }
  return(shp)
}

#' Preprocess adjacency arguments
#'
#' @param shp sf or SpatialPolygonsDataFrame
#' @param adj adjacency list for shp
#'
#' @return adjacency list
#'
#' @examples \dontrun{
#' data("fl25")
#' preproc.adj(fl25,NULL)
#' }
#' 
#' @noRd
preproc.adj <- function(shp, adj){
  if(is.null(adj)){
    adj <- redist.adjacency(shp)
  } else if(nrow(shp) != length(adj)){
    stop('Dimension of shp and adj do not match.')
  }
  return(adj)
}



#' Return a random subgraph of a shape
#'
#' @description `random.subgraph` returns a random subset of the shp provided
#'
#' @details Snowball sampling with backtracking
#'
#' @param shp sf object or SpatialPolygonsDataFrame
#' @param n number of edges to sample. n must be a positive integer.
#' @param adj Optional. zero indexed adjacency list.
#'
#' @return sf dataframe with n rows
#' @export
#' @importFrom dplyr union setdiff slice %>%
#'
#'
redist.random.subgraph <- function(shp, n, adj = NULL){
  # Check input:
  shp <- preproc.shp(shp)
  adj <- preproc.adj(shp, adj)

  if(n < 1){
    stop('Please provide "n" as a positive integer')
  }
  if(n > nrow(shp)){
    stop('"n" has more entries than shp. Please provide smaller "n"')
  }

  # create helper objects
  index <- rep(NA_real_, n)
  index[1] <- sample.int(nrow(shp), 1)
  candidates <- adj[[index[1]]]+1

  i <- 1
  while(i < n){
    # increment
    i <- i + 1
    # pick one from connected objects
    index[i] <- sample(candidates, 1)

    # add new candidate options and ignore existing ones
    candidates <- dplyr::union(candidates, adj[[index[i]]]+1) %>% dplyr::setdiff(index)
  }

  rlist <- list(shp = shp %>% dplyr::slice(sort(index)),
                keep_rows = sort(index)

  )
  return(rlist)

}

