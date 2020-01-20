##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/01/20
## Date Modified: 2020/01/20
## Purpose: R function to compute compactness
##############################################


#' Calculate compactness measures for a set of districts
#' 
#' \code{redist.compactness} is used to compute different compactness statistics for a 
#' shapefile. It currently computes the Polsby-Popper score.
#'
#' @usage redist.compactness(shp, district_membership, measure = "polsbypopper")
#' 
#' @param shp A SpatialPolygonsDataFrame or sf object
#' @param district_membership A string that is the name of a column of district 
#' identifiers in the shp object.
#' @param measure A vector with strings for each measure desired. Only "polsbypopper"
#' is currently implemented.
#' 
#' @details This function computes specified compactness scores for a map.  If there is more than one
#' shape specified for a single district, it combines them and computes one score.
#' 
#' @return A tibble with a column that specifies the district and a column for 
#' each specified measure.
#' 
#' @example
#' \dontrun{
#' library(sf)
#' library(lwgeom)
#' library(redist)
#' library(tidyverse)
#' #Create (or load) a shapefile, in this case the unit square 
#' box <-  rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)) %>% list() %>%  
#' st_polygon() %>% st_sfc() %>% st_as_sf()
#' # Index the congressional districts
#' box <- box %>% mutate(cds = 1)
#' # Run redist.compactness
#' redist.compactness(box, "cds")
#' }
#' 
#' @export
redist.compactness <- function(shp, district_membership, measure = "polsbypopper"){
  # Check for depends
  require(sf)
  require(lwgeom)
  require(tidyverse)
  require(redist)
  # Check Inputs
  if('SpatialPolygonsDataFrame' %in% class(shp)){
    shp <- shp %>%  st_as_sf()
  } else if(!('sf' %in% class(shp))){
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }
  
  if(class(district_membership) != 'character'){
    stop('Please provide "district_membership" as a string.')
  }
  
  if(!(district_membership %in% names(shp))){
    stop('Please provide "district_membership" as the name of a column in "shp".')
  }
  match.arg(measure, choices = c('polsbypopper'))
  
  # Compute compactness scores
  dists <- unique(shp[[district_membership]])
  nd <-  length(dists)
  
  # Initialize object
  comp <- tibble(districts = dists, polspop = rep(NA_real_, nd))
  
  # Compute Polsby Poppers
  for (i in 1:nd){
    united <- st_union(shp[shp[[district_membership]] == dists[i],]) 
    comp[i,2] <- 4*pi*(st_area(united))/(st_perimeter(united))^2
  }
  
  # Return results
  return(comp)
}
