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
#' shapefile. It currently computes the Polsby-Popper, Schwartzberg score, Length Width Ratio,
#' Convex Hull score, and Reock score.
#'
#' @usage redist.compactness(shp, district_membership, 
#' measure = c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock"))
#' 
#' @param shp A SpatialPolygonsDataFrame or sf object
#' @param district_membership A string that is the name of a column of district 
#' identifiers in the shp object.
#' @param measure A vector with a string for each measure desired. "PolsbyPopper", 
#' "Schwartzberg", "LengthWidth", "ConvexHull", and "Reock" are implemented.
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
#' 
#' @import sf
#' @import lwgeom
#' @import tidyverse
#' @export
redist.compactness <- function(shp, 
                               district_membership, 
                               measure = c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock")){

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
  match.arg(measure)
  
  # Compute compactness scores
  dists <- unique(shp[[district_membership]])
  nd <-  length(dists)
  
  # Initialize object
  comp <- tibble(districts = dists, PolsbyPopper = rep(NA_real_, nd), Schwartzberg = rep(NA_real_, nd),
                 LengthWidth = rep(NA_real_, nd),  ConvexHull = rep(NA_real_, nd),
                 Reock = rep(NA_real_,nd)) %>% 
    select('districts', measure)
  
  # Compute Specified Scores for provided districts
  for (i in 1:nd){
    united <- st_union(shp[shp[[district_membership]] == dists[i],])
    area <- st_area(united)
    perim <- st_perimeter(united)
    if('PolsbyPopper' %in% measure){  
      comp[['PolsbyPopper']][i] <- 4*pi*(area)/(perim)^2
    }
    if('Schwartzberg' %in% measure){
      comp[['Schwartzberg']][i] <- 1/(perim/(2*pi*sqrt(area/pi)))
    }
    if('LengthWidth' %in% measure){
      bbox <- st_bbox(united)
      ratio <- unname((bbox$xmax - bbox$xmin)/(bbox$ymax - bbox$ymin))
      comp[['LengthWidth']][i] <- ifelse(ratio>1, 1/ratio, ratio) 
    }
    if('ConvexHull' %in% measure){
      cvh <- st_area(st_convex_hull(united))
      comp[['ConvexHull']][i] <- area/cvh
    }
    if('Reock' %in% measure){
      mbc <- st_area(st_minimum_bounding_circle(united))
      comp[['Reock']][i] <- area/mbc
    }
  }
  
  # Return results
  return(comp)
}
