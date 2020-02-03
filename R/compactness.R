##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/01/20
## Date Modified: 2020/01/22
## Purpose: R function to compute compactness
##############################################


#' Calculate compactness measures for a set of districts
#' 
#' \code{redist.compactness} is used to compute different compactness statistics for a 
#' shapefile. It currently computes the Polsby-Popper, Schwartzberg score, Length Width Ratio,
#' Convex Hull score, and Reock score.
#'
#' @usage redist.compactness(shp, district_membership, 
#' measure = c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark"))
#' 
#' @param shp A SpatialPolygonsDataFrame or sf object
#' @param district_membership A string that is the name of a column of district 
#' identifiers in the shp object.
#' @param measure A vector with a string for each measure desired. "PolsbyPopper", 
#' "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", and "BoyceClark" are implemented. Defaults to 
#' all implemented measures.
#' @param population A numeric vector with the population for every observation. Defaults to NULL. Is
#' only necessary when "FryerHolden" is used for measure. Defaults to NULL.
#' @param nloop A numeric to specify loop number. Defaults to NA_real_.
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
#' box <- box %>% mutate(cds = 1, pop = 10)
#' # Run redist.compactness
#' redist.compactness(box, "cds")
#' }
#' 
#' @import sf
#' @import lwgeom
#' @import tidyverse
#' @export
redist.compactness <- function(shp, 
                               district_membership, 
                               measure = c("PolsbyPopper", "Schwartzberg", "LengthWidth", "ConvexHull", "Reock", "BoyceClark"),
                               population = NULL, nloop = NA_real_){

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

  if('FryerHolden' %in% measure & !(population %in% names(shp))){
    stop('Please provide "population" as the name of a column in "shp".')
  }
  
  if(class(nloop) != 'numeric'){
    stop('Please provide "nloop" as a numeric.')
  }
  
  # Compute compactness scores
  dists <- unique(shp[[district_membership]])
  nd <-  length(dists)
  
  # Initialize object
  comp <- tibble(districts = dists, PolsbyPopper = rep(NA_real_, nd), Schwartzberg = rep(NA_real_, nd),
                 LengthWidth = rep(NA_real_, nd),  ConvexHull = rep(NA_real_, nd),
                 Reock = rep(NA_real_,nd), BoyceClark = rep(NA_real_, nd), FryerHolden = rep(NA_real_, nd), 
                 nloop = rep(nloop, nd)) %>% 
    select(districts, measure, nloop)
  
  # Compute Specified Scores for provided districts
  for (i in 1:nd){
    united <- st_union(shp[shp[[district_membership]] == dists[i],])
    area <- st_area(united)
    
    if(is.null(st_crs(united$EPSG))){
      perim <- st_length(st_cast(st_cast(united, 'POLYGON'),'LINESTRING'))
    } else if (st_is_longlat(united)){
      perim <- st_length(united)
    } else {
      perim <- st_perimeter(united)
    }
    
    #perim <- ifelse(is.null(st_crs(united)$EPSG), , st_perimeter(united))
    if('PolsbyPopper' %in% measure){  
      comp[['PolsbyPopper']][i] <- 4*pi*(area)/(perim)^2
    }
    if('Schwartzberg' %in% measure){
      comp[['Schwartzberg']][i] <- (perim/(2*pi*sqrt(area/pi)))
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
    if('BoyceClark' %in% measure){
      suppressWarnings(center <- st_coordinates(st_centroid(united)))
      suppressWarnings(if(!st_within(center, united, sparse = F)[[1]]){
        center <- st_point_on_surface(united)
      })
      bbox <- st_bbox(united)
      max_dist <- sqrt((bbox$ymax-bbox$ymin)^2+(bbox$xmax-bbox$xmin)^2)
      st_crs(united) <- NA
      
      x_list <- center[1] + max_dist*cos(seq(0,15)*pi/8)
      y_list <- center[2] + max_dist*sin(seq(0,15)*pi/8)
      
      plot <- ggplot(united) + geom_sf()
      radials <- rep(NA_real_, 16)
      for(angle in 1:16){
        line <- data.frame(x = c(x_list[angle],center[1]), y = c(y_list[angle], center[2])) %>% 
          st_as_sf(coords = c('x','y'))  %>% st_coordinates() %>% st_linestring()
        radials[angle] <- max(0, dist(st_intersection(line, united)))
        plot <- plot + geom_sf(data = line)
      }
     comp[['BoyceClark']][i] <- 1 - (sum(abs(radials/sum(radials)*100-6.25))/200) 
    }
    if('FryerHolden' %in% measure){
      suppressWarnings(shp_subset <- st_coordinates(st_centroid(shp[[shp[shp[[district_membership]] == dists[i],]]]))
      dist_sqr <- as.matrix(dist(shp_subset))^2
      pop <- shp[[population]][shp[[district_membership == i,]]]
      pop <- pop*t(matrix(rep(pop,nd),nd))
      comp[['FryerHolden']][i] <- sum(pop*dist_sqr)
    }
  }
  
  # Return results
  return(comp)
}
