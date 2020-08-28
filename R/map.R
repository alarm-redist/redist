##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/06/20
## Date Modified: 2020/06/20
## Purpose: R function to make map plot
##############################################

#' Creates a map with optional graph overlay
#'
#' @param shp  A SpatialPolygonsDataFrame or sf object. Required.
#' @param district_membership A numeric vector with one row for each precinct in shp. 
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param centroids A logical indicating if centroids should be plotted. Default is \code{TRUE}.
#' @param edges A logical indicating if edges should connect adjacent centroids. Default is \code{TRUE}.
#' @param drop A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.
#' @param title A string title of plot. Defaults to empty string. Optional.
#'
#' @return ggplot map
#' 
#' @importFrom spdep nb2lines
#' @importFrom ggplot2 ggplot geom_sf theme_minimal theme labs aes
#' @importFrom dplyr filter .data
#' @importFrom sf st_centroid st_coordinates st_as_sf
#' 
#' @examples
#' \dontrun{
#' library(redist)
#' data("fl25")
#' data("algdat.p10")
#' cds <- algdat.p10$cdmat[,100]
#' redist.map(shp = fl25, district_membership = cds)
#' }
#' 
#' @export
redist.map <- function(shp = NULL, district_membership = NULL, centroids = TRUE, 
                       edges = TRUE, drop = FALSE, title = ""){
  
  # Check inputs
  if(is.null(shp)){
    stop('Please provide an argument to "shp".')
  }

  if('SpatialPolygonsDataFrame' %in% class(shp)){
    shp <- shp %>%  st_as_sf()
  } else if(!('sf' %in% class(shp))){
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }
  
  if(!any(class(district_membership) %in% c('numeric', 'integer', 'character'))){
    stop('Please provide "district_membership" as a vector.')
  }
  
  if(nrow(shp) != length(district_membership)){
    stop('Arguments "district_membership" and "shp" do not have same number of precincts.')
  }
  
  if(!edges & drop){
    warning('edges FALSE while drop TRUE, assumed edges should be TRUE.')
    edges <- TRUE
  }
  
  if(drop & is.null(district_membership)){
    stop('drop is TRUE but no districts supplied')
  }
  
  
  # Extract Centers
  if(edges | centroids){
    suppressWarnings(centers <- st_centroid(shp))
  }
  
  # Extract Edges
  if(edges){
    nb <- spdep::poly2nb(shp, queen = FALSE)
    nb <- spdep::nb2lines(nb, coords = sf::st_coordinates(centers))
    nb <- sf::st_as_sf(nb)
  }
  
  # Drop Edges that cross District Boundaries
  if(drop&!is.null(district_membership)){
    nb <- nb %>% 
      filter(district_membership[i] == district_membership[j])
  }
  
  # Create Plot
  if(!is.null(district_membership)){
    district_membership <- as.character(district_membership)
    plot <- shp %>% ggplot() +
      geom_sf(aes(fill = district_membership)) +
      theme_minimal()     +
      labs(fill = 'District Membership', x = 'Longitude', y = 'Latitude', title = title) +
      theme(legend.position = "bottom")
  } else {
    plot <- shp %>% ggplot() +
      geom_sf() +
      theme_minimal() +
      labs(x = 'Longitude', y = 'Latitude', title = title)
  }
  
  if(centroids){
    plot <- plot + 
      geom_sf(data = centers)
  }
  
  if(edges){
    plot <- plot +
      geom_sf(data = nb)
  }

  # return plot
  return(plot)
}



