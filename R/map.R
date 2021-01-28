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
#' @param adjacency A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied. Default is NULL.
#' @param district_membership A numeric vector with one row for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param centroids A logical indicating if centroids should be plotted. Default is \code{TRUE}.
#' @param edges A logical indicating if edges should connect adjacent centroids. Default is \code{TRUE}.
#' @param drop A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.
#' @param title A string title of plot. Defaults to empty string. Optional.
#'
#' @return ggplot map
#'
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
redist.map <- function(shp = NULL, adjacency = NULL, district_membership = NULL, centroids = TRUE,
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

  if(!is.null(district_membership)){
    if(!any(class(district_membership) %in% c('numeric', 'integer', 'character'))){
      stop('Please provide "district_membership" as a vector.')
    }
    if(nrow(shp) != length(district_membership)){
      stop('Arguments "district_membership" and "shp" do not have same number of precincts.')
    }

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
    st_crs(centers) <- st_crs(shp)
  }

  # Extract Edges
  if(edges){
    if(!is.null(adjacency)){
      nb <- lapply(adjacency, function(x){x+1L})
      class(nb) <- 'nb'
    } else{
      adjacency <- redist.adjacency(shp)
      nb <- lapply(adjacency, function(x){x+1L})
      class(nb) <- 'nb'
    }

    #nb <- spdep::nb2lines(nb, coords = sf::st_coordinates(centers))
    #suppressWarnings(nb <- sf::st_as_sf(nb))
    #st_crs(nb) <- st_crs(shp)
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


#' Creates a Choropleth
#'
#' @details Creates a basic choropleth for a provided shp with value. Recommended to
#' normalize data to avoid absolute values, in most use cases.
#'
#' @param shp  A SpatialPolygonsDataFrame or sf object. Required.
#' @param fill A numeric/integer vector with values to color the plot with. Optional.
#' @param fill_label A string title of plot. Defaults to empty string. Optional.
#' @param title A string title of plot. Defaults to empty string. Optional.
#' @param limit_colors A length two string vector with two colors, either hex or ggplot color names.
#' @param grad Number of colors to make a gradient with. Accepts values of 1 or 2.
#' @param lwd Line width. Defaults to 0
#'
#' @return ggplot map
#'
#' @importFrom ggplot2 ggplot geom_sf theme_minimal theme labs aes scale_fill_gradient2
#' @importFrom ggplot2 theme_void scale_fill_gradient
#' @importFrom dplyr filter .data
#'
#' @examples
#' \dontrun{
#' data("fl25")
#' redist.choropleth(shp = fl25, fill = fl25$BlackPop/fl25$TotPop)
#' DVS <- fl25$obama/(fl25$mccain+fl25$obama)
#' redist.choropleth(shp = fl25, fill = DVS, limit_colors = c('red', 'blue'))
#' }
#'
#' @export
redist.choropleth <- function(shp, fill = NULL, fill_label = "", title = "",
                              limit_colors = NULL, grad = 1, lwd = 0){

  # Check inputs
  if(missing(shp)){
    stop('Please provide an argument to "shp".')
  }
  if('SpatialPolygonsDataFrame' %in% class(shp)){
    shp <- shp %>%  st_as_sf()
  } else if(!('sf' %in% class(shp))){
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }
  if(!is.null(fill)){
    if(nrow(shp) != length(fill)){
      stop('Arguments "fill" and "shp" do not have same number of precincts.')
    }
  }

  plot <- shp %>% ggplot(aes(fill = fill)) +
    geom_sf(lwd = 0, color = NA) +
    theme_void() +
    labs(fill = fill_label, x = 'Longitude', y = 'Latitude', title = title) +
    theme(legend.position = "bottom")

  if(!is.null(fill)){
    l1 <- min(fill)
    l2 <- max(fill)
    mp <- mean(fill)

    if(l1 >= 0 & l2 <= 1){
      l1 <- 0
      l2 <- 1
      mp <- 0.5
    }

    if(is.null(limit_colors)){
      limit_colors <- c('#ffffff', '#08306b')
    }

    if(grad == 1){
      plot <- plot +
        scale_fill_gradient(low = limit_colors[1], high = limit_colors[2],
                            limits = c(l1,l2))
    } else if(grad == 2){
      plot <- plot +
        scale_fill_gradient2(low = limit_colors[1], high = limit_colors[3],
                             midpoint = mp, mid = limit_colors[2], limits = c(l1,l2))
    }

  }

  return(plot)
}

