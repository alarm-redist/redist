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
#' @param adj A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied. Default is NULL.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param centroids A logical indicating if centroids should be plotted. Default is \code{TRUE}.
#' @param edges A logical indicating if edges should connect adjacent centroids. Default is \code{TRUE}.
#' @param boundaries A logical indicating if precinct boundaries should be plotted.
#' @param drop A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.
#' @param title A string title of plot. Defaults to empty string. Optional.
#' @param adjacency Deprecated, use adj. A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied. Default is NULL.
#' @param district_membership Deprecated, use plan. A numeric vector with one row for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#'
#' @return ggplot map
#'
#' @importFrom ggplot2 ggplot geom_sf theme_minimal theme labs aes theme_void
#' @importFrom dplyr filter .data
#' @importFrom sf st_centroid st_coordinates st_as_sf st_linestring st_sfc
#'
#' @examples
#' \dontrun{
#' data(fl25)
#' data(fl25_enum)
#'
#' cds <- fl25_enum$plans[, 5118]
#' redist.map(shp = fl25, plan = cds)
#' }
#'
#' @concept plot
#' @export
redist.map <- function(shp = NULL, adj = NULL, plan = NULL, centroids = TRUE,
                       edges = TRUE, boundaries = TRUE, drop = FALSE,
                       title = '', adjacency, district_membership) {
  .Deprecated(msg = 'Please use redist.plot.map or redist.plot.adj.')
  if (!missing(adjacency)) {
    .Deprecated(new = 'adj', old = 'adjacency')
    adj <- adjacency
  }
  if (!missing(district_membership)) {
    .Deprecated(new = 'plan', old = 'district_membership')
    plan <- district_membership
  }

  # Check inputs
  if (is.null(shp)) {
    stop('Please provide an argument to "shp".')
  }

  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }

  if (!is.null(plan)) {
    if (!any(class(plan) %in% c('numeric', 'integer', 'character'))) {
      stop('Please provide "plan" as a vector.')
    }
    if (nrow(shp) != length(plan)) {
      stop('Arguments "plan" and "shp" do not have same number of precincts.')
    }
  }

  if (!edges & drop) {
    warning('edges FALSE while drop TRUE, assumed edges should be TRUE.')
    edges <- TRUE
  }

  if (drop & is.null(plan)) {
    stop('drop is TRUE but no districts supplied')
  }

  # Extract Centers
  if (edges | centroids) {
    suppressWarnings(centers <- st_centroid(shp))
    st_crs(centers) <- st_crs(shp)
  }

  # Extract Edges
  if (edges) {
    if (!is.null(adj)) {
      nb <- lapply(adj, function(x) {
        x + 1L
      })
    } else {
      adj <- redist.adjacency(shp)
      nb <- lapply(adj, function(x) {
        x + 1L
      })
    }

    edgedf <- tibble(
      start = rep(1:length(nb), lengths(nb)),
      finish = unlist(nb)
    )
    edgedf <- edgedf %>%
      rowwise() %>%
      mutate(i = min(start, finish), j = max(start, finish)) %>%
      select(i, j)
    edgedf <- edgedf[!duplicated(edgedf), ]

    edgedf <- edgedf %>%
      rowwise() %>%
      mutate(geometry = st_sfc(st_linestring(matrix(
        c(
          as.numeric(centers$geometry[[i]]),
          as.numeric(centers$geometry[[j]])
        ),
        nrow = 2,
        byrow = TRUE
      ))))

    suppressWarnings(nb <- sf::st_as_sf(edgedf))
    st_crs(nb) <- st_crs(shp)
  }


  # Drop Edges that cross District Boundaries
  if (drop & !is.null(plan)) {
    nb <- nb %>%
      filter(plan[i] == plan[j])
  }

  # Create Plot
  if (!is.null(plan)) {
    plan <- as.factor(plan)

    if (!is.null(adj)) {
      plan <- as.factor(color_graph(adj, as.integer(plan)))
    }

    plot <- ggplot(shp) +
      geom_sf(aes(fill = plan),
        size = 0.3 * boundaries,
        color = '#444444'
      ) +
      theme_minimal() +
      labs(
        fill = 'District Membership', x = 'Longitude',
        y = 'Latitude', title = title
      ) +
      theme(legend.position = 'bottom')

    if (inherits(shp, 'redist_map')) {
      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      plot <- plot + ggplot2::guides(fill = FALSE) +
        ggplot2::scale_fill_manual(values = PAL)
    }
  } else {
    plot <- ggplot(shp) +
      geom_sf() +
      theme_minimal() +
      labs(x = 'Longitude', y = 'Latitude', title = title)
  }

  if (centroids) {
    plot <- plot + geom_sf(data = centers)
  }

  if (edges) {
    plot <- plot + geom_sf(data = nb)
  }

  # return plot
  return(plot)
}


#' Plot a Map
#'
#' Create a ggplot map. It fills by plan or argument fill. If both are supplied,
#' plan is used as the color and fill as the alpha parameter.
#'
#' @param shp  A SpatialPolygonsDataFrame, sf object, or redist_map. Required.
#' @param adj A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied and needed for coloring. Default is NULL.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param fill A numeric/integer vector with values to color the plot with. Optional.
#' @param fill_label A string title of plot. Defaults to empty string.
#' @param boundaries A logical indicating if precinct boundaries should be plotted.
#' @param title A string title of plot. Defaults to empty string. Optional.
#'
#' @return ggplot map
#'
#' @importFrom ggplot2 ggplot geom_sf theme_minimal theme labs aes theme_void
#' @importFrom dplyr filter .data
#' @importFrom sf st_centroid st_coordinates st_as_sf st_linestring st_sfc
#'
#' @examples
#' \dontrun{
#' data(iowa)
#' redist.plot.map(shp = iowa, plan = iowa$cd_2010)
#' 
#' iowa %>% redist_map(existing_plan = cd_2010) %>% redist.plot.map(shp = ., plan = get_existing(.))
#' }
#'
#' @concept plot
#' @export
redist.plot.map <- function(shp, adj, plan = NULL, fill = NULL, fill_label = '', boundaries = TRUE, title = '') {

  # Check inputs
  if (missing(shp)) {
    stop('Please provide an argument to "shp".')
  }

  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }

  plan <- eval_tidy(enquo(plan), shp)
  if (!is.null(plan)) {
    if (!any(class(plan) %in% c('numeric', 'integer', 'character'))) {
      stop('Please provide "plan" as a vector.')
    }
    if (nrow(shp) != length(plan)) {
      stop('Arguments "plan" and "shp" do not have same number of precincts.')
    }
  }
  fill <- eval_tidy(enquo(fill), shp)

  # Create Plot
  if (!is.null(plan)) {
    if (inherits(shp, 'redist_map') & is.null(fill)) {
      plan <- as.factor(plan)
      adj <- get_adj(shp)
      plan <- as.factor(color_graph(adj, as.integer(plan)))


      plot <- ggplot(shp) +
        geom_sf(aes(fill = plan), size = 0.3 * boundaries, color = '#444444') +
        theme_void() +
        labs(fill = 'District Membership', title = title) +
        theme(legend.position = 'bottom')


      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      plot <- plot + ggplot2::guides(fill = FALSE) +
        ggplot2::scale_fill_manual(values = PAL) 
    } else if(inherits(shp, 'redist_map')) {
      plan <- as.factor(plan)
      adj <- get_adj(shp)
      plan <- as.factor(color_graph(adj, as.integer(plan)))
      
      if(max(fill, na.rm = TRUE) > 1){
        fill <- fill/max(fill)
      }
      
      plot <- ggplot(shp) +
        geom_sf(aes(fill = plan, alpha = fill), size = 0.3 * boundaries, color = '#444444') +
        theme_void() +
        labs(alpha = fill_label, title = title) +
        theme(legend.position = 'bottom')
      
      
      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      plot <- plot + ggplot2::guides(fill = FALSE) +
        ggplot2::scale_fill_manual(values = PAL)  +
        ggplot2::scale_alpha_continuous(range = c(0,1))
    } else {
      if(is.null(fill)){ # plan but no fill
        plot <- ggplot(shp) +
          geom_sf(aes(fill = as.character(plan)), size = 0.3 * boundaries, color = '#444444') +
          theme_void() +
          labs(fill = 'District Membership', title = title) +
          theme(legend.position = 'bottom')
      } else { # plan and fill
        if(max(fill, na.rm = TRUE) > 1){
          fill <- fill/max(fill)
        }
        
        plot <- ggplot(shp) +
          geom_sf(aes(fill = as.character(plan), alpha = fill), size = 0.3 * boundaries, color = '#444444') +
          theme_void() +
          labs(fill = 'District Membership', alpha = fill_label, title = title) +
          theme(legend.position = 'bottom')
        
        if (min(fill, na.rm = T) >= 0 & max(fill, na.rm = T) <= 1) {
          plot <- plot + lims(alpha = c(0, 1))
        }
        
      }

    }
  } else if(!is.null(fill)){ # no plan but fill
    plot <- plot <- ggplot(shp) +
      geom_sf(aes(fill = fill), size = 0.3 * boundaries, color = '#444444') +
      theme_void() +
      labs(fill = fill_label, title = title) +
      theme(legend.position = 'bottom')
    
    if(min(fill, na.rm = TRUE) >= 0 & max(fill, na.rm = TRUE) <= 1){
      plot <- plot + scale_fill_gradient(low = '#ffffff', high = '#08306b',
                          limits = c(0, 1))
    }
    
  } else   {
    plot <- ggplot(shp) +
      geom_sf() +
      theme_void() +
      labs(title = title)
  }

  # return plot
  return(plot)
}

#' Creates a Graph Overlay 
#'
#' @param shp  A SpatialPolygonsDataFrame or sf object. Required.
#' @param adj A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied. Default is NULL.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to remove edges that cross boundaries. Default is \code{NULL}.  Optional.
#' @param centroids A logical indicating if centroids should be plotted. Default is \code{TRUE}.
#' @param drop A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.
#' @param title A string title of plot. Defaults to empty string. Optional.
#'
#' @return ggplot map
#'
#' @importFrom ggplot2 ggplot geom_sf theme_minimal theme labs aes theme_void
#' @importFrom dplyr filter .data
#' @importFrom sf st_centroid st_coordinates st_as_sf st_linestring st_sfc
#' @importFrom rlang eval_tidy enquo
#'
#' @examples
#' \dontrun{
#' data(fl25)
#' data(fl25_enum)
#'
#' cds <- fl25_enum$plans[, 5118]
#' redist.map(shp = fl25, plan = cds)
#' }
#'
#' @concept plot
#' @export
redist.plot.adj <- function(shp = NULL, adj = NULL, plan = NULL, centroids = TRUE,
                       drop = FALSE, title = '') {

  # Check inputs
  if (is.null(shp)) {
    stop('Please provide an argument to "shp".')
  }
  
  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }
  
  plan <- eval_tidy(enquo(plan), shp)
  if (!is.null(plan)) {
    if (!any(class(plan) %in% c('numeric', 'integer', 'character'))) {
      stop('Please provide "plan" as a vector.')
    }
    if (nrow(shp) != length(plan)) {
      stop('Arguments "plan" and "shp" do not have same number of precincts.')
    }
  }
  
  if( inherits(shp, 'redist_map') ){
    if(missing(adj)) {
      adj <- get_adj(shp)
    }
  } else if (missing(adj)){
    adj <- redist.adjacency(shp)
  }

  
  if (drop & is.null(plan)) {
    stop('drop is TRUE but no plan supplied')
  }
  
  # Extract Centers
    suppressWarnings(centers <- st_centroid(shp))
    st_crs(centers) <- st_crs(shp)

  
  # Extract Edges
  
  nb <- lapply(adj, function(x) {
    x + 1L
  })

    
  edgedf <- tibble(
    start = rep(1:length(nb), lengths(nb)),
    finish = unlist(nb)
  )
  edgedf <- edgedf %>%
    rowwise() %>%
    mutate(i = min(start, finish), j = max(start, finish)) %>%
    select(i, j)
  edgedf <- edgedf[!duplicated(edgedf), ]
  
  edgedf <- edgedf %>%
    rowwise() %>%
    mutate(geometry = st_sfc(st_linestring(matrix(
      c(
        as.numeric(centers$geometry[[i]]),
        as.numeric(centers$geometry[[j]])
      ),
      nrow = 2,
      byrow = TRUE
    ))))
  
  suppressWarnings(nb <- sf::st_as_sf(edgedf))
  st_crs(nb) <- st_crs(shp)

  
  
  # Drop Edges that cross District Boundaries
  if (drop) {
    nb <- nb %>%
      filter(plan[i] == plan[j])
  }
  
  # Create Plot
  plot <- ggplot(nb) +
    geom_sf() +
    theme_void()
  
  if (centroids) {
    plot <- plot + geom_sf(data = centers)
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
#' data('fl25')
#' redist.choropleth(shp = fl25, fill = BlackPop / TotPop)
#' DVS <- fl25$obama / (fl25$mccain + fl25$obama)
#' redist.choropleth(shp = fl25, fill = DVS, limit_colors = c('red', 'blue'))
#' }
#'
#' @concept plot
#' @export
redist.choropleth <- function(shp, fill = NULL, fill_label = '', title = '',
                              grad = 1, lwd = 0) {
  .Deprecated(new = 'redist.plot.map')
  # Check inputs
  if (missing(shp)) {
    stop('Please provide an argument to "shp".')
  }
  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }

  fill <- rlang::eval_tidy(rlang::enquo(fill), shp)
  if (!is.null(fill) && !is.numeric(fill) && inherits(shp, 'redist_map')) {
    fill <- as.factor(color_graph(get_adj(shp), as.integer(as.factor(fill))))
  }

  plot <- ggplot(shp, aes(fill = fill)) +
    geom_sf(size = 0, color = NA) +
    theme_void() +
    labs(fill = fill_label, x = 'Longitude', y = 'Latitude', title = title) +
    theme(legend.position = 'bottom')

  if (!is.null(fill)) {
    if (is.numeric(fill)) {
      if (min(fill, na.rm = T) >= 0 & max(fill, na.rm = T) <= 1) {
        plot <- plot + lims(fill = c(0, 1))
      }
    } else {
      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      plot <- plot + ggplot2::scale_fill_manual(values = PAL) +
        ggplot2::guides(fill = FALSE)
    }
  }

  return(plot)
}

globalVariables(c('start', 'finish'))
