##############################################
## Author: Christopher T Kenny
## Institution: Harvard University
## Date Created: 2020/06/20
## Date Modified: 2020/06/20
## Purpose: R function to make map plot
##############################################


#' Plot a Map
#'
#' Create a ggplot map. It fills by plan or argument fill. If both are supplied,
#' plan is used as the color and fill as the alpha parameter.
#'
#' @param shp  A SpatialPolygonsDataFrame, sf object, or redist_map. Required.
#' @param adj A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied and needed for coloring. Default is NULL.
#' @param plan \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} A numeric
#'   vector with one entry for each precinct in shp. Used to color the
#'   districts. Default is \code{NULL}.  Optional.
#' @param fill \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} A
#'   numeric/integer vector with values to color the plot with. Optional.
#' @param fill_label A string title of plot. Defaults to the empty string
#' @param zoom_to \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} An
#'   indexing vector of units to zoom the map to.
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
#' data(iowa)
#' redist.plot.map(shp = iowa, plan = iowa$cd_2010)
#'
#' iowa_map = redist_map(iowa, existing_plan = cd_2010)
#' redist.plot.map(iowa_map, fill=dem_08/tot_08, zoom_to=(cd_2010 == 1))
#'
#' @concept plot
#' @export
redist.plot.map <- function(shp, adj, plan = NULL, fill = NULL, fill_label = '',
                            zoom_to=NULL, boundaries=is.null(fill), title='') {
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
      ndists = length(unique(plan))
      if (ndists > 6)
          plan <- as.factor(color_graph(adj, as.integer(plan)))


      plot <- ggplot(shp) +
        geom_sf(aes(fill = plan), size = 0.3 * boundaries, color = '#444444') +
        theme_void() +
        labs(fill = 'District', title = title) +
        theme(legend.position = 'bottom')


      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      if (ndists > 6) {
          plot <- plot + ggplot2::guides(fill = 'none') +
            ggplot2::scale_fill_manual(values = PAL)
      } else {
          plot <- plot + ggplot2::scale_fill_manual(values = PAL)
      }
    } else if(inherits(shp, 'redist_map')) {
      plan <- as.factor(plan)
      adj <- get_adj(shp)
      ndists = length(unique(plan))
      if (ndists > 6)
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
      plot <- plot + ggplot2::guides(fill = 'none') +
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
  } else if (!is.null(fill)) { # no plan but fill
    recolor <- FALSE
    if (inherits(shp, "redist_map") && (is.character(fill) || is.factor(fill))) {
      adj <- get_adj(shp)
      nlevels <- length(unique(fill))
      fill <- as.factor(fill)
      recolor <- TRUE
      if (nlevels > 6) {
        fill <- as.factor(color_graph(adj, as.integer(fill)))
      }
    }

    plot <- ggplot(shp) +
      geom_sf(aes(fill = fill), size = 0.3 * boundaries, color = '#444444') +
      theme_void() +
      labs(fill = fill_label, title = title) +
      theme(legend.position = 'bottom')

    if (recolor) {
      PAL <- c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
      plot <- plot + ggplot2::guides(fill = 'none') +
        ggplot2::scale_fill_manual(values = PAL)
    }

    if(is.numeric(fill) && min(fill, na.rm = TRUE) >= 0 &&
           max(fill, na.rm = TRUE) <= 1){
        plot <- plot + ggplot2::scale_fill_gradient(low = '#ffffff',
                                                    high = '#08306b',
                                                    limits = c(0, 1))
    }

  } else   {
    plot <- ggplot(shp) +
      geom_sf() +
      theme_void()
  }

  zoom_to = eval_tidy(enquo(zoom_to), shp)
  if (!is.null(zoom_to)) {
      bbox = sf::st_bbox(sf::st_geometry(shp)[zoom_to])
      plot <- plot + ggplot2::coord_sf(xlim=c(bbox$xmin, bbox$xmax),
                                       ylim=c(bbox$ymin, bbox$ymax))
  }


  plot <- plot + labs(title = title)
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
#' @param plot_shp A logical indicating if the shp should be plotted under the
#'   graph. Default is \code{TRUE}.
#' @param zoom_to \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} An
#'   indexing vector of units to zoom the map to.
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
#' data(iowa)
#' redist.plot.adj(shp = iowa, plan = iowa$cd_2010)
#'
#' @concept plot
#' @export
redist.plot.adj <- function(shp = NULL, adj = NULL, plan = NULL, centroids = TRUE,
                       drop = FALSE, plot_shp=TRUE, zoom_to=NULL, title = '') {

  # Check inputs
  if (is.null(shp)) {
    stop('Please provide an argument to "shp".')
  }

  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
  }

  plan_to_plot <- eval_tidy(enquo(plan), shp)
  if (!is.null(plan_to_plot)) {
    if (!any(class(plan_to_plot) %in% c('numeric', 'integer', 'character'))) {
      stop('Please provide "plan" as a vector.')
    }
    if (nrow(shp) != length(plan_to_plot)) {
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

  edge_cntr <- edge_center_df(shp, adj)

  nb <- edge_cntr$nb
  centers <- edge_cntr$centers

  # Drop Edges that cross District Boundaries
  if (drop) {
    nb <- nb %>%
      filter(plan_to_plot[i] == plan_to_plot[j])
  }

  # Create Plot
  if (plot_shp) {
      if(!is.null(plan)){
          plot <- ggplot(shp) +
              geom_sf(aes(fill = as.character(plan_to_plot)), size = 0.1) +
              theme_void() +
              theme(legend.position = 'none') +
              geom_sf(data = nb)
      } else{
          plot <- ggplot(shp) +
              geom_sf(size = 0.1) +
              theme_void() +
              geom_sf(data = nb)
      }
  } else {
      plot <- ggplot(nb) +
        geom_sf() +
        theme_void()
    }

  if (centroids) {
    if(!is.null(plan) & !plot_shp){
      plot <- plot + geom_sf(data = centers, aes(color = as.character(plan_to_plot)), size = 2) +
        theme(legend.position = 'none')
    } else {
      plot <- plot + geom_sf(data = centers)
    }
  }

  zoom_to = eval_tidy(enquo(zoom_to), shp)
  if (!is.null(zoom_to)) {
      bbox = sf::st_bbox(sf::st_geometry(shp)[zoom_to])
      plot <- plot + ggplot2::coord_sf(xlim=c(bbox$xmin, bbox$xmax),
                                       ylim=c(bbox$ymin, bbox$ymax))
  }


  plot <- plot + labs(title = title)
  # return plot
  return(plot)
}





edge_center_df <- function(shp, adj){
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
        as.numeric(sf::st_geometry(centers)[[i]]),
        as.numeric(sf::st_geometry(centers)[[j]])
      ),
      nrow = 2,
      byrow = TRUE
    ))))

  suppressWarnings(nb <- sf::st_as_sf(edgedf))
  suppressWarnings(st_crs(nb) <- st_crs(shp))

  return(list(nb = nb, centers = centers))
}



globalVariables(c('start', 'finish'))
