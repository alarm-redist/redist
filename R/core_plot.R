#' Plot Cores
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to color the districts. Required.
#' @param core Required. integer vector produced by \code{redist.identify.cores()}.
#' @param lwd Line width. Defaults to 2.
#'
#' @importFrom dplyr if_else
#'
#' @return ggplot
#' @export
#' @concept plot
redist.plot.cores <- function(shp, plan = NULL, core = NULL, lwd = 2) {
  if (missing(shp)) {
    stop('Please provide an argument to shp.')
  }
  

  plan <- eval_tidy(enquo(plan), shp)
  if (is.null(plan)) {
    if(inherits(shp, 'redist_map')){
     plan <- get_existing(shp) 
    } else {
      stop('Please provide an argument to plan.')
    }
  }
  
  core <- eval_tidy(enquo(core), shp)
  if (missing(core)) {
    stop('Please provide an argument to core.')
  }

  shp$plan <- plan
  shp$core <- core

  shp_un <- shp %>%
      group_by(plan) %>%
      summarize(geometry = st_union(geometry),
                .groups = 'drop') %>%
      suppressMessages()

  shp_cores <- shp %>%
      group_by(plan, core) %>%
      summarize(ct = n(),
                geometry = st_union(geometry),
                .groups = 'drop') %>%
      mutate(ct = if_else(.data$ct == 1, NA_integer_, .data$ct)) %>%
      suppressMessages()

  shp_cores %>%
    ggplot() +
    geom_sf(aes(fill = .data$ct)) +
    ggplot2::scale_fill_distiller(direction = 1, na.value = 'white') +
    geom_sf(fill = NA, data = shp_un, color = 'black', lwd = lwd) +
    labs(fill = 'Number of Units in Core') +
    theme_void() +
    theme(legend.position = 'bottom')
}
