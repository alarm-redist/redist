#' Plot Cores
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to color the districts. Required.
#' @param core Required. integer vector produced by \code{redist.identify.cores()}.
#'
#' @return ggplot
#' @export
#' @concept plot
redist.plot.cores <- function(shp, plan, core) {
  if (missing(shp)) {
    stop('Please provide an argument to shp.')
  }
  if (missing(plan)) {
    stop('Please provide an argument to plan.')
  }
  if (missing(core)) {
    stop('Please provide an argument to core.')
  }

  shp$plan <- plan
  shp$core <- core

  suppressMessages(
  shp_un <- shp %>%
    group_by(plan) %>%
    summarize(geometry = st_union(geometry),
              .groups = 'drop'))

  suppressMessages(
  shp_cores <- shp %>%
    group_by(plan, core) %>%
    summarize(
      ct = n(),
      geometry = st_union(geometry),
      .groups = 'drop'
    ) %>%
    mutate(ct = ifelse(ct == 1, NA_integer_, ct)))

  shp_cores %>%
    ggplot() +
    geom_sf(aes(fill = ct)) +
    scale_fill_distiller(direction = 1, na.value = 'white') +
    geom_sf(fill = NA, data = shp_un, color = 'black', lwd = 2) +
    labs(fill = 'Number of Units in Core') +
    theme_void() + 
    theme(legend.position = 'bottom')
}
