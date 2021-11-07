#' Plot Weighted Border Adjacency
#'
#' Plots the weighted adjacency graph by how often precincts coocur. If an argument
#' to counties is provided, it subsets the edges to plot to those that cross over the county boundary.
#'
#' @param shp A SpatialPolygonsDataFrame, sf object, or redist_map. Required.
#' @param plans A `redist_plans` object or matrix of redistricting plans, where each
#' column indicates a plan and each
#' @param counties unquoted name of a column in `shp` or a vector of county assignments. Subsets to
#' edges which cross this boundary if supplied.
#' @param ref Plot reference map? Defaults to TRUE which gets the existing plan from
#' @param adj A zero-indexed adjacency list. Extracted from `shp` if `shp` is a `redist_map`.
#' Otherwise created with redist.adjacency if not supplied. Default is NULL.
#' @param plot_shp Should the shapes be plotted? Default is TRUE.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' data(iowa)
#' shp <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
#' plans <- redist_smc(shp, 100)
#' redist.plot.wted.adj(shp, plans = plans, counties = region)
redist.plot.wted.adj <- function(shp = NULL, plans = NULL, counties = NULL,
                                 ref = TRUE, adj = NULL, plot_shp = TRUE) {

  # Check inputs ----
  if (is.null(shp)) {
    stop('Please provide an argument to `shp`.')
  }

  if ('SpatialPolygonsDataFrame' %in% class(shp)) {
    shp <- shp %>% st_as_sf()
  } else if (!('sf' %in% class(shp))) {
    stop('Please provide `shp` as a SpatialPolygonsDataFrame or sf object.')
  }

  if (is.null(plans)) {
    stop('`plans` is required.')
  } else if (inherits(plans, 'redist_plans')) {
    plans <- plans %>%
      subset_sampled() %>%
      get_plans_matrix()
  }

  if (inherits(shp, 'redist_map')) {
    if (missing(adj)) {
      adj <- get_adj(shp)
    }
  } else if (missing(adj)) {
    adj <- redist.adjacency(shp)
  }

  counties <- rlang::eval_tidy(rlang::enquo(counties), shp)

  # Build the sf graph ----
  edge_cntr <- edge_center_df(shp, adj)
  nb <- edge_cntr$nb
  # centers <- edge_cntr$centers

  # Create Plot ----
  p <- ggplot() +
    theme_void()

  if (plot_shp) {
    p <- p +
      geom_sf(data = shp, size = 0.1, fill = NA)


    if (!is.null(counties)) {
      cty <- shp %>%
        mutate(counties_input = counties) %>%
        group_by(counties_input) %>%
        summarize(geometry = st_union(geometry))

      p <- p +
        geom_sf(data = cty, size = 0.4, fill = NA)

      nb <- nb %>%
        filter(counties[i] != counties[j])
    }


    if (is.logical(ref)) {
      if (ref) {
        ref <- get_existing(shp)

        distr <- shp %>%
          mutate(ref_input = ref) %>%
          group_by(ref_input) %>%
          summarize(geometry = st_union(geometry))
        p <- p +
          geom_sf(data = distr, size = 1, fill = NA)
      }
    } else {
      if (length(ref == nrow(shp))) {
        distr <- shp %>%
          mutate(ref_input = ref) %>%
          group_by(ref_input) %>%
          summarize(geometry = st_union(geometry))
        p <- p +
          geom_sf(data = distr, size = 1, fill = NA)
      }
    }
  }

  # Add weighted adj ----
  cooc <- prec_cooccur(plans, seq_len(ncol(plans)))
  nb <- nb %>%
    mutate(wt = cooc[i, j])

  p <- p +
    geom_sf(data = nb, aes(color = wt), lwd = 1) +
    ggplot2::scale_color_distiller(palette = 'Reds', direction = 1) +
    labs(color = 'Coocurrence\nProportion')

  # return ----
  p
}
