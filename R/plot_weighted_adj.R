redist.plot.wted.adj <- function(shp = NULL, adj = NULL, plans = NULL,
                                 centroids = TRUE,
                                 drop = FALSE,
                                 plot_shp = TRUE) {

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
        plans <- get_plans_matrix(plans)
    }

    if (inherits(shp, 'redist_map') ){
        if (missing(adj)) {
            adj <- get_adj(shp)
        }
    } else if (missing(adj)){
        adj <- redist.adjacency(shp)
    }

    # Build the sf graph ----
    edge_cntr <- edge_center_df(shp, adj)
    nb <- edge_cntr$nb
    centers <- edge_cntr$centers

    # Create Plot ----
    plot <- ggplot() +
        theme_void()

    if (plot_shp) {
        plot <- plot +
            geom_sf(data = shp, size = 0.1)
    }

    # Add weighted adj


    # return ----
    plot
}
