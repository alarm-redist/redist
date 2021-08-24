#' Display an interactive map
#'
#' Plots an interactive Leaflet map of a \code{\link{redist_map}} object,
#' optionally colored by a quantity of interest. May also be accessed by setting
#' \code{interactive=TRUE} in \code{\link{plot.redist_map}}.
#'
#' If \code{leafgl} is installed, will use its faster rendering functions to
#' plot the map, which may be useful for larger maps.
#'
#' @param map the \code{\link{redist_map}} object
#' @param fill \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} If
#'   present, will be used to color the map units.
#' @param scale the color scale to use, for numeric \code{fill}.
#' @param limits the color scale limits. Defaults to the range of the data.
#' @param useGL if \code{TRUE} and \code{leafgl} is installed, use WebGL for
#'   faster plotting.
#'
#' @returns a Leaflet object
#'
#' @export
redist.plot.interactive = function(map, fill=NULL,
                                   scale=ggplot2::scale_fill_viridis_c,
                                   limits=NULL, useGL=FALSE) {
    if (!requireNamespace("leaflet", quietly=TRUE))
        stop("Package `leaflet` required for interactive plotting.\n",
             "Recommend `leafgl` as well, for faster plotting.")
    stopifnot(inherits(map, "redist_map"))

    fill <- eval_tidy(enquo(fill), map)

    iplot = st_transform(map, 4326) %>%
        st_cast("POLYGON") %>%
        suppressWarnings() %>%
        leaflet::leaflet() %>%
        leaflet::addMapPane(name = "polygons", zIndex = 410) %>%
        leaflet::addMapPane(name = "maplabels", zIndex = 420) %>%
        leaflet::addProviderTiles(leaflet::providers$Stamen.TonerBackground)

    if (useGL && requireNamespace("leafgl", quietly=TRUE)) {
        add_poly = function(map, ...) {
            leafgl::addGlPolygons(map, ..., data=leaflet::getMapData(iplot),
                                  options=leaflet::leafletOptions(pane="polygons")) %>%
                leafgl::addGlPolylines(..., data=st_cast(leaflet::getMapData(iplot), "LINESTRING"),
                                       options=leaflet::leafletOptions(pane="polygons")) %>%
                suppressWarnings()
        }
    } else {
        add_poly = function(map, ...) {
            leaflet::addPolygons(map, ..., data=leaflet::getMapData(iplot),
                                 options=leaflet::leafletOptions(pane="polygons"))
        }
    }

    existing = get_existing(map)
    PAL = c('#6D9537', '#364B7F', '#DCAD35', '#9A9BB9', '#2A4E45', '#7F4E28')
    if (is.null(fill)) {
        if (is.null(existing)) {
            iplot = add_poly(map=iplot, weight=0.5+0.5*!useGL, fillOpacity=0.0,
                             opacity=0.9, color=PAL[2])
        } else {
            ndists = attr(map, "ndists")
            adj = get_adj(map)
            plan = as.integer(as.factor(existing))
            if (ndists > 6)
                plan = color_graph(adj, plan)

            iplot = add_poly(map=iplot, weight=0.5+0.5*!useGL,
                             fillOpacity=0.8, opacity=0.9,
                             fillColor=PAL[plan], color=PAL[plan])
        }
    } else {
        if (is.character(fill) || is.factor(fill)) {
            adj = get_adj(map)
            nlevels = length(unique(fill))
            fill = as.integer(as.factor(fill))
            if (nlevels > 6)
                fill = color_graph(adj, fill)

            iplot = add_poly(map=iplot, weight=0.5+0.5*!useGL,
                             fillOpacity=0.8, opacity=0.9,
                             fillColor=PAL[fill], color=PAL[fill])
        } else {
            pal = scale()$palette
            if (is.null(limits))
                limits = range(fill, na.rm=TRUE)
            else
                limits = dplyr::coalesce(limits, range(fill, na.rm=TRUE))
            fill = (fill - limits[1]) / diff(limits)

            iplot = add_poly(map=iplot, weight=0.5+0.5*!useGL,
                             fillOpacity=0.8, opacity=0.9,
                             fillColor=pal(fill), color=pal(fill))
        }
    }


    leaflet::addProviderTiles(iplot, leaflet::providers$Stamen.TonerLabels,
                              options=leaflet::leafletOptions(pane="maplabels"))
}
