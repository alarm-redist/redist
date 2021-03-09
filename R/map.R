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
#' @param adjacency Deprecated, use adj. A zero-indexed adjacency list. Created with redist.adjacency
#' if not supplied. Default is NULL.
#' @param plan A numeric vector with one entry for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param district_membership Deprecated, use plan. A numeric vector with one row for each precinct in shp.
#' Used to color the districts. Default is \code{NULL}.  Optional.
#' @param centroids A logical indicating if centroids should be plotted. Default is \code{TRUE}.
#' @param edges A logical indicating if edges should connect adjacent centroids. Default is \code{TRUE}.
#' @param boundaries A logical indicating if precinct boundaries should be plotted.
#' @param drop A logical indicating if edges that cross districts should be dropped. Default is \code{FALSE}.
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
#' library(redist)
#' data("fl25")
#' data("algdat.p10")
#' cds <- algdat.p10$cdmat[,100]
#' redist.map(shp = fl25, plan = cds)
#' }
#'
#' @concept plot
#' @export
redist.map <- function(shp = NULL, adj = NULL, adjacency, plan = NULL,
                       district_membership, centroids = TRUE, edges = TRUE,
                       boundaries = TRUE, drop = FALSE, title = "") {
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
        shp <- shp %>%  st_as_sf()
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
            nb <- lapply(adj, function(x) { x + 1L })
        } else{
            adj <- redist.adjacency(shp)
            nb <- lapply(adj, function(x) { x + 1L })
        }

        edgedf <- tibble(start = rep(1:length(nb), lengths(nb)),
                   finish = unlist(nb))
        edgedf <- edgedf %>%
            rowwise() %>%
            mutate(i = min(start, finish), j = max(start, finish)) %>%
            select(i, j)
        edgedf <- edgedf[!duplicated(edgedf), ]

        edgedf <- edgedf %>% rowwise() %>%
            mutate(geometry = st_sfc(st_linestring(matrix(
                c(as.numeric(centers$geometry[[i]]),
                  as.numeric(centers$geometry[[j]])),
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
            plan = as.factor(color_graph(adj, as.integer(plan)))
        }

        plot <- ggplot(shp) +
            geom_sf(aes(fill = plan),
                    size = 0.3 * boundaries,
                    color = "#444444") +
            theme_minimal()     +
            labs(fill = 'District Membership', x = 'Longitude',
                y = 'Latitude', title = title) +
            theme(legend.position = "bottom")

        if (inherits(shp, "redist_map")) {
            PAL = c("#6D9537", "#364B7F", "#DCAD35", "#9A9BB9", "#2A4E45", "#7F4E28")
            plot = plot + ggplot2::guides(fill = FALSE) +
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
#' redist.choropleth(shp = fl25, fill = BlackPop/TotPop)
#' DVS <- fl25$obama/(fl25$mccain+fl25$obama)
#' redist.choropleth(shp = fl25, fill = DVS, limit_colors = c('red', 'blue'))
#' }
#'
#' @concept plot
#' @export
redist.choropleth <- function(shp, fill = NULL, fill_label = "", title = "",
                              limit_colors = NULL, grad = 1, lwd = 0) {

    # Check inputs
    if(missing(shp)) {
        stop('Please provide an argument to "shp".')
    }
    if('SpatialPolygonsDataFrame' %in% class(shp)) {
        shp <- shp %>%  st_as_sf()
    } else if(!('sf' %in% class(shp))) {
        stop('Please provide "shp" as a SpatialPolygonsDataFrame or sf object.')
    }

    fill <- rlang::eval_tidy(rlang::enquo(fill), shp)
    if (!is.numeric(fill) && inherits(shp, "redist_map")) {
        fill = as.factor(color_graph(get_adj(shp), as.integer(as.factor(fill))))
    }

    plot <- ggplot(shp, aes(fill = fill)) +
        geom_sf(size = 0, color = NA) +
        theme_void() +
        labs(fill = fill_label, x = 'Longitude', y = 'Latitude', title = title) +
        theme(legend.position = "bottom")

    if (!is.null(fill)) {
        if (is.numeric(fill)) {
            l1 <- min(fill, na.rm=T)
            l2 <- max(fill, na.rm=T)
            mp <- mean(fill, na.rm=T)

            if (l1 >= 0 & l2 <= 1) {
                l1 <- 0
                l2 <- 1
                mp <- 0.5
            }

            if (is.null(limit_colors)) {
                limit_colors <- c('#ffffff', '#08306b')
            }

            if (grad == 1) {
                plot <- plot +
                    scale_fill_gradient(low = limit_colors[1], high = limit_colors[2],
                                        limits = c(l1, l2))
            } else if (grad == 2) {
                plot <- plot +
                    scale_fill_gradient2(low = limit_colors[1], high = limit_colors[3],
                                         midpoint = mp, mid = limit_colors[2], limits = c(l1, l2))
            }
        } else {
            PAL = c("#6D9537", "#364B7F", "#DCAD35", "#9A9BB9", "#2A4E45", "#7F4E28")
            plot <- plot + ggplot2::scale_fill_manual(values=PAL) +
               ggplot2::guides(fill=FALSE)
        }
    }

    return(plot)
}

globalVariables(c('start', 'finish'))
