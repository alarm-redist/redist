##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

#######################
# constructors and reconstructors

# Main internal constructor
new_redist_map = function(data, graph, n_distr, pop_bounds, pop_col="pop",
                          graph_col="graph", add_graph=TRUE, existing_col=NULL) {
    if (add_graph) {
        stopifnot(!is.null(graph))
        if (!is.null(data[[graph_col]]))
            stop("Column `", graph_col, "` already present in data. Specify an alternate graph column.")

        data[[graph_col]] = graph
    }

    stopifnot(is.integer(n_distr))
    stopifnot(is.numeric(pop_bounds))
    stopifnot(length(pop_bounds) == 3)

    data = reconstruct.redist_map(data)
    attr(data, "n_distr") = n_distr
    attr(data, "pop_bounds") = pop_bounds
    attr(data, "pop_col") = pop_col
    attr(data, "graph_col") = graph_col
    attr(data, "existing_col") = existing_col

    data
}

validate_redist_map = function(data, check_contig=T) {
    if (!is.data.frame(data)) stop("Not a data frame")
    if (!inherits(data, "redist_map")) stop("Not a `redist_map` object")

    col = attr(data, "graph_col")
    if (is.null(col)) stop("No graph column found")
    if (!is.list(data[[col]]))
        stop("Graph column not a properly formatted adjacency list.")

    if (check_contig && !is_contiguous(data))
        stop("Graph not contiguous.")

    stopifnot(!is.null(attr(data, "pop_col")))
    stopifnot(!is.null(attr(data, "n_distr")))

    pop_bounds = attr(data, "pop_bounds")
    stopifnot(!is.null(pop_bounds))
    if (!all(diff(pop_bounds) > 0))
        stop("`pop_bounds` must satisfy lower < target < upper.")

    data
}

reconstruct.redist_map = function(data, old) {
    classes = c("tbl_df", "tbl", "data.frame")

    if (inherits(data, "grouped_df"))
        classes = c("grouped_df", classes)
    if (inherits(data, "sf"))
        classes = c("sf", classes)

    if (!missing(old)) {
        if (attr(old, "pop_col") %in% colnames(data))
            attr(data, "pop_col") = attr(old, "pop_col")
        if (attr(old, "graph_col") %in% colnames(data))
            attr(data, "graph_col") = attr(old, "graph_col")
        if (isTRUE((exist_col <- attr(old, "existing_col")) %in% colnames(data))) {
            attr(data, "existing_col") = exist_col
            attr(data, "n_distr") = length(unique(data[[exist_col]]))
        } else {
            attr(data, "n_distr") = attr(old, "n_distr")
        }

        attr(data, "pop_bounds") = attr(old, "pop_bounds")
    }

    class(data) = c("redist_map", classes)
    data
}

#' Create a \code{redist_map} object.
#'
#' Sets up a redistricting problem.
#'
#' A \code{redist_map} object is a \code{\link{tibble}} which contains an
#' adjacency list and additional information about the number of districts and
#' population bounds.  It supports all of the \code{dplyr} generics, and will
#' adjust the adjacency list and attributes according to these functions; i.e.,
#' if we \code{filter} to a subset of units, the graph will change to subset to
#' these units, and the population bounds will adjust accordingly.  If an
#' existing map is also attached to the object, the number of districts will
#' also adjust.  Subsetting with \code{`[`} and \code{`[[`} does not recompute
#' graphs or attributes.
#'
#' @param ... column elements to be bound into a \code{redist_map} object or a
#'   single \code{list} or \code{data.frame}.  These will be passed on to the
#'   \code{\link{tibble}} constructor.
#' @param n_distr \code{\link[dplyr]{<data-masking>}} the integer number of
#'   districts to partition the map into
#' @param pop_tol \code{\link[dplyr]{<data-masking>}} the population tolerance.
#'   The percentage deviation from the average population will be constrained to
#'   be no more than this number.
#' @param pop_bounds \code{\link[dplyr]{<data-masking>t} more specific
#'   population bounds, in the form of \code{c(lower, target, upper)}.
#' @param pop_col \code{\link[dplyr]{<tidy-select>}} the name of the population
#'   vector column
#' @param graph the adjacency graph for the object. Defaults to being computed
#'     from the data if it is coercible to a shapefile.
#' @param graph_col the name of the adjacency graph column
#' @param existing_col \code{\link[dplyr]{<tidy-select>}} the name of a column
#'   with existing district assignment
#' @param planarize a number, indicating the CRS to project the shapefile to if
#'   it is latitude-longitude based. Set to NULL to avoid planarizing.
#'
#' @examples
#' d = redist_map(fl25, n_distr=3, pop_tol=0.05)
#' dplyr::filter(d, pop >= 10e3)
#'
#' @export
redist_map = function(..., n_distr, pop_tol=0.01, pop_bounds=NULL, pop_col="pop",
                      graph=NULL, graph_col="graph", existing_col=NULL, planarize=3857) {
    x = tibble(...)
    is_sf = any(vapply(x, function(x) inherits(x, "sfc"), TRUE))
    if (is_sf) {
        x = sf::st_sf(x)

        if (is.na(sf::st_crs(x))) {
            warning("Missing CRS, assuming NAD83 (4269).")
            sf::st_crs(x) = 4269
        }

        if (isTRUE(sf::st_is_longlat(sf::st_geometry(x)))) {
            if (!is.null(planarize)) {
                message("Projecting to CRS ", planarize)
                x = sf::st_transform(x, planarize)
            } else {
                warning("Using latitude and longitude coordinates, ",
                        "which may cause problems with geometric operations")
            }
        }
    }

    if (is_sf && is.null(graph))
        graph = redist.adjacency(x)

    n_distr = as.integer(rlang::eval_tidy(rlang::enquo(n_distr), x))
    pop_tol = rlang::eval_tidy(rlang::enquo(pop_tol), x)

    if (is.null(pop_bounds)) {
        stopifnot(!is.null(pop_tol))

        target = sum(x[[pop_col]]) / n_distr
        pop_bounds = target * c(1 - pop_tol, 1, 1 + pop_tol)
    } else {
        pop_bounds = rlang::eval_tidy(rlang::enquo(pop_bounds), x)
    }

    pop_col = names(tidyselect::eval_select(rlang::enquo(pop_col), x))
    existing_col = names(tidyselect::eval_select(rlang::enquo(existing_col), x))
    if (length(existing_col) == 0)
        existing_col = NULL

    validate_redist_map(
        new_redist_map(x, graph, n_distr, pop_bounds, pop_col, graph_col,
                   add_graph=T, existing_col)
    )
}

#' @rdname redist_map
#' @export
as_redist_map = function(x) {
    reconstruct.redist_map(x)
}

# extract graph
get_graph = function(x) {
    stopifnot(inherits(x, "redist_map"))

    x[[attr(x, "graph_col")]]
}

# extract graph
get_existing = function(x) {
    stopifnot(inherits(x, "redist_map"))

    exist_col = attr(x, "existing_col")
    if (is.null(exist_col)) NULL else x[[exist_col]]
}


#######################
# generics

#' @method dplyr_row_slice redist_map
#' @export
dplyr_row_slice.redist_map = function(data, i, ...) {
    if (is.logical(i)) i = which(i)

    # reduce adj. graph
    y = vctrs::vec_slice(data, i)
    gr_col = attr(data, "graph_col")
    y[[gr_col]] = redist.reduce.adjacency(data[[gr_col]], i)

    # fix n_distr if existing_col exists
    exist_col = attr(data, "existing_col")
    new_distr = attr(data, "n_distr")
    if (!is.null(exist_col))
        new_distr = length(unique(y[[exist_col]]))
    attr(y, "n_distr") = new_distr

    # fix pop. bounds
    bounds = attr(data, "pop_bounds")
    bounds[2] = sum(y[[attr(data, "pop_col")]]) / new_distr
    attr(y, "pop_bounds") = bounds

    y
}

# 'template' is the old df
#' @method dplyr_reconstruct redist_map
#' @export
dplyr_reconstruct.redist_map = function(data, template) {
    reconstruct.redist_map(data, template)
}

#' @export
#' @importFrom dplyr summarise
summarise.redist_map = function(.data, ...) {
    ret = NextMethod()

    # rebuild the graph if need be
    graph_col = attr(.data, "graph_col")
    if (!(graph_col %in% colnames(ret))) {
        ret[[graph_col]] = collapse_adj(get_graph(.data),
                                        dplyr::group_indices(.data) - 1)
    }

    attr(ret, "merge_idx") = dplyr::group_indices(.data)

    reconstruct.redist_map(ret, .data)
}


#' @method print redist_map
#' @export
print.redist_map = function(x, ...) {
    cli::cat_line("A redist_map object with ", nrow(x),
                  " units and ", ncol(x), " fields")

    bounds = attr(x, "pop_bounds")
    cli::cat_line("To be partitioned into ", attr(x, "n_distr"),
                  " districts with population between ",
                  format(bounds[2], nsmall=0, big.mark=","), " - ",
                  format(100 - 100*bounds[1]/bounds[2], nsmall=1), "% and ",
                  format(bounds[2], nsmall=0, big.mark=","), " + ",
                  format(100*bounds[3]/bounds[2] - 100, nsmall=1), "%")

    merge_idx = attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cli::cat_line("Merged from another map with reindexing: ", capture_output(str(merge_idx)))

    if (inherits(x, "sf")) {
        geom = st_geometry(x)

        cli::cat_line("With geometry:")
        bb = signif(attr(geom, "bbox"), options("digits")$digits)
        cli::cat_line("    bbox:           ",
                      paste(paste(names(bb), bb[], sep = ": "), collapse = " "))

        crs = st_crs(geom)
        if (is.na(crs)) {
            cat(paste0("    CRS:            NA\n"))
        } else {
            p = sf:::crs_parameters(crs)
            if (p$Name == "unknown") {
                if (!is.character(crs$input) || is.na(crs$input))
                    cat(paste0("proj4string:    ", crs$proj4string,
                               "\n"))
                else cat(paste0("CRS:            ", crs$input, "\n"))
            }
            else if (p$IsGeographic)
                cat(paste0("    geographic CRS: ", p$Name, "\n"))
            else
                cat(paste0("    projected CRS:  ", p$Name, "\n"))
        }
    }

    getS3method("print", "tbl")(x)

    invisible(x)
}

#' Plot a \code{redist_map}
#'
#' @param x the \code{redist_map} object
#' @param y \code{\link[dplyr]{<data-masking>}} the optional value to use to
#'   color the units. If absent, \code{\link{redist.map}} will be called;
#'   otherwise \code{\link{redist.choropleth}} will be called.
#' @param ... passed on to underlying functions
#'
#' @examples
#' d = redist_map(fl25, n_distr=3, pop_tol=0.05)
#' plot(d)
#' plot(d, adjacency=F)
#' plot(d, BlackPop/pop)
#'
#' @method plot redist_map
#' @export
plot.redist_map = function(x, y, ...) {
    if (!inherits(x, "sf")) stop("Plotting requires a shapefile.")

    if (missing(y)) {
        orig_col = attr(x, "existing_col")
        if (!is.null(orig_col)) {
            redist.map(x, get_graph(x), x[[orig_col]], ...) +
                ggplot2::theme_void()
        } else {
            redist.map(x, get_graph(x), ...) +
                ggplot2::theme_void()
        }
    } else {
        redist.choropleth(x, !!rlang::enquo(y), ...)
    }
}

