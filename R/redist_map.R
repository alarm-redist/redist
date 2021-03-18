
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

#######################
# constructors and reconstructors

# Main internal constructor
new_redist_map = function(data, adj, ndists, pop_bounds, pop_col="pop",
                          adj_col="adj", add_adj=TRUE, existing_col=NULL) {
    if (add_adj) {
        stopifnot(!is.null(adj))

        data[[adj_col]] = adj
    }

    stopifnot(is.integer(ndists))
    stopifnot(is.numeric(pop_bounds))
    stopifnot(length(pop_bounds) == 3)

    data = reconstruct.redist_map(data)
    attr(data, "ndists") = ndists
    attr(data, "pop_bounds") = pop_bounds
    attr(data, "pop_col") = pop_col
    attr(data, "adj_col") = adj_col
    attr(data, "existing_col") = existing_col

    data
}

validate_redist_map = function(data, check_contig=T) {
    if (!is.data.frame(data)) stop("Not a data frame")
    if (!inherits(data, "redist_map")) stop("Not a `redist_map` object")

    col = attr(data, "adj_col")
    if (is.null(col)) stop("No adjacency graph column found")
    if (!is.list(data[[col]]))
        stop("Adjacency graph column not a properly formatted adjacency list.")

    if (check_contig && !is_contiguous(data)) {
        components = contiguity(get_adj(data), rep(1, nrow(data)))
        disconn = which(components != which.max(table(components)))
        stop("Adjacency graph not contiguous.\n",
             "Try manually editing the output of `redist.adjacency`.\n",
             "Disconnected precincts: c(", paste0(disconn, collapse=", "), ")")
    }

    stopifnot(!is.null(attr(data, "pop_col")))
    stopifnot(!is.null(attr(data, "ndists")))

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
        if (attr(old, "adj_col") %in% colnames(data))
            attr(data, "adj_col") = attr(old, "adj_col")
        if (is.null(attr(data, "merge_idx")))
            attr(data, "merge_idx") = attr(old, "merge_idx")

        if (isTRUE((exist_col <- attr(old, "existing_col")) %in% colnames(data))) {
            attr(data, "existing_col") = exist_col
            attr(data, "ndists") = length(unique(data[[exist_col]]))
        } else {
            attr(data, "ndists") = attr(old, "ndists")
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
#' Other useful methods for \code{redist_map} objects:
#' * \code{\link{merge_by}}
#' * \code{\link{get_adj}}
#' * \code{\link{plot.redist_map}}
#'
#' @param ... column elements to be bound into a \code{redist_map} object or a
#'   single \code{list} or \code{data.frame}.  These will be passed on to the
#'   \code{\link{tibble}} constructor.
#' @param ndists \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the integer number of
#'   districts to partition the map into
#' @param pop_tol \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the population tolerance.
#'   The percentage deviation from the average population will be constrained to
#'   be no more than this number.
#' @param pop_bounds \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} more specific
#'   population bounds, in the form of \code{c(lower, target, upper)}.
#' @param total_pop \code{\link[tidyr:tidyr_tidy_select]{<tidy-select>}} the vector
#'   of precinct populations. Defaults to the \code{pop} column, if one exists.
#' @param adj the adjacency graph for the object. Defaults to being computed
#'     from the data if it is coercible to a shapefile.
#' @param adj_col the name of the adjacency graph column
#' @param existing_plan \code{\link[tidyr:tidyr_tidy_select]{<tidy-select>}} the
#'   existing district assignment.
#' @param planarize a number, indicating the CRS to project the shapefile to if
#'   it is latitude-longitude based. Set to NULL to avoid planarizing.
#'
#' @examples
#' data(fl25)
#' d = redist_map(fl25, ndists=3, pop_tol=0.05)
#' dplyr::filter(d, pop >= 10e3)
#'
#' @concept prepare
#' @md
#' @export
redist_map = function(..., ndists=NULL, pop_tol=0.01, pop_bounds=NULL, total_pop="pop",
                      adj=NULL, adj_col="adj", existing_plan=NULL, planarize=3857) {
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

    if (is_sf && is.null(adj)) {
        if (!is.null(x[[adj_col]]))
            stop("Column `", adj_col, "` already present in data. ",
                 "Specify an alternate adj column.")

        adj = redist.adjacency(x)
    }

    pop_col = names(tidyselect::eval_select(rlang::enquo(total_pop), x))
    existing_col = names(tidyselect::eval_select(rlang::enquo(existing_plan), x))
    if (length(existing_col) == 0)
        existing_col = NULL

    if (is.null(ndists))  {
        if (!is.null(existing_col))
            ndists = length(unique(x[[existing_col]]))
        else
            stop("Must specify `ndists` if `existing_plan` is not supplied")
    } else {
        ndists = as.integer(rlang::eval_tidy(rlang::enquo(ndists), x))
    }

    pop_tol = rlang::eval_tidy(rlang::enquo(pop_tol), x)

    if (is.null(pop_bounds)) {
        stopifnot(!is.null(pop_tol))
        stopifnot(pop_tol > 0)

        target = sum(x[[pop_col]]) / ndists
        pop_bounds = target * c(1 - pop_tol, 1, 1 + pop_tol)
    } else {
        pop_bounds = rlang::eval_tidy(rlang::enquo(pop_bounds), x)
    }


    validate_redist_map(
        new_redist_map(x, adj, ndists, pop_bounds, pop_col, adj_col,
                   add_adj=T, existing_col)
    )
}

#' @rdname redist_map
#' @param x an object to be coerced
#' @export
as_redist_map = function(x) {
    reconstruct.redist_map(x)
}

#' Get and set the adjacency graph from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns a zero-indexed adjacency list (\code{get_adj})
#' @concept prepare
#' @export
get_adj = function(x) {
    stopifnot(inherits(x, "redist_map"))

    x[[attr(x, "adj_col")]]
}

#' @rdname get_adj
#' @param adj a new adjacency list.
#' @returns the modified \code{redist_map} object (\code{set_adj})
#' @export
set_adj = function(x, adj) {
    stopifnot(inherits(x, "redist_map"))
    stopifnot(is.list(adj))
    # zero-index if need be
    if ((min_idx <- min(sapply(adj, min))) != 0) {
        adj = lapply(adj, function(x) x - min_idx)
    }

    x[[attr(x, "adj_col")]] = adj
    # contiguity check etc.
    validate_redist_map(x)
}

#' Extract the existing district assignment from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns an integer vector of district numbers
#' @concept prepare
#' @export
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
    gr_col = attr(data, "adj_col")
    y[[gr_col]] = redist.reduce.adjacency(data[[gr_col]], keep_rows=i)

    # fix ndists if existing_col exists
    exist_col = attr(data, "existing_col")
    new_distr = attr(data, "ndists")
    if (!is.null(exist_col))
        new_distr = length(unique(y[[exist_col]]))
    attr(y, "ndists") = new_distr

    # fix merge_idx
    merge_idx = attr(data, "merge_idx")
    if (!is.null(merge_idx))
        merge_idx = as.integer(as.factor(merge_idx[i]))
    attr(y, "merge_idx") = merge_idx

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
summarise.redist_map = function(.data, ..., .groups=NULL) {
    ret = NextMethod()

    # rebuild the graph if need be
    adj_col = attr(.data, "adj_col")
    if (!(adj_col %in% colnames(ret))) {
        ret[[adj_col]] = collapse_adj(get_adj(.data),
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
    cli::cat_line("To be partitioned into ", attr(x, "ndists"),
                  " districts with population between ",
                  format(bounds[2], nsmall=0, big.mark=","), " - ",
                  format(100 - 100*bounds[1]/bounds[2], nsmall=1), "% and ",
                  format(bounds[2], nsmall=0, big.mark=","), " + ",
                  format(100*bounds[3]/bounds[2] - 100, nsmall=1), "%")

    merge_idx = attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cli::cat_line("Merged from another map with reindexing:",
                      utils::capture.output(str(merge_idx, vec.len=2)))

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
            if (crs$Name == "unknown") {
                if (!is.character(crs$input) || is.na(crs$input))
                    cat(paste0("proj4string:    ", crs$proj4string,
                               "\n"))
                else cat(paste0("CRS:            ", crs$input, "\n"))
            }
            else if (crs$IsGeographic)
                cat(paste0("    geographic CRS: ", crs$Name, "\n"))
            else
                cat(paste0("    projected CRS:  ", crs$Name, "\n"))
        }
    }

    utils::getS3method("print", "tbl")(x)

    invisible(x)
}

#' Plot a \code{redist_map}
#'
#' @param x the \code{redist_map} object
#' @param val \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the optional value to use to
#'   color the units. If absent, \code{\link{redist.map}} will be called;
#'   otherwise \code{\link{redist.choropleth}} will be called.
#' @param ... passed on to underlying functions
#'
#' @examples
#' data(fl25)
#' d = redist_map(fl25, ndists=3, pop_tol=0.05)
#' plot(d)
#' plot(d, edges=FALSE)
#' plot(d, BlackPop/pop)
#'
#' @method plot redist_map
#' @concept prepare
#' @concept plot
#' @export
plot.redist_map = function(x, val=NULL, ...) {
    if (!inherits(x, "sf")) stop("Plotting requires a shapefile.")

    val = rlang::enquo(val)
    if (rlang::quo_is_null(val)) {
        existing = get_existing(x)
        if (!is.null(existing)) {
            redist.map(x, get_adj(x), plan=existing, ...) +
                ggplot2::theme_void()
        } else {
            redist.map(x, get_adj(x), plan=NULL, ...) +
                ggplot2::theme_void()
        }
    } else {
        redist.choropleth(shp = x, fill = !!val, ...)
    }
}

