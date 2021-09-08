
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

    exist_col = attr(data, "existing_col")
    if (!is.null(exist_col))  stopifnot(is.numeric(data[[exist_col]]))

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

        if (is.null(attr(data, "pop_bounds")))
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
#' @param existing_plan \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the
#'   existing district assignment. Must be numeric or convertable to numeric.
#' @param pop_tol \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the population tolerance.
#'   The percentage deviation from the average population will be constrained to
#'   be no more than this number. If `existing_plan` is provided, defaults to
#'   the parity of that plan; otherwise, defaults to 0.01.
#' @param total_pop \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the vector
#'   of precinct populations. Defaults to the \code{pop}, \code{population}, or
#'   \code{total_pop} columns, if one exists.
#' @param ndists \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the integer number of
#'   districts to partition the map into. Must be specified if `existing_plan` is not supplied.
#' @param pop_bounds \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} more specific
#'   population bounds, in the form of \code{c(lower, target, upper)}.
#' @param adj the adjacency graph for the object. Defaults to being computed
#'     from the data if it is coercible to a shapefile.
#' @param adj_col the name of the adjacency graph column
#' @param planarize a number, indicating the CRS to project the shapefile to if
#'   it is latitude-longitude based. Set to NULL or FALSE to avoid planarizing.
#'
#' @return A redist_map object
#'
#' @examples
#' data(fl25)
#' d = redist_map(fl25, ndists=3, pop_tol=0.05, total_pop = pop)
#' dplyr::filter(d, pop >= 10e3)
#'
#' @concept prepare
#' @md
#' @export
redist_map = function(..., existing_plan=NULL, pop_tol=NULL,
                      total_pop=c("pop", "population", "total_pop", "POP100"),
                      ndists=NULL, pop_bounds=NULL,
                      adj=NULL, adj_col="adj", planarize=3857) {
    x = tibble(...)
    is_sf = any(vapply(x, function(x) inherits(x, "sfc"), TRUE))
    if (is_sf) {
        x = sf::st_sf(x)

        if (is.na(sf::st_crs(x))) {
            warning("Missing CRS, assuming NAD83 (4269).")
            sf::st_crs(x) = 4269
        }

        if (isTRUE(sf::st_is_longlat(sf::st_geometry(x)))) {
            if (!is.null(planarize) && !isFALSE(planarize)) {
                message("Projecting to CRS ", planarize)
                x = sf::st_transform(x, planarize)
            } else {
                warning("Using latitude and longitude coordinates, ",
                        "which may cause problems with geometric operations")
            }
        }
    }

    pop_col = names(x)[tidyselect::eval_select(rlang::enquo(total_pop), x,
                                               strict=FALSE)]
    if (length(pop_col) == 0) {
        names = rlang::as_label(rlang::enquo(total_pop))
        stop("Population column `", names, "` not found. ",
             "Population must be specified in the `total_pop` argument.")
    } else if (length(pop_col) > 1) {
        pop_col = pop_col[1]
        warning("Multiple potential population columns found, using `", pop_col,
                "`.\nConsider specifying `total_pop` manually.")
    }

    existing_col = names(tidyselect::eval_select(rlang::enquo(existing_plan), x))
    if (length(existing_col) == 0)
      existing_col = NULL

    if (!is.null(existing_col)) {
      if (!is.numeric(x[[existing_col]])) {
        temp_col <- NULL
        suppressWarnings({temp_col <- as.numeric(x[[existing_col]])})
        if (!any(is.na(temp_col))) {
          x[[existing_col]] <- temp_col
        } else {
          stop('`existing_col` was not numeric and could not be converted to numeric.')
        }
      }
    }

    if (is.null(ndists)) {
        if (!is.null(existing_col)) {
            ndists = length(unique(x[[existing_col]]))
        } else {
          stop("Must specify `ndists` if `existing_plan` is not supplied.")
        }
    } else {
        ndists = as.integer(rlang::eval_tidy(rlang::enquo(ndists), x))
    }

    pop_tol = eval_tidy(enquo(pop_tol), x)
    if (is.null(pop_tol) && is.null(pop_bounds)) {
        if (!is.null(existing_col)) {
            pop_tol = redist.parity(x[[existing_col]], x[[pop_col]])
            if (pop_tol <= 0.001)
                message("`pop_tol` calculated from existing plan is \u2264 0.1%")
        } else {
            pop_tol = 0.01
            warning("`pop_tol` not provided; defaulting to 1%")
        }
    }

    if (is.null(pop_bounds)) {
        stopifnot(!is.null(pop_tol))
        stopifnot(pop_tol > 0)

        target = sum(x[[pop_col]]) / ndists
        pop_bounds = target * c(1 - pop_tol, 1, 1 + pop_tol)
    } else {
        pop_bounds = rlang::eval_tidy(rlang::enquo(pop_bounds), x)
    }

    if (is_sf && is.null(adj)) {
        if (!is.null(x[[adj_col]]))
            stop("Column `", adj_col, "` already present in data. ",
                 "Specify an alternate adj column.")

        adj = redist.adjacency(x)
    }

    validate_redist_map(
        new_redist_map(x, adj, ndists, pop_bounds, pop_col, adj_col,
                       add_adj=TRUE, existing_col)
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

#' Extract the target district population from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns a single numeric value, the target population
#' @concept prepare
#' @export
get_target = function(x) {
    stopifnot(inherits(x, "redist_map"))

    attr(x, "pop_bounds")[2]
}

#' Get and set the population tolerance from a \code{redist_map} object
#'
#' @param map the \code{\link{redist_map}} object
#'
#' @return For \code{get_pop_tol}, a single numeric value, the population
#'   tolerance
#'
#' @concept  prepare
#' @export
get_pop_tol <- function(map) {
  stopifnot(inherits(map, 'redist_map'))

  bot <- 1 - attr(map, 'pop_bounds')[1] / attr(map, 'pop_bounds')[2]
  top <- attr(map, 'pop_bounds')[3] / attr(map, 'pop_bounds')[2] - 1

  if (!isTRUE(all.equal(bot, top))) {
    warning('Population bounds were not symmetric, using the smaller tolerance.')
  }

  return(min(bot, top))
}

#' @param pop_tol the population tolerance
#'
#' @rdname get_pop_tol
#' @return For \code{seet_pop_tol}, an updated \code{\link{redist_map}} object
#'
#' @concept  prepare
#' @export
set_pop_tol <- function(map, pop_tol) {
  stopifnot(inherits(map, 'redist_map'))

  target <- get_target(map)
  bot <- (1 - pop_tol) * target
  top <- (1 + pop_tol) * target

  attr(map, 'pop_bounds') <- c(bot, target, top)

  validate_redist_map(map, check_contig = FALSE)
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

    if(bounds[1] > bounds[2] || bounds[3] < bounds[1]){
      warning('Your subset was not based on districts. Please use `set_pop_tol(map)`
              to update your `redist_map` object or create a new `redist_map`
              object with the correct number of districts.')
    }

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

#' @export
#' @importFrom dplyr rename
rename.redist_map <- function(.data, ...) {
    ret = NextMethod()

    cols = tidyselect::eval_rename(rlang::expr(c(...)), .data)
    if (!is.na(colnum <- match(attr(.data, "adj_col"), names(.data)[cols]))) {
        attr(.data, "adj_col") = names(cols)[colnum]
    }
    if (!is.na(colnum <- match(attr(.data, "pop_col"), names(.data)[cols]))) {
        attr(.data, "pop_col") = names(cols)[colnum]
    }
    if (!is.na(colnum <- match(attr(.data, "existing_col"), names(.data)[cols]))) {
        attr(.data, "existing_col") = names(cols)[colnum]
    }

    reconstruct.redist_map(ret, .data)
}

#' @export
#' @importFrom dplyr select
select.redist_map <- function(.data, ...) {
    ret = NextMethod()

    cols = tidyselect::eval_select(rlang::expr(c(...)), .data)
    if (!is.na(colnum <- match(attr(.data, "adj_col"), names(.data)[cols]))) {
        attr(.data, "adj_col") = names(cols)[colnum]
    } else {
        stop("Must keep `", attr(.data, "adj_col"), "` column, ",
             "or convert to a tibble with `as_tibble()`.")
    }
    if (!is.na(colnum <- match(attr(.data, "pop_col"), names(.data)[cols]))) {
        attr(.data, "pop_col") = names(cols)[colnum]
    } else {
        stop("Must keep `", attr(.data, "pop_col"), "` column, ",
             "or convert to a tibble with `as_tibble()`.")
    }
    if (!is.na(colnum <- match(attr(.data, "existing_col"), names(.data)[cols]))) {
        attr(.data, "existing_col") = names(cols)[colnum]
    } else {
        stop("Must keep `", attr(.data, "existing_col"), "` column, ",
             "or convert to a tibble with `as_tibble()`.")
    }

    reconstruct.redist_map(ret, .data)
}


#' Generic to print redist_map
#' @param x redist_map
#' @param \dots additional argumentss
#' @method print redist_map
#' @return Prints to console and returns input redist_map
#' @export
print.redist_map = function(x, ...) {
    cat("A redist_map object with", nrow(x), "units and", ncol(x), "fields\n")

    bounds = attr(x, "pop_bounds")
    cat("To be partitioned into ", attr(x, "ndists"),
        " districts with population between ",
        format(bounds[2], nsmall=0, big.mark=","), " - ",
        format(100 - 100*bounds[1]/bounds[2], nsmall=1), "% and ",
        format(bounds[2], nsmall=0, big.mark=","), " + ",
        format(100*bounds[3]/bounds[2] - 100, nsmall=1), "%\n",
        sep="")

    merge_idx = attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cat("Merged from another map with reindexing:",
            utils::capture.output(str(merge_idx, vec.len=2)), "\n", sep="")

    if (inherits(x, "sf")) {
        geom = st_geometry(x)

        cat("With geometry:\n")
        bb = signif(attr(geom, "bbox"), options("digits")$digits)
        cat("    bbox:           ",
            paste(paste(names(bb), bb[], sep = ": "), collapse = " "),
            "\n", sep="")

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
#' @param fill \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} If
#'   present, will be used to color the map units. If using data masking, may
#'   need to explicitly name argument \code{fill=...} in non-interactive
#'   contexts to avoid S3 generic issues.
#' @param by_distr if \code{TRUE} and \code{fill} is not missing and, color by
#'   district and indicate the \code{fill} variable by shading.
#' @param adj if \code{TRUE}, force plotting the adjacency graph. Overrides
#'   \code{by_distr}.
#' @param interactive if \code{TRUE}, show an interactive map in the viewer
#'   rather than a static map. Ignores \code{adj} and \code{by_distr}.
#' @param ... passed on to \code{\link{redist.plot.map}} (or
#'   \code{\link{redist.plot.adj}} if \code{adj=TRUE}, or
#'   \code{\link{redist.plot.interactive}} if \code{interactive=TRUE}).
#'   Useful parameters may include \code{zoom_to}, \code{boundaries}, and
#'   \code{title}.
#'
#' @examples
#' data(fl25)
#' d = redist_map(fl25, ndists=3, pop_tol=0.05)
#' plot(d)
#' plot(d, BlackPop/pop)
#'
#' data(fl25_enum)
#' fl25$dist <- fl25_enum$plans[, 5118]
#' d <- redist_map(fl25, existing_plan = dist)
#' plot(d)
#'
#'
#' @return ggplot2 object
#' @method plot redist_map
#' @concept prepare
#' @concept plot
#' @export
plot.redist_map = function(x, fill=NULL, by_distr=FALSE, adj=FALSE,
                           interactive=FALSE, ...) {
    if (!inherits(x, "sf"))
        stop("Plotting requires a shapefile.\n  ",
             "If you've just used `merge_by`, consider passing `drop_geom=FALSE`.")

    if (interactive) {
        if (adj) warning("`adj` ignored when `interactive=TRUE`")
        if (by_distr) warning("`by_distr` ignored when `interactive=TRUE`")
        return(redist.plot.interactive(x, !!enquo(fill), ...))
    }

    fill = rlang::enquo(fill)
    existing = get_existing(x)
    if (rlang::quo_is_null(fill)) {
        if (!is.null(existing) && isFALSE(adj)) {
            redist.plot.map(shp = x, adj = get_adj(x), plan=existing, ...)
        } else {
            redist.plot.adj(shp = x, adj = get_adj(x), ...)
        }
    } else {
        fill_name = rlang::quo_text(fill)
        if (!is.null(existing) && isTRUE(by_distr)) {
            redist.plot.map(shp = x, adj = get_adj(x), plan=existing,
                            fill = !!fill, fill_label=fill_name, ...)
        } else {
            redist.plot.map(shp = x, fill = !!fill, fill_label=fill_name, ...)
        }
    }
}
