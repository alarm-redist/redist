##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################

#######################
# constructors and reconstructors

# Main internal constructor
new_redist_map <- function(
        data, adj, ndists, total_seats, district_seat_sizes,
        pop_bounds, pop_col = "pop",
        adj_col = "adj", add_adj = TRUE, existing_col = NULL) {
    if (add_adj) {
        stopifnot(!is.null(adj))

        data[[adj_col]] <- adj
    }

    if (!is.integer(ndists)) cli_abort("{.arg ndists} must be an integer.")
    if (!is.numeric(pop_bounds) || length(pop_bounds) != 3)
        cli_abort("{.arg pop_bounds} must be a numeric vector of length 3.")
    if(ndists > total_seats){
        cli::cli_abort("{.arg ndists} must be less than or equal to {.arg total_seats}")
    }

    data <- reconstruct.redist_map(data)
    attr(data, "ndists") <- ndists
    attr(data, "total_seats") <- total_seats
    attr(data, "pop_bounds") <- pop_bounds
    attr(data, "district_seat_sizes") <- district_seat_sizes
    attr(data, "pop_col") <- pop_col
    attr(data, "adj_col") <- adj_col
    attr(data, "existing_col") <- existing_col
    # set the districting scheme
    if(assertthat::is.scalar(district_seat_sizes) && district_seat_sizes == 1){
        attr(data, "districting_scheme") <- "SMD"
    }else{
        attr(data, "districting_scheme") <- "MMD"
    }

    data
}

validate_redist_map <- function(data, check_contig = TRUE, call = parent.frame()) {
    if (!inherits(data, "redist_map"))
        cli_abort("Not a {.cls redist_map}", call = call)
    if (!is.data.frame(data))
        cli_abort("Not a data frame", call = call)

    col <- attr(data, "adj_col")
    if (is.null(col)) cli_abort("No adjacency graph column found.", call = call)
    if (!is.list(data[[col]]))
        cli_abort("Adjacency graph column must be a properly formatted adjacency list.", call = call)

    if (check_contig && !is_contiguous(data)) {
        components <- contiguity(get_adj(data), rep(1, nrow(data)))
        disconn <- which(components != which.max(table(components)))
        cli_abort(c("Adjacency graph not contiguous.",
            ">" = "Try manually editing the output of {.fun redist.adjacency}.",
            "i" = "Disconnected precincts: c({paste0(disconn, collapse=', ')})"),
        call = call)
    }

    stopifnot(!is.null(attr(data, "pop_col")))
    stopifnot(!is.null(attr(data, "ndists")))
    if(is.null(attr(data, "total_seats"))){
        attr(data, "total_seats") <- attr(data, "ndists")
    }
    if(is.null(attr(data, "district_seat_sizes"))){
        attr(data, "district_seat_sizes") <- 1L
    }
    stopifnot(!is.null(attr(data, "total_seats")))
    stopifnot(!is.null(attr(data, "district_seat_sizes")))
    district_seat_sizes <- attr(data, "district_seat_sizes")

    if(is.null(attr(data, "districting_scheme"))){
        # set the districting scheme
        if(assertthat::is.scalar(district_seat_sizes) && district_seat_sizes == 1){
            attr(data, "districting_scheme") <- "SMD"
        }else{
            attr(data, "districting_scheme") <- "MMD"
        }
    }
    stopifnot(!is.null(attr(data, "districting_scheme")))


    exist_col <- attr(data, "existing_col")
    if (!is.null(exist_col) && !is.numeric(data[[exist_col]]))
        cli_abort("Existing plan {.field {exist_col}} must be a numeric vector.", call = call)

    pop_bounds <- attr(data, "pop_bounds")
    stopifnot(!is.null(pop_bounds))
    if (!all(diff(pop_bounds) > 0))
        cli_abort("{.arg pop_bounds} must satisfy lower < target < upper.", call = call)

    data
}

reconstruct.redist_map <- function(data, old) {
    classes <- c("tbl_df", "tbl", "data.frame")

    if (inherits(data, "grouped_df"))
        classes <- c("grouped_df", classes)
    if (inherits(data, "rowwise_df"))
        classes <- c("rowwise_df", classes)
    if (inherits(data, "sf"))
        classes <- c("sf", classes)

    if (!missing(old)) {
        if (attr(old, "pop_col") %in% names(data))
            attr(data, "pop_col") <- attr(old, "pop_col")
        if (attr(old, "adj_col") %in% names(data))
            attr(data, "adj_col") <- attr(old, "adj_col")
        if (is.null(attr(data, "merge_idx")))
            attr(data, "merge_idx") <- attr(old, "merge_idx")

        if (isTRUE((exist_col <- attr(old, "existing_col")) %in% names(data))) {
            attr(data, "existing_col") <- exist_col
            attr(data, "ndists") <- length(unique(data[[exist_col]]))
        } else {
            attr(data, "ndists") <- attr(old, "ndists")
        }

        if (is.null(attr(data, "pop_bounds")))
            attr(data, "pop_bounds") <- attr(old, "pop_bounds")

        others <- setdiff(names(attributes(old)), names(attributes(data)))
        if (length(others) > 1) {
            for (i in seq_len(length(others))) {
                attr(data, others[i]) <- attr(old, others[i])
            }
        }
    }

    class(data) <- c("redist_map", classes)
    data
}

#' Create a \code{redist_map} object.
#'
#' Sets up a redistricting problem.
#'
#' A \code{redist_map} object is a \code{\link{tibble}} which contains an
#' adjacency list and additional information about the number of districts, seats, and
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
#' single \code{list} or \code{data.frame}.  These will be passed on to the
#' \code{\link{tibble}} constructor.
#' @param existing_plan \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the
#' existing district assignment. Must be numeric or convertible to numeric.
#' @param pop_tol \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the population tolerance.
#' The percentage deviation from the average population will be constrained to
#' be no more than this number. If `existing_plan` is provided, defaults to
#' the parity of that plan; otherwise, defaults to 0.01.
#' @param total_pop \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the vector
#' of precinct populations. Defaults to the \code{pop}, \code{population}, or
#' \code{total_pop} columns, if one exists.
#' @param ndists \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the integer number of
#' districts to partition the map into. Must be specified if `existing_plan` is not supplied.
#' @param total_seats \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} the total number of
#' legislative seats in the map. For single-member districting schemes this is equal to `ndists`.
#' If not provided then defaults to `ndists`.
#' @param district_seat_sizes Sizes that a valid district is allowed to be. Sizes
#' here refers to the number of seats contained in a district. For single member
#' districting schemes this is always 1.
#' @param pop_bounds \code{\link[dplyr:dplyr_data_masking]{<data-masking>}} more specific
#' population bounds, in the form of \code{c(lower, target, upper)}.
#' @param adj the adjacency graph for the object. Defaults to being computed
#' from the data if it is coercible to a shapefile.
#' @param adj_col the name of the adjacency graph column
#' @param planarize a number, indicating the CRS to project the shapefile to if
#' it is latitude-longitude based. Set to NULL or FALSE to avoid planarizing.
#'
#' @return A redist_map object
#'
#' @examples
#' data(fl25)
#' d <- redist_map(fl25, ndists = 3, pop_tol = 0.05, total_pop = pop)
#' dplyr::filter(d, pop >= 10e3)
#'
#' @concept prepare
#' @md
#' @export
redist_map <- function(..., existing_plan = NULL, pop_tol = NULL,
                       total_pop = c("pop", "population", "total_pop", "POP100"),
                       ndists = NULL, total_seats = NULL, district_seat_sizes = NULL,
                       pop_bounds = NULL,
                       adj = NULL, adj_col = "adj", planarize = 3857) {
    x <- tibble(...)
    is_sf <- any(vapply(x, function(x) inherits(x, "sfc"), TRUE))
    if (is_sf) {
        x <- sf::st_sf(x)

        if (is.na(sf::st_crs(x))) {
            cli_warn("Missing CRS, assuming NAD83 (4269).")
            sf::st_crs(x) <- 4269
        }

        if (isTRUE(sf::st_is_longlat(sf::st_geometry(x)))) {
            if (!is.null(planarize) && !isFALSE(planarize)) {
                message("Projecting to CRS ", planarize)
                x <- sf::st_transform(x, planarize)
            } else {
                cli_warn("Using latitude and longitude coordinates,
                          which may cause problems with geometric operations")
            }
        }
    }

    pop_col <- names(x)[tidyselect::eval_select(rlang::enquo(total_pop), x,
        strict = FALSE)]
    if (length(pop_col) == 0) {
        names <- rlang::as_label(rlang::enquo(total_pop))
        cli_abort(c("Population column {.field {names}} not found.",
            ">" = "Population must be specified in the {.arg total_pop} argument."))
    } else if (length(pop_col) > 1) {
        pop_col <- pop_col[1]
        cli_warn(c("Multiple potential population columns found, using {.field {pop_col}}.",
            ">" = "Consider specifying {.arg total_pop} manually."))
    }
    if (any(is.na(x[[pop_col]]))) {
        cli_abort("The population column {.field {pop_col}} must have no missing values.")
    }

    existing_col <- names(tidyselect::eval_select(rlang::enquo(existing_plan), x))
    if (length(existing_col) == 0)
        existing_col <- NULL

    if (!is.null(existing_col)) {
        if (!is.numeric(x[[existing_col]])) {
            temp_col <- NULL
            suppressWarnings({temp_col <- as.numeric(x[[existing_col]])})
            if (!any(is.na(temp_col))) {
                x[[existing_col]] <- temp_col
            } else {
                cli_abort("Existing plan {.field {exist_col}} must be a numeric vector.")
            }
        }
    }

    if (is.null(ndists)) {
        if (!is.null(existing_col)) {
            ndists <- length(unique(x[[existing_col]]))
        } else {
            cli_abort("Must specify {.arg ndists} if {.arg existing_plan} is not supplied.")
        }
    } else {
        ndists <- as.integer(rlang::eval_tidy(rlang::enquo(ndists), x))
    }

    if(is.null(total_seats)){
        total_seats <- ndists
    }
    # if no district size passed throw error if MMD
    if(is.null(district_seat_sizes)){
        if(ndists == total_seats){
            district_seat_sizes <- 1L
        }else{
            cli::cli_abort("Must specify {.arg district_seat_sizes} if multi-member districting scheme is being used.")
        }
    }else{
        # check the district sizes are numbers
        if(!is.numeric(district_seat_sizes)){
            cli::cli_abort("{.arg district_seat_sizes} must be integers.")
        }
        # check they are integers
        if(!all(as.integer(district_seat_sizes) == district_seat_sizes)){
            cli::cli_abort("{.arg district_seat_sizes} must be integers.")
        }
        # check they are positive
        if(!all(district_seat_sizes > 0)){
            cli::cli_abort("{.arg district_seat_sizes} must be positive.")
        }
        # check its not bigger than the total number of seats
        if(!all(district_seat_sizes <= total_seats)){
            cli::cli_abort("{.arg district_seat_sizes} must all be less than the total number of seats.")
        }
    }

    pop_tol <- eval_tidy(enquo(pop_tol), x)
    if (is.null(pop_tol) && is.null(pop_bounds)) {
        if (!is.null(existing_col)) {
            pop_tol <- redist.parity(x[[existing_col]], x[[pop_col]])
            if (pop_tol <= 0.001)
                cli_inform("{.arg pop_tol} calculated from existing plan is \u2264 0.1%")
        } else {
            pop_tol <- 0.01
            cli_warn("{.arg pop_tol} not provided; defaulting to 1.0%")
        }
    }

    if (is.null(pop_bounds)) {
        stopifnot(!is.null(pop_tol))
        stopifnot(pop_tol > 0)

        target <- sum(x[[pop_col]])/total_seats
        pop_bounds <- target*c(1 - pop_tol, 1, 1 + pop_tol)
    } else {
        pop_bounds <- rlang::eval_tidy(rlang::enquo(pop_bounds), x)
    }

    if (is_sf && is.null(adj)) {
        if (!is.null(x[[adj_col]]))
            cli_abort(c("Column {.field adj_col} already present in data. ",
                ">" = "Specify an alternate adj column."))

        adj <- redist.adjacency(x)
    }

    validate_redist_map(
        new_redist_map(x, adj, ndists, total_seats, district_seat_sizes,
                       pop_bounds, pop_col, adj_col,
            add_adj = TRUE, existing_col)
    )
}

#' @rdname redist_map
#' @param x an object to be coerced
#' @export
as_redist_map <- function(x) {
    reconstruct.redist_map(x)
}

#' Get and set the adjacency graph from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns a zero-indexed adjacency list (\code{get_adj})
#' @concept prepare
#' @export
get_adj <- function(x) {
    if (!inherits(x, "redist_map")) cli_abort("Not a {.cls redist_map}")

    x[[attr(x, "adj_col")]]
}

#' @rdname get_adj
#' @param adj a new adjacency list.
#' @returns the modified \code{redist_map} object (\code{set_adj})
#' @export
set_adj <- function(x, adj) {
    if (!inherits(x, "redist_map")) cli_abort("Not a {.cls redist_map}")
    if (!is.list(adj)) cli_abort("{.arg adj} must be a list")

    # zero-index if need be
    if ((min_idx <- min(sapply(adj, min))) != 0) {
        adj <- lapply(adj, function(x) x - min_idx)
    }

    x[[attr(x, "adj_col")]] <- adj
    # contiguity check etc.
    validate_redist_map(x)
}

#' Extract the existing district assignment from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns an integer vector of district numbers
#' @concept prepare
#' @export
get_existing <- function(x) {
    if (!inherits(x, "redist_map")) cli_abort("Not a {.cls redist_map}")

    exist_col <- attr(x, "existing_col")
    if (is.null(exist_col)) NULL else x[[exist_col]]
}

#' Extract the target district population from a \code{redist_map} object
#'
#' @param x the \code{redist_map} object
#' @returns a single numeric value, the target population
#' @concept prepare
#' @export
get_target <- function(x) {
    if (!inherits(x, "redist_map")) cli_abort("Not a {.cls redist_map}")

    attr(x, "pop_bounds")[2]
}

#' Get and set the population tolerance from a \code{redist_map} object
#'
#' @param map the \code{\link{redist_map}} object
#'
#' @return For \code{get_pop_tol}, a single numeric value, the population
#' tolerance
#'
#' @concept  prepare
#' @export
get_pop_tol <- function(map) {
    if (!inherits(map, "redist_map")) cli_abort("Not a {.cls redist_map}")

    bot <- 1 - attr(map, "pop_bounds")[1]/attr(map, "pop_bounds")[2]
    top <- attr(map, "pop_bounds")[3]/attr(map, "pop_bounds")[2] - 1

    if (!isTRUE(all.equal(bot, top))) {
        cli_warn("Population bounds were not symmetric, using the smaller tolerance.")
    }

    min(bot, top)
}

#' @param pop_tol the population tolerance
#'
#' @rdname get_pop_tol
#' @return For \code{seet_pop_tol}, an updated \code{\link{redist_map}} object
#'
#' @concept  prepare
#' @export
set_pop_tol <- function(map, pop_tol) {
    if (!inherits(map, "redist_map")) cli_abort("Not a {.cls redist_map}")

    target <- get_target(map)
    bot <- (1 - pop_tol)*target
    top <- (1 + pop_tol)*target

    attr(map, "pop_bounds") <- c(bot, target, top)

    validate_redist_map(map, check_contig = FALSE)
}



#######################
# generics

#' @method dplyr_row_slice redist_map
#' @export
dplyr_row_slice.redist_map <- function(data, i, ...) {
    if (is.logical(i)) i <- which(i)

    # reduce adj. graph
    y <- vctrs::vec_slice(data, i)
    gr_col <- attr(data, "adj_col")
    y[[gr_col]] <- redist.reduce.adjacency(data[[gr_col]], keep_rows = i)

    # fix ndists if existing_col exists
    exist_col <- attr(data, "existing_col")
    new_distr <- attr(data, "ndists")
    if (!is.null(exist_col))
        new_distr <- length(unique(y[[exist_col]]))
    attr(y, "ndists") <- new_distr

    # fix merge_idx
    merge_idx <- attr(data, "merge_idx")
    if (!is.null(merge_idx))
        merge_idx <- vctrs::vec_group_id(merge_idx[i])
    attr(y, "merge_idx") <- merge_idx

    # fix pop. bounds
    bounds <- attr(data, "pop_bounds")
    attr(y, "pop_bounds") <- bounds
    new_tgt <- sum(y[[attr(data, "pop_col")]])/new_distr

    if (new_distr > 0) {
        if (bounds[1] > new_tgt || bounds[3] < new_tgt) {
            cli_warn(c("Your subset was not based on districts.",
                       ">" = "Please use {.fn set_pop_tol} to update your
                        {.cls redist_map} or create a new {.cls redist_map}
                        with the correct number of districts."))
        }
    }

    y
}

# 'template' is the old df
#' @method dplyr_reconstruct redist_map
#' @export
dplyr_reconstruct.redist_map <- function(data, template) {
    reconstruct.redist_map(data, template)
}

#' @export
#' @importFrom dplyr summarise
summarise.redist_map <- function(.data, ..., .groups = NULL) {
    ret <- NextMethod()

    # rebuild the graph if need be
    adj_col <- attr(.data, "adj_col")
    if (!(adj_col %in% names(ret))) {
        ret[[adj_col]] <- collapse_adj(get_adj(.data),
            dplyr::group_indices(.data) - 1)
    }

    attr(ret, "merge_idx") <- dplyr::group_indices(.data)

    reconstruct.redist_map(ret, .data)
}

#' @export
#' @importFrom dplyr rename
rename.redist_map <- function(.data, ...) {
    ret <- NextMethod()

    cols <- tidyselect::eval_rename(rlang::expr(c(...)), .data)
    if (!is.na(colnum <- match(attr(.data, "adj_col"), names(.data)[cols]))) {
        attr(.data, "adj_col") <- names(cols)[colnum]
    }
    if (!is.na(colnum <- match(attr(.data, "pop_col"), names(.data)[cols]))) {
        attr(.data, "pop_col") <- names(cols)[colnum]
    }
    if (!is.na(colnum <- match(attr(.data, "existing_col"), names(.data)[cols]))) {
        attr(.data, "existing_col") <- names(cols)[colnum]
    }

    reconstruct.redist_map(ret, .data)
}

#' @export
#' @importFrom dplyr select
select.redist_map <- function(.data, ...) {
    ret <- NextMethod()

    cols <- tidyselect::eval_select(rlang::expr(c(...)), .data)
    if (!is.na(colnum <- match(attr(.data, "adj_col"), names(.data)[cols]))) {
        attr(.data, "adj_col") <- names(cols)[colnum]
    } else {
        cli_abort("Must keep {.field {attr(.data, 'adj_col')} column,
                   or convert to a tibble with {.fun as_tibble}.")
    }
    if (!is.na(colnum <- match(attr(.data, "pop_col"), names(.data)[cols]))) {
        attr(.data, "pop_col") <- names(cols)[colnum]
    } else {
        cli_abort("Must keep {.field {attr(.data, 'pop_col')} column,
                   or convert to a tibble with {.fun as_tibble}.")
    }
    if (!is.null(exist_col <- attr(.data, "existing_col"))) {
        if (!is.na(colnum <- match(exist_col, names(.data)[cols]))) {
            attr(.data, "existing_col") <- names(cols)[colnum]
        } else {
            cli_abort("Must keep {.field {attr(.data, 'existing_col')} column,
                       or convert to a tibble with {.fun as_tibble}.")
        }
    }

    reconstruct.redist_map(ret, .data)
}


#' Generic to print redist_map
#' @param x redist_map
#' @param \dots additional arguments
#' @method print redist_map
#' @return Prints to console and returns input redist_map
#' @export
print.redist_map <- function(x, ...) {
    cli_text("A {.cls redist_map} with {nrow(x)} units and {ncol(x)} fields")

    bounds <- attr(x, "pop_bounds")
    cat("To be partitioned into ", attr(x, "ndists"),
        " districts with population between ",
        format(bounds[2], nsmall = 0, big.mark = ","), " - ",
        format(100 - 100*bounds[1]/bounds[2], nsmall = 1), "% and ",
        format(bounds[2], nsmall = 0, big.mark = ","), " + ",
        format(100*bounds[3]/bounds[2] - 100, nsmall = 1), "%\n",
        sep = "")

    merge_idx <- attr(x, "merge_idx")
    if (!is.null(merge_idx))
        cat("Merged from another map with reindexing:",
            utils::capture.output(str(merge_idx, vec.len = 2)), "\n", sep = "")

    if (inherits(x, "sf")) {
        geom <- st_geometry(x)

        cat("With geometry:\n")
        bb <- signif(attr(geom, "bbox"), options("digits")$digits)
        cat("    bbox:           ",
            paste(paste(names(bb), bb[], sep = ": "), collapse = " "),
            "\n", sep = "")

        crs <- st_crs(geom)
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
#' present, will be used to color the map units. If using data masking, may
#' need to explicitly name argument \code{fill=...} in non-interactive
#' contexts to avoid S3 generic issues.
#' @param by_distr if \code{TRUE} and \code{fill} is not missing and, color by
#' district and indicate the \code{fill} variable by shading.
#' @param adj if \code{TRUE}, force plotting the adjacency graph. Overrides
#' \code{by_distr}.
#' @param ... passed on to \code{\link{redist.plot.map}} (or
#' \code{\link{redist.plot.adj}} if \code{adj=TRUE}).
#' Useful parameters may include \code{zoom_to}, \code{boundaries}, and
#' \code{title}.
#'
#' @examples
#' data(fl25)
#' d <- redist_map(fl25, ndists = 3, pop_tol = 0.05)
#' plot(d)
#' plot(d, BlackPop/pop)
#'
#' data(fl25_enum)
#' fl25$dist <- fl25_enum$plans[, 5118]
#' d <- redist_map(fl25, existing_plan = dist)
#' plot(d)
#'
#' @return ggplot2 object
#' @method plot redist_map
#' @concept prepare
#' @concept plot
#' @export
plot.redist_map <- function(x, fill = NULL, by_distr = FALSE, adj = FALSE, ...) {
    if (!inherits(x, "sf"))
        cli_abort(c("Plotting requires a shapefile.",
            ">" = "If you've just used {.fun merge_by},
                    consider passing {.arg drop_geom = FALSE}."))

    fill <- rlang::enquo(fill)
    existing <- get_existing(x)
    if (rlang::quo_is_null(fill)) {
        if (!is.null(existing) && isFALSE(adj)) {
            redist.plot.map(shp = x, adj = get_adj(x), plan = existing, ...)
        } else {
            redist.plot.adj(shp = x, adj = get_adj(x), ...)
        }
    } else {
        fill_name <- rlang::quo_text(fill)
        if (!is.null(existing) && isTRUE(by_distr)) {
            redist.plot.map(shp = x, adj = get_adj(x), plan = existing,
                fill = !!fill, fill_label = fill_name, ...)
        } else {
            redist.plot.map(shp = x, fill = !!fill, fill_label = fill_name, ...)
        }
    }
}
