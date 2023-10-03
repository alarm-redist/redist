

#' Check that a \code{redist_map} object is contiguous
#'
#' @param x the object
#'
#' @return \code{TRUE} if contiguous.
#' @concept prepare
#'
#' @export
is_contiguous <- function(x) {
    if (!inherits(x, "redist_map")) cli_abort("{.arg x} must be a {.cls redist_map}")
    all(contiguity(get_adj(x), rep(1, nrow(x))) == 1)
}


## merge helpers
# checks if is a proportion/pct
is_prop <- function(x) is.double(x) && min(x, na.rm = TRUE) >= 0 && max(x, na.rm = TRUE) <= 1
# checks if is not a  proportion/pct
is_nonprop <- function(x) is.numeric(x) && (min(x, na.rm = TRUE) < 0 || max(x, na.rm = TRUE) > 1)
# checks if x is constant within levels of rel
is_const_rel <- function(rel) {
    function(x) {
        !is.numeric(x) && all(tapply(x, rel, FUN = function(y) length(unique(y))) == 1)
    }
}

#' Merge map units
#'
#' In performing a county-level or cores-based analysis it is often necessary to
#' merge several units together into a larger unit.  This function performs this
#' operation, modifying the adjacency graph as needed and attempting to properly
#' aggregate other data columns.
#'
#' @param .data a \code{\link{redist_map}} object
#' @param ... \code{\link[dplyr:dplyr_tidy_select]{<tidy-select>}} the column(s) to merge by
#' @param by_existing if an existing assignment is present, whether to also group by it
#' @param drop_geom whether to drop the geometry column. Recommended, as
#' otherwise a costly geometric merge is required.
#' @param collapse_chr if \code{TRUE}, preserve character columns by collapsing
#' their values. For example, a county name column in Iowa might be merged and
#' have entries such as "Cedar~Clinton~Des Moines". Set to \code{FALSE} to
#' drop character columns instead.
#'
#' @returns A merged \code{\link{redist_map}} object
#'
#' @concept prepare
#' @export
merge_by <- function(.data, ..., by_existing = TRUE, drop_geom = TRUE, collapse_chr = TRUE) {
    .data <- as_redist_map(.data)

    dots <- rlang::enquos(...)
    key_val <- rlang::eval_tidy(dots[[1]], .data)
    pop_col <- attr(.data, "pop_col")

    if (drop_geom && inherits(.data, "sf"))
        .data <- sf::st_drop_geometry(.data)

    col <- attr(.data, "existing_col")
    unique_chr <- function(x) paste(unique(x), collapse = "~")
    is_col_chr <- if (collapse_chr) is.character else (function(x) FALSE)
    if (!is.null(col) && by_existing) {
        dplyr::group_by(.data, dplyr::across(dplyr::all_of(col)), !!!dots) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w = .data[[pop_col]], na.rm = TRUE)),
                             dplyr::across(where(is.numeric), sum, na.rm = TRUE),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is_col_chr), unique_chr),
                             .groups = "drop")
    } else {
        dplyr::group_by(.data, !!!dots) %>%
            dplyr::summarize(dplyr::across(where(is_prop),
                                           ~ weighted.mean(., w = .data[[pop_col]], na.rm = TRUE)),
                             dplyr::across(where(is.numeric), sum, na.rm = TRUE),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is_col_chr), unique_chr),
                             .groups = "drop") %>%
            `attr<-`("existing_col", NULL)
    }
}

#' @rdname redist.identify.cores
#' @order 1
#'
#' @param .data a \code{\link{redist_map}} object
#' @param boundary the core is defined to be at least this number of steps within
#' district boundaries
#'
#' @concept prepare
#' @export
make_cores <- function(.data = cur_map(), boundary = 1, focus = NULL) {
    if (is.null(.data))
        cli_abort("Must provide {.arg .data} if not called within a {.pkg dplyr} verb")
    if (!inherits(.data, "redist_map")) cli_abort("{.arg .data} must be a {.cls redist_map}")

    existing <- get_existing(.data)
    if (is.null(existing))
        cli_abort(c("No existing plan found from which to compute cores.",
                    ">" = "Add one using the {.arg existing_plan} argument to {.fun redist_map}"))

    redist.identify.cores(adj = get_adj(.data),
                          plan = as.integer(as.factor(existing)),
                          boundary = boundary, focus = focus, simplify = TRUE)
}
