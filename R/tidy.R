##############################################
## Author: Cory McCartan
## Institution: Harvard University
## Date Created: 2021/01/28
## Purpose: redist functions for a tidy workflow
##############################################



#######################
# new redist functions


#' Check that a \code{redist_map} object is contiguous
#'
#' @param x the object
#'
#' @return \code{TRUE} if contiguous.
is_contiguous = function(x) {
    stopifnot(inherits(x, "redist_map"))
    all(contiguity(get_graph(x), rep(1, nrow(x))) == 1)
}


## merge helpers
# checks if is a proportion/pct
is_prop = function(x) is.double(x) && min(x, na.rm=T) >= 0 && max(x, na.rm=T) <= 1
# checks if is not a  proportion/pct
is_nonprop = function(x) is.numeric(x) && (min(x, na.rm=T) < 0 || max(x, na.rm=T) > 1)
# checks if x is constant within levels of rel
is_const_rel = function(rel) {
    function(x) {
        !is.numeric(x) && all(tapply(x, rel, FUN=function(y) length(unique(y))) == 1)
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
#' @param key \code{\link[tidyr:tidyr_tidy_select]{<tidy-select>}} the column to merge by
#' @param by_existing if an existing assignment is present, whether to also group by it
#' @param drop_geom whether to drop the geometry column. Recommended, as
#'   otherwise a costly geometric merge is required.
#'
#' @returns A merged \code{\link{redist_map}} object
#'
#' @export
merge_by = function(.data, key, by_existing=TRUE, drop_geom=TRUE) {
    .data = as_redist_map(.data)

    key_val = rlang::eval_tidy(rlang::enquo(key), .data)
    pop_col = attr(.data, "pop_col")

    if (drop_geom)
        .data = sf::st_drop_geometry(.data)

    col = attr(.data, "existing_col")
    unique_chr = function(x) paste(unique(x), collapse="~")
    if (!is.null(col) && by_existing) {
        dplyr::group_by(.data, {{ key }}, dplyr::across(dplyr::all_of(col))) %>%
            dplyr::summarize(dplyr::across(where(is_nonprop), sum, na.rm=T),
                             dplyr::across(where(is_prop), ~ weighted.mean(., {{ pop_col }}, na.rm=T)),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is.character), unique_chr))
    } else {
        dplyr::group_by(.data, {{ key }}) %>%
            dplyr::summarize(dplyr::across(where(is_nonprop), sum, na.rm=T),
                             dplyr::across(where(is_prop), ~ weighted.mean(., {{ pop_col }}, na.rm=T)),
                             dplyr::across(where(is_const_rel(key_val)), ~ .[1]),
                             dplyr::across(where(is.character), unique_chr)) %>%
            `attr<-`("existing_col", NULL)
    }
}

#' @rdname redist.identify.cores
#'
#' @param .data a \code{\link{redist_map}} object
#' @param within the core is defined to be at least this number of steps within
#'   district boundaries
#'
#' @export
make_cores = function(.data=get0(".", parent.frame()), within=1, focus=NULL) {
    if (is.null(.data))
        stop("Must provide `.data` if not called within a pipe")

    redist.identify.cores(get_graph(.data), get_existing(.data), within, focus, simplify=TRUE)
}


