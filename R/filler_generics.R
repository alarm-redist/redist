#' @export
#' @importFrom dplyr mutate
mutate.redist_map <- function(.data, ..., env=rlang::caller_env()) {
    dots = rlang::quos(...)
    for (i in seq_along(dots)) {
        dot_call = rlang::call_standardise(rlang::quo_get_expr(dots[[i]]))
        if (!(".data" %in% names(dot_call))) {
            dot_call$.data = .data
        }
        dots[[i]] = rlang::quo_set_expr(dots[[i]], dot_call)
    }

    .data_ex = rlang::enexpr(.data)
    mut_call = rlang::expr(getS3method("mutate", "data.frame")(!!.data_ex, !!!dots))
    res = eval(mut_call, env)

    if (inherits(.data, "sf")) {
        agr = sf::st_agr(.data)
        sf_column_name = attr(.data, "sf_column")
        attr(res, "sf_column") = sf_column_name
        attr(res, "agr") = agr
        class(res) = c("sf", class(res))
    }

    reconstruct.redist_map(res, .data)

}

#' @export
#' @importFrom dplyr transmute
transmute.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr filter
filter.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr arrange
arrange.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr distinct
distinct.redist_map <- function(.data, ..., .keep_all = FALSE) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr full_join
full_join.redist_map <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_map(NextMethod(), x)
}

#' @export
#' @importFrom dplyr inner_join
inner_join.redist_map <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_map(NextMethod(), x)
}

#' @export
#' @importFrom dplyr left_join
left_join.redist_map <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_map(NextMethod(), x)
}

#' @export
#' @importFrom dplyr right_join
right_join.redist_map <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_map(NextMethod(), x)
}

#' @export
#' @importFrom dplyr select
select.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr slice
slice.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr group_by
group_by.redist_map <- function(.data, ..., add = FALSE) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr ungroup
ungroup.redist_map <- function(x, ...) {
    reconstruct.redist_map(NextMethod(), x)
}
