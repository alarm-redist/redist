#' @export
#' @importFrom dplyr mutate
mutate.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
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
