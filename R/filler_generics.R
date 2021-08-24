#' @export
dplyr::filter

# redist_map generics ------------

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
filter.redist_map <- function(.data, ..., .preserve=FALSE) {
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
#' @importFrom dplyr slice
slice.redist_map <- function(.data, ...) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr group_by
group_by.redist_map <- function(.data, ..., .add = FALSE) {
    reconstruct.redist_map(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr ungroup
ungroup.redist_map <- function(x, ...) {
    reconstruct.redist_map(NextMethod(), x)
}


# redist_plans generics ------------

#' @export
#' @importFrom dplyr mutate
mutate.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr transmute
transmute.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr summarise
summarise.redist_plans = function(.data, ..., .groups=NULL) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr filter
filter.redist_plans <- function(.data, ..., .preserve=FALSE) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr arrange
arrange.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr distinct
distinct.redist_plans <- function(.data, ..., .keep_all = FALSE) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr full_join
full_join.redist_plans <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_plans(NextMethod(), x)
}

#' @export
#' @importFrom dplyr inner_join
inner_join.redist_plans <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_plans(NextMethod(), x)
}

#' @export
#' @importFrom dplyr left_join
left_join.redist_plans <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_plans(NextMethod(), x)
}

#' @export
#' @importFrom dplyr right_join
right_join.redist_plans <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"), ...) {
    reconstruct.redist_plans(NextMethod(), x)
}

#' @export
#' @importFrom dplyr select
select.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr rename
rename.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr slice
slice.redist_plans <- function(.data, ...) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr group_by
group_by.redist_plans <- function(.data, ..., .add = FALSE) {
    reconstruct.redist_plans(NextMethod(), .data)
}

#' @export
#' @importFrom dplyr ungroup
ungroup.redist_plans <- function(x, ...) {
    reconstruct.redist_plans(NextMethod(), x)
}
