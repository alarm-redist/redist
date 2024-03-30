
#' Combine factors
#'
#' @param ... factors to combine
#'
#' @return a factor
#'
#' @noRd
factor_combine <- function(...) {
    f <- rlang::list2(...)
    lvls <- lapply(f, levels)
    lvls <- purrr::reduce(lvls, union)
    factor(unlist(f), levels = lvls)
}
