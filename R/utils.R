
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
    lvls <- lapply(seq_along(f), function(i) {
        if (is.null(lvls[[i]])) {
            unique(f[[i]])
        } else {
            lvls[[i]]
        }
    })
    lvls <- Reduce(f = union, x = lvls, init = c())
    factor(unlist(f), levels = lvls)
}
