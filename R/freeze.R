#' Freeze Parts of a Map
#'
#'
#' @param adj Required, zero indexed adjacency list.
#' @param freeze_row Required, logical vector where TRUE freezes and FALSE lets
#' a precinct stay free or a vector of indices to freeze
#' @param plan A vector of district assignments, which if provided will create
#' separate groups by district. Recommended. In \code{freeze} defaults to the
#' existing plan, if one exists.
#'
#' @return integer vector to group by
#' @export
#' @concept prepare
#' @examples
#' library(redist)
#' library(dplyr)
#' data(fl25)
#' data(fl25_enum)
#' data(fl25_adj)
#' plan <- fl25_enum$plans[, 5118]
#' freeze_id <- redist.freeze(adj = fl25_adj, freeze_row = (plan == 2),
#'     plan = plan)
#'
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.02)
#' map <- map %>% merge_by(freeze(cd_2010 == 1, .data = .))
#'
redist.freeze <- function(adj, freeze_row, plan = rep(1, length(adj))) {
    if (missing(adj)) {
        cli_abort("Please provide an object to {.arg adj}.")
    }

    if (!"list" %in% class(adj)) {
        cli_abort("{.arg adj} must be an adjacency list.")
    }

    if (missing(freeze_row)) {
        cli_abort("Please provide a vector to {.arg freeze_row}.")
    }

    if (is.numeric(freeze_row)) {
        if (max(freeze_row) < length(adj) & all(freeze_row > 0)) {
            freeze_row <- 1:length(adj) %in% freeze_row
        }
    } else if (is.logical(freeze_row)) {
        if (length(freeze_row) != length(adj)) {
            cli_abort("{.arg freeze_row} is logical but does not match {.arg adj} length.")
        }
    } else {
        stop("{.arg freeze_row} must be a logical vector or have numeric/integer indices.")
    }

    tb <- tibble(plan = plan, freeze_row = freeze_row) %>%
        group_by(plan, freeze_row) %>%
        mutate(id = cur_group_id()) %>%
        ungroup()
    cont <- contiguity(adj = adj, group = tb$id)
    tb <- tb %>% mutate(cont = cont) %>%
        mutate(rn = row_number()) %>%
        group_by(plan, cont, freeze_row) %>%
        mutate(rn = ifelse(freeze_row, min(rn), rn)) %>%
        ungroup() %>%
        group_by(rn) %>%
        mutate(gid = cur_group_id())

    return(tb$gid)
}

#' @rdname redist.freeze
#' @order 1
#'
#' @param .data a \code{\link{redist_map}} object
#'
#' @concept prepare
#' @export
freeze <- function(freeze_row, plan, .data = cur_map()) {
    if (is.null(.data))
        cli_abort("Must provide {.arg .data} if not called within a {.pkg dplyr} verb")
    if (!inherits(.data, "redist_map"))
        cli_abort("{.arg .data} must be a {.cls redist_map}")
    if (missing(freeze_row))
        cli_abort("{.arg freeze_row} cannot be missing.")

    adj <- get_adj(.data)
    if (missing(plan))
        plan <- as.integer(get_existing(.data))
    if (isFALSE(plan) || is.null(plan))
        plan <- rep(1, length(adj))

    freeze_row <- eval_tidy(enquo(freeze_row), .data)
    if (is.logical(freeze_row[1]))
        freeze_row <- as.logical(freeze_row)
    else
        freeze_row <- as.numeric(freeze_row)

    redist.freeze(adj = adj, freeze_row = freeze_row, plan = plan)
}

globalVariables("rn")
