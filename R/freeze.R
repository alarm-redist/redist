#' Freeze Parts of a Map
#'
#'
#' @param adj Required, zero indexed adjacency list.
#' @param freeze_row Required, logical vector where TRUE freezes and FALSE lets
#'   a precinct stay free or a vector of indices to freeze
#' @param plan A vector of district assignments, which if provided will create
#'   separate groups by district. Recommended. In \code{freeze} defaults to the
#'   existing plan, if one exists.
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
#' plan = fl25_enum$plans[, 5118]
#' freeze_id <- redist.freeze(adj = fl25_adj, freeze_row = (plan == 2),
#'                            plan = plan)
#' 
#' data(iowa)
#' map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.02)
#' map <- map %>% merge_by(freeze(cd_2010 == 1, .data = .))
#'
redist.freeze <- function(adj, freeze_row, plan = rep(1, length(adj))){
  if(missing(adj)){
    stop('Please provide an object to adj.')
  }

  if(! 'list' %in% class(adj)){
    stop('adj must be an adj list.')
  }

  if(missing(freeze_row)){
    stop('Please provide a vector to freeze_row.')
  }

  if (is.numeric(freeze_row)) {
    if(max(freeze_row) < length(adj) & all(freeze_row > 0)){
      freeze_row <- 1:length(adj) %in% freeze_row
    }
  } else if (is.logical(freeze_row)) {
    if (length(freeze_row) != length(adj)) {
      stop('freeze_row is logical but does not match adj length.')
    }
  } else {
    stop('freeze_row must be a logical vector or have numeric/integer indices.')
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
freeze <- function(freeze_row, plan, .data=cur_map()) {
    if (is.null(.data))
        stop('Must provide `.data` if not called within a dplyr verb')
    if (!inherits(.data, 'redist_map'))
        stop('`.data` must be a `redist_map` object')
    if (missing(freeze_row))
        stop('Please provide an argument to freeze_row.')

    adj <- get_adj(.data)
    if (missing(plan))
        plan <- as.integer(get_existing(.data))
    if (isFALSE(plan) || is.null(plan))
        plan <- rep(1, length(adj))

    freeze_row <- eval_tidy(enquo(freeze_row), .data)
    if (is.logical(freeze_row[1]))
        freeze_row = as.logical(freeze_row)
    else
        freeze_row = as.numeric(freeze_row)

    redist.freeze(adj = adj, freeze_row = freeze_row, plan = plan)
}

globalVariables('rn')
