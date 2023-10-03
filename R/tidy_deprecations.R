
#' @rdname redist.compactness
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.compactness}
#'
#' @concept analyze
#' @export
distr_compactness <- function(map, measure = "FracKept", .data = cur_plans(), ...) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)

    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    redist.compactness(shp = map, plans = get_plans_matrix(.data),
        measure = measure, total_pop = map[[attr(map, "pop_col")]],
        adj = get_adj(map), ...)[[measure]]
}

#' @rdname redist.segcalc
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
segregation_index <- function(map, group_pop, total_pop = map[[attr(map, "pop_col")]],
                              .data = cur_plans()) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)
    group_pop <- rlang::eval_tidy(rlang::enquo(group_pop), map)
    total_pop <- rlang::eval_tidy(rlang::enquo(total_pop), map)
    plan_m <- get_plans_matrix(.data)
    rep(as.numeric(redist.segcalc(plans = plan_m, group_pop = group_pop,
        total_pop = total_pop)),
    each = attr(map, "ndists"))
}

#' @rdname redist.metrics
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.metrics}
#'
#' @concept analyze
#' @export
partisan_metrics <- function(map, measure, rvote, dvote, ...,
                             .data = cur_plans()) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)
    # districts not in ascending order
    if (length(unique(diff(as.integer(.data$district)))) > 2)
        cli_warn("Districts not sorted in ascending order; output may be incorrect.")

    rvote <- rlang::eval_tidy(rlang::enquo(rvote), map)
    dvote <- rlang::eval_tidy(rlang::enquo(dvote), map)
    as.numeric(redist.metrics(plans = get_plans_matrix(.data),
        measure = measure, rvote = rvote, dvote = dvote, ...)[[measure]])
}

#' @rdname redist.competitiveness
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
competitiveness <- function(map, rvote, dvote, .data = cur_plans()) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)
    rvote <- rlang::eval_tidy(rlang::enquo(rvote), map)
    dvote <- rlang::eval_tidy(rlang::enquo(dvote), map)
    redist.competitiveness(plans = get_plans_matrix(.data),
        rvote = rvote, dvote = dvote)
}

#' @rdname redist.splits
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
county_splits <- function(map, counties, .data = cur_plans()) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)
    counties <- rlang::eval_tidy(rlang::enquo(counties), map)
    redist.splits(plans = get_plans_matrix(.data), counties = counties)
}


#' @rdname redist.muni.splits
#' @order 1
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#'
#' @concept analyze
#' @export
muni_splits <- function(map, munis, .data = cur_plans()) {
    .Deprecated("redistmetrics")
    check_tidy_types(map, .data)
    idxs <- unique(as.integer(.data$draw))
    munis <- rlang::eval_tidy(rlang::enquo(munis), map)
    redist.muni.splits(plans = get_plans_matrix(.data)[, idxs, drop = FALSE], munis = munis)
}
