#' Calculates Population Parity
#'
#' \code{redist.parity} computes the population parity of a matrix of maps.
#' @param plans A matrix with one row for each precinct and one column for each
#'   map. Required.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' @param district_membership Deprecated, use plans. A matrix with one row
#' for each precinct and one column for each map. Required.
#' @param population Deprecated, use total_pop. A numeric vector with the population for every precinct.
#'
#' @importFrom foreach %do% %dopar% foreach
#' @return numeric vector with the population parity for each column
#'
#' @concept analyze
#' @export
redist.parity <- function(plans, total_pop, ncores = 1, district_membership, population) {
    if (!missing(population)) {
        .Deprecated(new = 'total_pop', old = 'population')
        total_pop <- population
    }
    if (!missing(district_membership)) {
        .Deprecated(new = 'plans', old = 'district_membership')
        plans <- district_membership
    }

    if (!any(class(total_pop) %in% c('numeric', 'integer'))) {
        stop('Please provide "total_pop" as a numeric vector.')
    }
    if (!is.matrix(plans)) {
        plans <- matrix(plans, ncol=1)
    }
    if (!any(class(plans) %in% c('numeric', 'matrix'))) {
        stop('Please provide "plans" as a matrix.')
    }

    if (length(total_pop) != nrow(plans)) {
        stop('Arguments "plans" and "total_pop" do not have same number of precincts.')
    }

    # parallelize as in fastLink package to avoid Windows/unix issues
    N = ncol(plans)
    nc <- min(ncores, max(1, floor(N/2)))

    if (nc == 1){
        `%oper%` <- `%do%`
    } else {
        `%oper%` <- `%dopar%`
        cl <- makeCluster(nc, setup_strategy = 'sequential')
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    if (min(plans[,1]) == 0)
        plans = plans + 1
    n_distr = max(plans[,1])


    chunks = split(1:N, rep(1:nc, each=ceiling(N/nc))[1:N])
    out = foreach(map=chunks, .combine = "c") %oper% {
        max_dev(plans[, map, drop = FALSE], total_pop, n_distr)
    }

    unlist(out)
}

#' @rdname redist.parity
#'
#' @param map a \code{\link{redist_map}} object
#' @param .data a \code{\link{redist_plans}} object
#' @param ... passed on to \code{redist.parity}
#'
#' @concept analyze
#' @export
plan_parity <- function(map, .data = get0('.', parent.frame()), ...) {
  if (!inherits(.data, 'redist_plans')) {
    stop('Must provide `.data` if not called within a pipe')
  }
  ndists <- attr(map, 'ndists')
  total_pop <- map[[attr(map, 'pop_col')]]
  stopifnot(!is.null(total_pop))

  rep(max_dev(get_plan_matrix(.data), total_pop, ndists), each = ndists)
}
