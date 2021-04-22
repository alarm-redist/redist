#' Calculates Maximum Deviation from Population Parity
#'
#' Computes the deviation from population parity from a plan.
#' Higher values indicate that (at least) a single district in the map deviates
#' from population parity. See Details.
#'
#' @details With a map with \code{pop} representing the populations of each district,
#'  the deviation from population parity is given as \code{max(abs(pop - parity) / parity)}
#'  where \code{parity = sum(pop)/length(pop)} is the population size for the
#'  average district.
#'  Therefore, the metric can be thought of as the maximum percent deviation from
#'  equal population. For example, a value of 0.03 in this metric indicates that
#'  all districts are within 3 percent of population parity.
#'
#' @param plans A matrix with one row for each precinct and one column for each
#'   map. Required.
#' @param total_pop A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#' @importFrom foreach %do% %dopar% foreach
#' @return numeric vector with the population parity for each column
#'
#' @concept analyze
#' @export
redist.parity <- function(plans, total_pop, ncores = 1) {
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
        cl <- makeCluster(nc, setup_strategy = 'sequential', methods=FALSE)
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
