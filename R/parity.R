#' Calculates Population Parity
#'
#' \code{redist.parity} computes the population parity of a matrix of maps.
#'
#' @param district_membership A matrix with one row
#' for each precinct and one column for each map. Required.
#' @param population A numeric vector with the population for every precinct.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#'
#' @importFrom foreach %do% %dopar% foreach
#' @return numeric vector with the population parity for each column
#'
#' @export
redist.parity <- function(district_membership, population, ncores = 1) {
    if(!any(class(population) %in% c('numeric', 'integer'))) {
        stop('Please provide "population" as a numeric vector.')
    }

    if(!any(class(district_membership) %in% c('numeric', 'matrix'))) {
        stop('Please provide "district_membership" as a matrix.')
    }

    if(!is.matrix(district_membership)) {
        district_membership <- as.matrix(district_membership)
    }

    if(length(population) != nrow(district_membership)) {
        stop('Arguments "district_membership" and "population" do not have same number of precincts.')
    }

    # parallize as in fastLink package to avoid Windows/unix issues
    nc <- min(ncores, ncol(district_membership))
    if (nc == 1) {
        `%oper%` <- `%do%`
    } else {
        `%oper%` <- `%dopar%`
        cl <- makeCluster(nc, , setup_strategy = 'sequential')
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    if (min(district_membership[,]) == 0)
        district_membership = district_membership + 1
    n_distr = max(district_membership[,1])
    foreach(map=1:ncol(district_membership), .combine = "cbind") %oper% {
        max_dev(district_membership[,map], population, n_distr)
    }
}
