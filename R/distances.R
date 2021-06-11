#' Compute Distance between Partitions
#'
#'
#' @param plans A matrix with one row for each precinct and one
#' column for each map. Required.
#' @param measure String vector indicating which distances to compute. Implemented
#' currently are "Hamming", "Manhattan", "Euclidean", and "variation of information",
#' Use "all" to return all implemented measures.  Not case sensitive, and
#' any unique substring is enough, e.g. "ham" for Hamming, or "info" for
#' variation of information.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' @param total_pop The vector of precinct populations. Used only if computing
#'   variation of information. If not provided, equal population of precincts
#'   will be assumed, i.e. the VI will be computed with respect to the precincts
#'   themselves, and not the population.
#'
#' @details
#' Hamming distance measures the number of different precinct assignments
#' between plans. Manhattan and Euclidean distances are the 1- and 2-norms for
#' the assignment vectors.  All three of the Hamming, Manhattan, and Euclidean
#' distances implemented here are not invariant to permutations of the district
#' labels; permuting will cause large changes in measured distance, and maps
#' which are identical up to a permutation may be computed to be maximally
#' distant.
#'
#' Variation of Information is a metric on population partitions (i.e.,
#' districtings) which is invariant to permutations of the district labels, and
#' arises out of information theory. It is calculated as \deqn{
#'     VI(\xi, \xi') = -\sum_{i=1}^n\sum_{j=1}^n pop(\xi_i \cap \xi'_j)/P
#'     (2log(pop(\xi_i \cap \xi'_j)) - log(pop(\xi_i)) - log(pop(\xi'_j)))
#' } where \eqn{\xi,\xi'} are the partitions, \eqn{\xi_i,\xi_j} the individual
#' districts, \eqn{pop(\cdot)} is the population, and \eqn{P} the total
#' population of the state. VI is also expressible as the difference between
#' the joint entropy and the mutual information (see references).
#'
#' @return a named list of distance matrices, one for each distance measure selected.
#'
#' @references
#' Cover, T. M. and Thomas, J. A. (2006). \emph{Elements of information theory.} John Wiley & Sons, 2 edition.
#'
#' @examples
#' data(fl25)
#' data(fl25_enum)
#'
#' plans_05 <- fl25_enum$plans[, fl25_enum$pop_dev <= 0.05]
#' distances <- redist.distances(plans_05)
#' distances$Hamming[1:5, 1:5]
#' 
#' @concept analyze
#' @export
redist.distances <- function(plans, measure = "Hamming",
                             ncores = 1, total_pop = NULL) {

    supported = c("all", "Hamming", "Manhattan", "Euclidean", "variation of information")
    # fuzzy matching
    measure = supported[agrep(measure, supported, max.distance=0, ignore.case=T)]
    #check inputs
    if (measure == "all") {
        measure <- supported[-1] # all but 'all'
    }

    # init vars
    distances <- list()
    name <- c()
    done <- 0

    # parallel setup
    nc <- min(ncores, ncol(plans))
    if (nc == 1) {
        `%oper%` <- `%do%`
    } else {
        `%oper%` <- `%dopar%`
        cl <- makeCluster(nc, setup_strategy = 'sequential', methods=FALSE)
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    # Compute Hamming Distance Metric
    if ("Hamming" %in% measure) {
        ham <- foreach(map = 1:ncol(plans), .combine = "cbind") %oper% {
            hamming(v = plans[,map], m = plans)
        }
        colnames(ham) <- NULL

        done = done + 1
        distances[[done]] <- ham
        names(distances)[done] <- "Hamming"
    }


    # Compute Manhattan Distance Metric
    if ("Manhattan" %in% measure) {
        man <- foreach(map = 1:ncol(plans), .combine = "cbind") %oper% {
            minkowski(v = plans[,map], m = plans, p = 1)
        }
        colnames(man) <- NULL

        done = done + 1
        distances[[done]] <- man
        names(distances)[done] <- "Manhattan"
    }

    # Compute Euclidean Distance Metric
    if ("Euclidean" %in% measure) {
        euc <- foreach(map = 1:ncol(plans), .combine = "cbind") %oper% {
            minkowski(v = plans[,map], m = plans, p = 2)
        }
        colnames(euc) <- NULL

        done = done + 1
        distances[[done]] <- euc
        names(distances)[done] <- "Euclidean"
    }

    if ("variation of information" %in% measure) {
        if (is.null(total_pop)) {
            warning("Population not provided, using default of equal population.")
            total_pop = rep(1, nrow(plans))
        }
        if (length(total_pop) != nrow(plans))
            stop("Mismatch: length of population vector does not match the number of precincts.")

        # 1-index in preparation
        if (min(plans) == 0)
            plans = plans + 1

        vi = foreach(map = 1:ncol(plans), .combine = "cbind") %oper% {
            var_info_mat(plans, map-1, total_pop) # 0-index
        }
        colnames(vi) <- NULL
        # copy over other half of matrix; we only computed upper triangle
        vi[lower.tri(vi)] = t(vi)[lower.tri(vi)]

        done = done + 1
        distances[[done]] <- vi
        names(distances)[done] <- "VI"
    }

    return(distances)
}

#' @rdname redist.distances
#' @order 1
#'
#' @param plans a \code{\link{redist_plans}} object.
#'
#' @returns \code{distance_matrix} returns a numeric distance matrix for the
#'   chosen metric.
#'
#' @export
plan_distances = function(plans, measure="variation of information", ncores=1) {
    choices = c("variation of information", "Hamming", "Manhattan", "Euclidean")
    measure = match.arg(measure, choices)
    pop = attr(plans, "prec_pop")
    if (is.null(pop))
        stop("Precinct population must be stored in `prec_pop` attribute of `plans` object")

    redist.distances(get_plans_matrix(plans), measure, ncores=ncores, total_pop=pop)[[1]]
}

#' Calculate the diversity of a set of plans
#'
#' Returns the off-diagonal elements of the variation of information distance
#' matrix for a sample of plans, which can be used as a diagnostic measure to
#' assess the diversity of a set of plans. If there are many values close to zero,
#' then the sample has many similar plans and may not be a good approximation
#' to the target distribution.
#'
#' @param plans a \code{\link{redist_plans}} object.
#' @param n_max the maximum number of plans to sample in computing the
#'   distances. Larger numbers will have less sampling error but will require
#'   more computation time.
#' @param ncores the number of cores to use in computing the distances.
#'
#' @return A numeric vector of off-diagonal variation of information distances.
#'
#' @examples
#' data(iowa)
#' ia <- redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01)
#' plans <- redist_smc(ia, 100, silent=TRUE)
#' hist(plans_diversity(plans))
#'
#' @concept analyze
#' @export
plans_diversity = function(plans, n_max=100, ncores=1) {
    m = get_plans_matrix(plans)
    n_eval = min(n_max, ncol(m))
    idx = sample.int(ncol(m), n_eval, replace=FALSE)

    ndists = max(m[, 1])
    pop = attr(plans, "prec_pop")
    if (is.null(pop))
        stop("Precinct population must be stored in `prec_pop` attribute of `plans` object")

    dists = redist.distances(m[, idx], "variation of information",
                             ncores=ncores, total_pop=pop)$VI
    0.5 * dists[upper.tri(dists)] / log(ndists)
}


utils::globalVariables(names = "map")
