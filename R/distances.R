#' Compute Distance between Partitions
#'
#' @param district_membership A matrix with one row for each precinct and one
#' column for each map. Required.
#' @param measure String vector indicating which distances to compute. Implemented
#' currently are "Hamming", "Manhattan", "Euclidean", and "variation of information",
#' Use "all" to return all implemented measures.  Not case sensitive, and
#' any unique substring is enough, e.g. "ham" for Hamming, or "info" for
#' variation of information.
#' @param ncores Number of cores to use for parallel computing. Default is 1.
#' @param pop The vector of precinct populations. Used only if computing
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
#' districts, \eqn{\text{pop(\cdot)}} is the population, and \eqn{P} the total
#' population of the state. VI is also expressible as the difference between
#' the joint entropy and the mutual information (see references).
#'
#' @return a named list of distance matrices, one for each distance measure selected.
#'
#' @references
#' Cover, T. M. and Thomas, J. A. (2006). \emph{Elements of information theory.} John Wiley & Sons, 2 edition.
#'
#' @examples \dontrun{
#' data("algdat.p10")
#' distances <- redist.distances(district_membership = algdat.p10$cdmat)
#' distances$Hamming[1:5,1:5]
#' }
#' @export
redist.distances <- function(district_membership, measure = "Hamming",
                             ncores = 1, pop=NULL) {

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
    nc <- min(ncores, ncol(district_membership))
    if (nc == 1) {
        `%oper%` <- `%do%`
    } else {
        `%oper%` <- `%dopar%`
        cl <- makeCluster(nc, setup_strategy = 'sequential')
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
    }

    # Compute Hamming Distance Metric
    if ("Hamming" %in% measure) {
        ham <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
            hamming(v = district_membership[,map], m = district_membership)
        }
        colnames(ham) <- NULL

        done = done + 1
        distances[[done]] <- ham
        names(distances)[done] <- "Hamming"
    }


    # Compute Manhattan Distance Metric
    if ("Manhattan" %in% measure) {
        man <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
            minkowski(v = district_membership[,map], m = district_membership, p = 1)
        }
        colnames(man) <- NULL

        done = done + 1
        distances[[done]] <- man
        names(distances)[done] <- "Manhattan"
    }

    # Compute Euclidean Distance Metric
    if ("Euclidean" %in% measure) {
        euc <- foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
            minkowski(v = district_membership[,map], m = district_membership, p = 2)
        }
        colnames(euc) <- NULL

        done = done + 1
        distances[[done]] <- euc
        names(distances)[done] <- "Euclidean"
    }

    if ("variation of information" %in% measure) {
        if (is.null(pop)) {
            warning("Population not provided, using default of equal population.")
            pop = rep(1, nrow(district_membership))
        }
        if (length(pop) != nrow(district_membership))
            stop("Mismatch: length of population vector does not match the number of precincts.")

        # 1-index in preparation
        if (min(district_membership) == 0)
            district_membership = district_membership + 1

        vi = foreach(map = 1:ncol(district_membership), .combine = "cbind") %oper% {
            var_info_mat(district_membership, map-1, pop) # 0-index
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

utils::globalVariables(names = "map")
