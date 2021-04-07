#' Subset a shp
#'
#' Subsets a shp object along with its adjacency. Useful for running smaller analyses
#' on pairs of districts. Provide population, ndists, pop_tol, and sub_ndists to get proper
#' population parity constraints on subsets.
#'
#' @param shp  An sf object
#' @param adj A zero-indexed adjacency list. Created with
#' \code{redist.adjacency} if not supplied.
#' @param keep_rows row numbers of precincts to keep. Random submap selected if not supplied.
#' @param total_pop numeric vector with one entry for the population of each precinct.
#' @param ndists integer, number of districts in whole map
#' @param pop_tol The strength of the hard population constraint.
#' @param sub_ndists integer, number of districts in subset map
#'
#' @return a list containing the following components:
#' \item{shp}{The subsetted shp object}
#' \item{adj}{The subsetted adjacency list for shp}
#' \item{keep_rows}{The indices of the rows kept.}
#' \item{sub_ndists}{The number of districts in the subset.}
#' \item{sub_pop_tol}{The new parity constraint for a subset.}
#'
#' @concept prepare
#' @export
redist.subset <- function(shp, adj, keep_rows, total_pop, ndists,
                          pop_tol, sub_ndists) {
    if (missing(shp)) {
        stop('Please provide an argument to "shp". ',
             'Use redist.reduce.adjacency to subset adjacency lists.')
    }
    if (!('sf' %in% class(shp))) {
        stop('Please provide "shp" as an sf object.')
    }

    if (missing(adj)) {
        adj <- redist.adjacency(shp)
    }

    if (missing(keep_rows)) {
        n <- sample(1:nrow(shp), 1)
        keep_rows <- redist.random.subgraph(shp, n, adj)$keep_rows
    }

    if (!missing(total_pop) &
        !missing(ndists) & !missing(pop_tol) & !missing(sub_ndists)) {
        parpop <- sum(total_pop) / ndists
        subparpop <- sum(total_pop[keep_rows]) / sub_ndists
        subdev <-
            min(abs(subparpop - parpop * (1 - pop_tol)),
                abs(subparpop - parpop * (1 + pop_tol)))
        sub_pop_tol <- subdev / subparpop
    } else{
        sub_ndists <- NA_real_
        sub_pop_tol <- NA_real_
    }


    rlist <- list(
        shp = shp %>% dplyr::slice(keep_rows),
        adj = redist.reduce.adjacency(adj, keep_rows = keep_rows),
        keep_rows = keep_rows,
        sub_ndists = sub_ndists,
        sub_pop_tol = sub_pop_tol
    )

    return(rlist)
}
