#' Adjacency List functionality for redist
#'
#' @param shp A SpatialPolygonsDataFrame or sf object. Required.
#' @param plan A numeric vector (if only one map) or matrix with one row
#'
#' @return Adjacency list
#' @description Creates an adjacency list that is zero indexed with no skips
#'
#'
#' @importFrom sf st_relate
#' @export
redist.adjacency <- function(shp, plan) {
    # Check input
    if (!any(c("sf", "SpatialPolygonsDataFrame") %in% class(shp))) {
        cli_abort("{.arg shp} must be a {.cls sf} or {.cls sp} object")
    }

    # Create the adjacency with spdep function
    # Get standard rooks contiguity
    adj <- suppressMessages(st_relate(shp, shp, pattern = "F***1****"))
    # items contained entirely within ~ even if validly "rooks" adjacent ~ do not meet this, you need:
    withinadj <- suppressMessages(st_relate(x = shp, pattern = "2121**2*2"))
    adj <- lapply(1:nrow(shp), function(x) c(adj[[x]], withinadj[[x]]))

    # Check for zero indexing
    zero <- min(unlist(adj)) == 0

    # Make zero indexed if not
    min <- min(unlist(adj))
    if (!zero) {
        adj <- lapply(adj, function(x) {x - min})
    }

    # Check that no numbers are skipped
    # low resolution shp files may result in skips, this fixes most issues
    correct_n <- nrow(shp) == length(unique(unlist(adj)))

    if (!correct_n) {
        cli_warn("At least one precinct had no adjacent precincts.")
    } else if (any(contiguity(adj = adj, group = rep(1, length(adj))) > 1)) {
        cli_warn("All precincts have at least one neighbor, but the graph is disconnected.")
    }


    if (!missing(plan)) {
        cont <- contiguity(adj, plan)
        if (any(cont > 1)) {
            cli_warn("District {unique(plan[cont>1])} was not contiguous.")
        }
    }

    # return a checked adjacency list
    adj
}

#' Reduce Adjacency List
#'
#' Tool to help reduce adjacency lists for analyzing subsets of maps.
#'
#' @param adj A zero-indexed adjacency list. Required.
#' @param keep_rows row numbers of precincts to keep
#'
#' @return zero indexed adjacency list with max value length(keep_rows) - 1
#'
#' @concept prepare
#' @export
#'
#' @examples
#' data(fl25_adj)
#' redist.reduce.adjacency(fl25_adj, c(2, 3, 4, 6, 21))
#'
redist.reduce.adjacency <- function(adj, keep_rows) {
    # Check inputs:
    if (!(class(keep_rows) %in% c("numeric", "integer"))) {
        cli_warn("{.arg keep_rows} must be a numeric or integer vector.")
    }
    if (min(unlist(adj)) != 0) {
        cli_abort("{.arg adj} must be a 0-indexed list.")
    }
    if (max(unlist(adj)) != (length(adj) - 1)) {
        cli_warn("{.arg adj} did not have typical values of 0:(length(adj)-1)")
    }

    prec_map <- match(seq_along(adj), keep_rows) - 1L
    prec_map[is.na(prec_map)] = -1L

    # Reduce!
    reduce_adj(adj, prec_map, length(keep_rows))
}


#' Coarsen Adjacency List
#'
#' @param adj A zero-indexed adjacency list. Required.
#' @param groups integer vector of elements of adjacency to group
#'
#' @return adjacency list coarsened
#'
#' @concept prepare
#' @export
redist.coarsen.adjacency <- function(adj, groups) {
    if (min(unlist(adj)) != 0) {
        cli_abort("{.arg adj} must be a 0-indexed list.")
    }
    if (max(unlist(adj)) != (length(adj) - 1)) {
        cli_warn("{.arg adj} did not have typical values of 0:(length(adj)-1)")
    }
    if (length(groups) != length(adj)) {
        cli_abort("{.arg groups} and {.arg adj} have different sizes.")
    }
    if (min(groups) != 0) {
        groups <- groups - min(groups)
    }

    groups <- as.integer(groups)

    coarsen_adjacency(adj, groups)
}

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
        cli_abort(c("{.arg shp} is required.",
            "i" = "Use {.fn redist.reduce.adjacency} to subset adjacency lists."))
    }
    if (!inherits(shp, "sf")) {
        cli_abort("{.arg shp} must be an {.cls sf} object.")
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
        parpop <- sum(total_pop)/ndists
        subparpop <- sum(total_pop[keep_rows])/sub_ndists
        subdev <-
            min(abs(subparpop - parpop*(1 - pop_tol)),
                abs(subparpop - parpop*(1 + pop_tol)))
        sub_pop_tol <- subdev/subparpop
    } else {
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

    rlist
}
