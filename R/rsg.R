#' Redistricting via Random Seed and Grow Algorithm
#'
#' \code{redist.rsg} generates redistricting plans using a random seed a grow
#' algorithm.  This is the non-compact districting algorithm described in Chen and
#' Rodden (2013).  The algorithm can provide start values for the other
#' redistricting routines in this package.
#'
#'
#' @param adj List of length N, where N is the number of precincts.
#' Each list element is an integer vector indicating which precincts that precinct
#' is adjacent to.  It is assumed that precinct numbers start at 0.
#' @param total_pop numeric vector of length N, where N is the number of precincts.
#' Each element lists the population total of the corresponding precinct, and is
#' used to enforce population constraints.
#' @param ndists  integer, the number of districts we want to partition the
#' precincts into.
#' @param pop_tol numeric, indicating how close district population targets have
#' to be to the target population before algorithm converges.  thresh=0.05 for
#' example means that all districts must be between 0.95 and 1.05 times the size
#' of target.pop in population size.
#' @param verbose  boolean, indicating whether the time to run the algorithm is
#' printed.
#' @param maxiter  integer, indicating maximum number of iterations to attempt
#' before convergence to population constraint fails.  If it fails once, it will
#' use a different set of start values and try again.  If it fails again,
#' redist.rsg() returns an object of all NAs, indicating that use of more
#' iterations may be advised.
#'
#' @return list, containing three objects containing the completed redistricting
#' plan.
#' \itemize{
#' \item{\code{plan}}{  A vector of length N, indicating the
#' district membership of each precinct.}
#' \item{\code{district_list}}{  A list of length Ndistrict.  Each list contains a
#' vector of the precincts in the respective district.}
#' \item{\code{district_pop}}{  A vector of length Ndistrict, containing the
#' population totals of the respective districts.}
#' }
#'
#' @author Benjamin Fifield, Department of Politics, Princeton University
#' \email{benfifield@@gmail.com}, \url{https://www.benfifield.com/}
#'
#' Michael Higgins, Department of Statistics, Kansas State University
#' \email{mikehiggins@@k-state.edu},
#' \url{http://www-personal.k-state.edu/~mikehiggins/}
#'
#' Kosuke Imai, Department of Politics, Princeton University
#' \email{imai@@harvard.edu}, \url{http://imai.fas.harvard.edu}
#'
#' James Lo, \email{jameslo@@princeton.edu}
#'
#' Alexander Tarr, Department of Electrical Engineering, Princeton University
#' \email{atarr@@princeton.edu}
#'
#' @references Jowei Chen and Jonathan Rodden (2013) ``Unintentional
#' Gerrymandering: Political Geography and Electoral Bias in Legislatures.''
#' Quarterly Journal of Political Science. 8(3): 239-269.
#'
#' @examples
#' ### Real data example from test set
#' data(fl25)
#' data(fl25_adj)
#'
#' res <- redist.rsg(adj = fl25_adj, total_pop = fl25$pop,
#'  ndists = 3, pop_tol = 0.05)
#'
#' @concept simulate
#' @export
redist.rsg <- function(adj, total_pop, ndists, pop_tol,
                       verbose = TRUE, maxiter=5000) {

    if(verbose){
        divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

        cat("\n")
        cat(divider)
        cat("redist.rsg(): Automated Redistricting Starts\n\n")
    }

    target.pop <- sum(total_pop) / ndists

    ## Main Call to Computation - if returning NA, break.
    ## If returning districts but not contiguous, repeat
    ## First attempt
    time <- system.time(ret <- .Call('_redist_rsg',
                                     PACKAGE = 'redist',
                                     adj,
                                     total_pop,
                                     ndists,
                                     target.pop,
                                     pop_tol,
                                     as.integer(maxiter)
                                     ))
    ## Make another call if stuck, but only do one more try
    ## because maxiter might be too low
    if(is.na(ret$plan[1])){
        time <- system.time(ret <- .Call('_redist_rsg',
                                         PACKAGE = 'redist',
                                         adj,
                                         total_pop,
                                         ndists,
                                         target.pop,
                                         pop_tol,
                                         as.integer(maxiter)
                                         ))
    }

    if(is.na(ret$plan[1])){

        stop("redist.rsg() failed to return a valid partition. Try increasing maxiterrsg")

    } else {
        ret$plan <- ret$plan + 1
    }

    if(verbose){
        cat(paste("\n\t", ndists, " districts built using ",
                  length(adj), " precincts in ",
                  round(time[3], digits=2), " seconds...\n\n", sep = ""), append = TRUE)
    }

    return(ret)

}

