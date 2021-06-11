#' Redistricting via Compact Random Seed and Grow Algorithm
#'
#' \code{redist.crsg} generates redistricting plans using a random seed a grow
#' algorithm.  This is the compact districting algorithm described in Chen and
#' Rodden (2013).
#'
#' @param adj List of length N, where N is the number of precincts.
#' Each list element is an integer vector indicating which precincts that precinct
#' is adjacent to.  It is assumed that precinct numbers start at 0.
#' @param total_pop numeric vector of length N, where N is the number of precincts.
#' Each element lists the population total of the corresponding precinct, and is
#' used to enforce pop_tol constraints.
#' @param shp An sf dataframe to compute area and centroids with.
#' @param ndists  integer, the number of districts we want to partition the
#' precincts into.
#' @param pop_tol numeric, indicating how close district population targets have
#' to be to the target population before algorithm converges.  pop_tol=0.05 for
#' example means that all districts must be between 0.95 and 1.05 times the size
#' of target.pop in population size.
#' @param verbose  boolean, indicating whether the time to run the algorithm is
#' printed.
#' @param maxiter  integer, indicating maximum number of iterations to attempt
#' before convergence to population constraint fails.  If it fails once, it will
#' use a different set of start values and try again.  If it fails again,
#' redist.rsg() returns an object of all NAs, indicating that use of more
#' iterations may be advised. Default is 5000.
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
#'
#' @references Jowei Chen and Jonathan Rodden (2013) ``Unintentional
#' Gerrymandering: Political Geography and Electoral Bias in Legislatures.''
#' Quarterly Journal of Political Science. 8(3): 239-269.
#'
#' @examples
#' data("fl25")
#' adj <- redist.adjacency(fl25)
#' redist.crsg(adj = adj, total_pop = fl25$pop, shp = fl25, ndists = 2, pop_tol = .1)
#'
#' @concept simulate
#' @export
redist.crsg <- function(adj, total_pop, shp,  ndists, pop_tol, verbose = TRUE,
                        maxiter = 5000){
  if (missing(shp)) {
    stop('An argument to shp is now required.')
  }

  suppressWarnings(coords <- st_coordinates(st_centroid(shp)))
  x_center <- coords[,1]
  y_center <- coords[,2]
  area <- st_area(shp)


  if(verbose){
    divider <- c(paste(rep("=", 20), sep = "", collapse = ""), "\n")

    cat("\n")
    cat(divider)
    cat("redist.crsg(): Automated Redistricting Starts\n\n")
  }

  target.pop <- sum(total_pop) / ndists

  # Main call to function - unlike rsg, uses RcppExporting not direct .Call
  time <- system.time(ret <- crsg(adj_list = adj,
                                  population = total_pop,
                                  area = area,
                                  x_center = x_center,
                                  y_center = y_center,
                                  Ndistrict = ndists,
                                  target_pop = target.pop,
                                  thresh = pop_tol,
                                  maxiter = maxiter))

  ## Make another call if NA, but beware this may be due to maxiter.
  ## This could also be if there are no valid moves.
  if(is.na(ret$plan[1])){
    time <- system.time(ret <- crsg(adj_list = adj,
                                    population = total_pop,
                                    area = area,
                                    x_center = x_center,
                                    y_center = y_center,
                                    Ndistrict = ndists,
                                    target_pop = target.pop,
                                    thresh = pop_tol,
                                    maxiter = maxiter))
  }


  if(is.na(ret$plan[1])){

    warning("redist.crsg() failed to return a valid partition. Try increasing maxiter")

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
