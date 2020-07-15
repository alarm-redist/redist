#' All Partitions of 25 Precincts into 3 Congressional Districts
#' (No Population Constraint)
#'
#' This data set contains demographic and geographic information about 25
#' contiguous precincts in the state of Florida. The data lists all possible
#' partitions of the 25 precincts into three contiguous congressional districts.
#'
#' @name algdat.pfull
#' @usage data("algdat.pfull")
#' @format A list with five entries:
#' \describe{
#' \item{\code{adjlist}}{An adjacency list for the 25 precincts.}
#' \item{\code{cdmat}}{A matrix containing every partition of the 25 precincts
#' into three contiguous congressional districts, with no population constraint.}
#' \item{\code{precinct.data}}{A matrix containing demographic information for
#' each of the 25 precincts.}
#' \item{\code{segregation.index}}{A matrix containing the dissimilarity index of
#' segregation (Massey and Denton 1987) for each congressional district map in
#' \code{cdmat}.}
#' \item{\code{distancemat}}{A square matrix containing the squared distance
#' between the centroids of any two precincts.}
#' }
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social Segregation".
#' Social Forces.
#'
#' @examples \dontrun{
#' data(algdat.pfull)
#' }
NULL

#' All Partitions of 25 Precincts into 3 Congressional Districts
#' (10\% Population Constraint)
#'
#' This data set contains demographic and geographic information about 25
#' contiguous precincts in the state of Florida. The data lists all possible
#' partitions of the 25 precincts into three contiguous congressional districts,
#' conditional on the congressional districts falling within 10\% of population
#' parity.
#'
#' @name algdat.p10
#' @usage data("algdat.p10")
#' @format A list with five entries:
#' \describe{
#' \item{\code{adjlist}}{An adjacency list for the 25 precincts.}
#' \item{\code{cdmat}}{A matrix containing every partition of the 25 precincts
#' into three contiguous congressional districts, with no population constraint.}
#' \item{\code{precinct.data}}{A matrix containing demographic information for
#' each of the 25 precincts.}
#' \item{\code{segregation.index}}{A matrix containing the dissimilarity index of
#' segregation (Massey and Denton 1987) for each congressional district map in
#' \code{cdmat}.}
#' \item{\code{distancemat}}{A square matrix containing the squared distance
#' between the centroids of any two precincts.}
#' }
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social Segregation".
#' Social Forces.
#'
#' @examples \dontrun{
#' data(algdat.p10)
#' }
NULL

#' All Partitions of 25 Precincts into 3 Congressional Districts
#' (20\% Population Constraint)
#'
#' This data set contains demographic and geographic information about 25
#' contiguous precincts in the state of Florida. The data lists all possible
#' partitions of the 25 precincts into three contiguous congressional districts,
#' conditional on the congressional districts falling within 20\% of population
#' parity.
#'
#' @name algdat.p20
#' @usage data("algdat.p20")
#' @format A list with five entries:
#' \describe{
#' \item{\code{adjlist}}{An adjacency list for the 25 precincts.}
#' \item{\code{cdmat}}{A matrix containing every partition of the 25 precincts
#' into three contiguous congressional districts, with no population constraint.}
#' \item{\code{precinct.data}}{A matrix containing demographic information for
#' each of the 25 precincts.}
#' \item{\code{segregation.index}}{A matrix containing the dissimilarity index of
#' segregation (Massey and Denton 1987) for each congressional district map in
#' \code{cdmat}.}
#' \item{\code{distancemat}}{A square matrix containing the squared distance
#' between the centroids of any two precincts.}
#' }
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social Segregation".
#' Social Forces.
#'
#' @examples \dontrun{
#' data(algdat.p20)
#' }
NULL

#' Florida 25 Precinct File
#'
#' This data set contains the 25 Precinct shapefile and related data for each precinct.
#'
#' @name fl25
#' @usage data("fl25")
#' @format sf data.frame containing columns for useful data related to the 
#' redistricting process, subsetted from real data in Florida, and sf geometry column.
#' \describe{
#' \item{\code{geoid}}{ Contains unique identifier for each precinct which can be matched to the full Florida dataset.}
#' \item{\code{pop}}{ Contains the population of each precinct.}
#' \item{\code{vap}}{ Contains the voting age population of each precinct.}
#' \item{\code{obama}}{ Contains the 2012 presidential vote for Obama.}
#' \item{\code{mccain}}{ Contains the 2012 presidential vote for McCain.}
#' \item{\code{TotPop}}{ Contains the population of each precinct. Identical to pop.}
#' \item{\code{BlackPop}}{Contains the black population of each precinct.}
#' \item{\code{HispPop}}{Contains the Hispanic population of each precinct.}
#' \item{\code{VAP}}{ Contains the voting age population of each precinct. Identical to vap.}
#' \item{\code{BlackVAP}}{ Contains the voting age population of black constituents of each precinct.}
#' \item{\code{HispVAP}}{ Contains the voting age population of hispanic constituents of each precinct.}
#' \item{\code{geometry}}{ Contains sf geometry of each precinct.}
#' }
#' 
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#' @examples \dontrun{
#' data(fl25)
#' }
NULL

#' Florida 70 Precinct File
#'
#' This data set contains the 70 Precinct shapefile and related data for each precinct.
#'
#' @name fl70
#' @usage data("fl70")
#' @format sf data.frame containing columns for useful data related to the 
#' redistricting process, subsetted from real data in Florida, and sf geometry column.
#' \describe{
#' \item{\code{geoid}}{ Contains unique identifier for each precinct which can be matched to the full Florida dataset.}
#' \item{\code{pop}}{ Contains the population of each precinct.}
#' \item{\code{vap}}{ Contains the voting age population of each precinct.}
#' \item{\code{obama}}{ Contains the 2012 presidential vote for Obama.}
#' \item{\code{mccain}}{ Contains the 2012 presidential vote for McCain.}
#' \item{\code{TotPop}}{ Contains the population of each precinct. Identical to pop.}
#' \item{\code{BlackPop}}{Contains the black population of each precinct.}
#' \item{\code{HispPop}}{Contains the Hispanic population of each precinct.}
#' \item{\code{VAP}}{ Contains the voting age population of each precinct. Identical to vap.}
#' \item{\code{BlackVAP}}{ Contains the voting age population of black constituents of each precinct.}
#' \item{\code{HispVAP}}{ Contains the voting age population of hispanic constituents of each precinct.}
#' \item{\code{geometry}}{ Contains sf geometry of each precinct.}
#' }
#' 
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#' @examples \dontrun{
#' data(fl70)
#' }
NULL

#' Florida 250 Precinct File
#'
#' This data set contains the 250 Precinct shapefile and related data for each precinct.
#'
#' @name fl250
#' @usage data("fl250")
#' @format sf data.frame containing columns for useful data related to the 
#' redistricting process, subsetted from real data in Florida, and sf geometry column.
#' \describe{
#' \item{\code{geoid}}{ Contains unique identifier for each precinct which can be matched to the full Florida dataset.}
#' \item{\code{pop}}{ Contains the population of each precinct.}
#' \item{\code{vap}}{ Contains the voting age population of each precinct.}
#' \item{\code{obama}}{ Contains the 2012 presidential vote for Obama.}
#' \item{\code{mccain}}{ Contains the 2012 presidential vote for McCain.}
#' \item{\code{TotPop}}{ Contains the population of each precinct. Identical to pop.}
#' \item{\code{BlackPop}}{Contains the black population of each precinct.}
#' \item{\code{HispPop}}{Contains the Hispanic population of each precinct.}
#' \item{\code{VAP}}{ Contains the voting age population of each precinct. Identical to vap.}
#' \item{\code{BlackVAP}}{ Contains the voting age population of black constituents of each precinct.}
#' \item{\code{HispVAP}}{ Contains the voting age population of hispanic constituents of each precinct.}
#' \item{\code{geometry}}{ Contains sf geometry of each precinct.}
#' }
#' 
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#' @examples \dontrun{
#' data(fl250)
#' }
NULL