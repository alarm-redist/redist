#' All Partitions of 25 Precincts into 3 Congressional Districts
#' (No Population Constraint)
#'
#' This data set contains demographic and geographic information about 25
#' contiguous precincts in the state of Florida. The data lists all possible
#' partitions of the 25 precincts into three contiguous congressional districts.
#' The 25-precinct shapefile may be found in \code{\link{fl25}}
#'
#' @name fl25_enum
#' @usage data("fl25_enum")
#' @format A list with two entries:
#' \describe{
#' \item{\code{plans}}{A matrix containing every partition of the 25 precincts
#' into three contiguous congressional districts, with no population constraint.}
#' \item{\code{pop_dev}}{A vector containing the maximum population deviation
#' across the three districts for each plan.}
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
#' @concept data
#' @examples 
#' data(fl25_enum)
#' 
NULL


#' Florida 25 Precinct Shape File
#'
#' This data set contains the 25-precinct shapefile and related data for each precinct.
#' All possible partitions of the 25 precincts into three contiguous
#' congressional districts are stored in \code{\link{fl25_enum}}, and the
#' corresponding adjacency graph is stored in \code{\link{fl25_adj}}.
#' This is generally useful for demonstrating basic algorithms locally.
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
#' @concept data
#' @examples 
#' data(fl25)
#' 
NULL

#' Florida 25 Precinct File
#'
#' This data set contains the 25-precinct shapefile and related data for each precinct.
#' All possible partitions of the 25 precincts into three contiguous
#' congressional districts are stored in \code{\link{fl25_enum}}, and the
#' corresponding adjacency graph is stored in \code{\link{fl25_adj}}.
#'
#' @name fl25_adj
#' @format A list storing the adjacency graph for the 25-precinct subset of Florida.
#'
#' @references Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr.
#' (2016) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo."
#' Working Paper. Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#' @concept data
#' @examples 
#' data(fl25_adj)
#' 
NULL




#' Florida 70 Precinct Shape File
#'
#' This data set contains the 70 Precinct shapefile and related data for each precinct.
#'
#' It is a random 70 precinct connected subset from Florida's precincts. This was introduced by
#' <doi:10.1080/2330443X.2020.1791773>
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
#' @references Benjamin Fifield, Kosuke Imai, Jun Kawahara & Christopher T. Kenny (2020)
#' The Essential Role of Empirical Validation in Legislative Redistricting Simulation,
#' Statistics and Public Policy, 7:1, 52-68, doi:10.1080/2330443X.2020.1791773
#' @concept data
#' @examples 
#' data(fl70)
#' 
NULL

#' Florida 250 Precinct Shape File
#'
#' This data set contains the 250 Precinct shapefile and related data for each precinct.
#'
#' It is a random 70 precinct connected subset from Florida's precincts. This was introduced by
#' <doi:10.1080/2330443X.2020.1791773>
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
#' @references Benjamin Fifield, Kosuke Imai, Jun Kawahara & Christopher T. Kenny (2020)
#' The Essential Role of Empirical Validation in Legislative Redistricting Simulation,
#' Statistics and Public Policy, 7:1, 52-68, doi:10.1080/2330443X.2020.1791773
#' @concept data
#' @examples
#' data(fl250)
#' 
NULL

#' Iowa County File
#'
#' This data contains geographic and demographic information on the 99 counties
#' of the state of Iowa.
#'
#' @name iowa
#' @usage data("iowa")
#' @format sf tibble containing columns for useful data related to the
#'   redistricting process
#' \describe{
#' \item{\code{fips}}{The FIPS code for the county.}
#' \item{\code{cd_2010}}{The 2010 congressional district assignments.}
#' \item{\code{pop}}{The total population of the precinct, according to the 2010 Census.}
#' \item{\code{white}}{The non-Hispanic white population of the precinct.}
#' \item{\code{black}}{The non-Hispanic Black population of the precinct.}
#' \item{\code{hisp}}{The Hispanic population (of any race) of the precinct.}
#' \item{\code{vap}}{The voting-age population of the precinct.}
#' \item{\code{wvap}}{The white voting-age population of the precinct.}
#' \item{\code{bvap}}{The Black voting-age population of the precinct.}
#' \item{\code{hvap}}{The Hispanic voting-age population of the precinct.}
#' \item{\code{tot_08}}{Number of total votes for president in the county in 2008.}
#' \item{\code{dem_08}}{Number of votes for Barack Obama in 2008.}
#' \item{\code{rep_08}}{Number of votes for John McCain in 2008.}
#' \item{\code{region}}{The 28E agency regions for counties.}
#' \item{\code{geometry}}{The sf geometry column containing the geographic information.}
#' }
#'
#' @concept data
#' @examples
#' data(iowa)
#' print(iowa)
NULL
