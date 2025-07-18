% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_helpers.R, R/cores.R
\name{make_cores}
\alias{make_cores}
\alias{redist.identify.cores}
\title{Identify Cores of a District (Heuristic)}
\usage{
make_cores(.data = cur_map(), boundary = 1, focus = NULL)

redist.identify.cores(adj, plan, boundary = 1, focus = NULL, simplify = TRUE)
}
\arguments{
\item{.data}{a \code{\link{redist_map}} object}

\item{boundary}{Number of steps to check for. Defaults to 1.}

\item{focus}{Optional. Integer. A single district to focus on.}

\item{adj}{zero indexed adjacency list.}

\item{plan}{An integer vector or matrix column of district assignments.}

\item{simplify}{Optional. Logical. Whether to return extra information or just grouping ID.}
}
\value{
integer vector (if simplify is false). Otherwise it returns a tibble with the grouping
variable as \code{group_id} and additional information on connected components.
}
\description{
Creates a grouping ID to unite geographies and perform analysis on a smaller
set of precincts. It identifies all precincts more than \code{boundary} edges
of a district district boundary. Each contiguous group of precincts more than
\code{boundary} steps away from another district gets it own group. Some
districts may have multiple, disconnected components that make up the core,
but each of these is assigned a separate grouping id so that a call to
\code{sf::st_union()} would produce only connected pieces.
}
\details{
This is a loose interpretation of the NCSL's summary
of redistricting criteria to preserve the cores of prior districts. Using the
adjacency graph for a given plan, it will locate the precincts on the
boundary of the district, within \code{boundary} steps of the edge. Each of
these  is given their own group. Each remaining entry that is not near the
boundary of the district is given an id that can be used to group the
remainder of the district by connected component. This portion is deemed the
core of the district.
}
\examples{
data(fl250)
fl250_map <- redist_map(fl250, ndists = 4, pop_tol = 0.01)
plan <- as.matrix(redist_smc(fl250_map, 20, silent = TRUE))
core <- redist.identify.cores(adj = fl250_map$adj, plan = plan)
redist.plot.cores(shp = fl250, plan = plan, core = core)

}
\seealso{
\code{\link[=redist.plot.cores]{redist.plot.cores()}} for a plotting function
}
\concept{prepare}
