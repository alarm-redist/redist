\name{algdat.pfull}
\alias{algdat.pfull}
\docType{data}
\title{
All Partitions of 25 Precincts into 3 Congressional Districts (No Population Constraint)
}
\description{
This data set contains demographic and geographic information about 25 contiguous precincts in the state of Florida. The data lists all possible partitions of the 25 precincts into three contiguous congressional districts. 
}
\usage{data("algdat.pfull")}
\format{
  A list with five entries.
  \describe{
    \item{\code{adjlist}}{An adjacency list for the 25 precincts.}
    \item{\code{cdmat}}{A matrix containing every partition of the 25 precincts into three contiguous congressional districts, with no population constraint.}
    \item{\code{precinct.data}}{A matrix containing demographic information for each of the 25 precincts.}
    \item{\code{segregation.index}}{A matrix containing the dissimilarity index of segregation (Massey and Denton 1987) for each congressional district map in \code{cdmat}.}
    \item{\code{distancemat}}{A square matrix containing the squared distance between the centroids of any two precincts.}
  }
}

\references{
  Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr. (2015) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo." Working 
  Paper. Available at \url{http://imai.princeton.edu/research/files/redist.pdf}.
  
  Massey, Douglas and Nancy Denton. (1987) "The Dimensions of Social Segregation". Social Forces.
  }
\examples{
\dontrun{
data(algdat.pfull)
}}
\keyword{datasets}
