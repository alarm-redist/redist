\name{redist.enumerate}
\alias{redist.enumerate}
\title{
Exact Redistricting Plan Enumerator
}
\description{
\code{redist.enumerate} uses a spanning-tree method to fully enumerate all valid redistricting plans with $n$ districts given a set of geographic units. \code{redist.enumerate} also allows suers to implement minimum and maximum numbers of geographic units per district, as well as population parity requirements.
}
\usage{
redist.enumerate(adjobj,
                 ndists = 2,
                 popvec = NULL,
                 nconstraintlow = NULL,
                 nconstrainthigh = NULL,
                 popcons = NULL,
		 contiguitymap = "rooks")
}
\arguments{
  \item{adjobj}{An adjacency list, matrix, or object of class \code{SpatialPolygonsDataFrame}.}
  \item{ndists}{The desired number of congressional districts. The default is 2.}
  \item{popvec}{A vector of geographic unit populations. The default is \code{NULL}.}
  \item{nconstraintlow}{Lower bound for number of geographic units to include in a district. The default is \code{NULL}.}
  \item{nconstrainthigh}{Lower bound for number of geographic units to include in a district. The default is \code{NULL}.}
  \item{popcons}{The strength of the hard population constraint. \code{popcons} = 0.05 means that any proposed swap that brings a district more than 5\% away from population parity will be rejected. The default is \code{NULL}.}
  \item{contiguitymap}{Use queens or rooks distance criteria for generating an adjacency list from a "SpatialPolygonsDataFrame" data type. Default is "rooks".} 
}
\details{
This function allows users to input a set of geographic units to generate all valid partitions of $n$ congressional districts. The function uses a set of spanning-tree methods to generate all valid, contiguous partitions, which makes it more efficient than brute-force methods. However, even with these methods, full redistricting problems quickly become intractable, necessitating the use of the MCMC-based methods implemented in \code{redist.mcmc}.
}
\value{
\code{redist.enumerate} returns an object of class "list". Each entry in the list is a vector of congressional district assignments, where the first entry in the vector corresponds to the congressional district assignment of the first geographic unit.
}
\references{
  Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr. (2015) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo." Working 
  Paper. Available at \url{http://imai.princeton.edu/research/files/redist.pdf}.
}

\examples{
\dontrun{
data(algdat.pfull)
test <- redist.enumerate(adjobj = algdat.pfull$adjlist)
}
}

\seealso{
\code{\link{redist.mcmc}} for scaleable Markov chain Monte Carlo methods for redistricting simulation.}

\keyword{tree}
\keyword{spanningtree}
\keyword{enumerate}
