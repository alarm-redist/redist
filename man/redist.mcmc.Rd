\name{redist.mcmc} 
\alias{redist.mcmc}
\title{MCMC Redistricting Simulator} 
\description{ 
\code{redist.mcmc} is used to simulate Congressional redistricting
plans using Markov Chain Monte Carlo methods. 
} 

\usage{
  redist.mcmc(adjobj, popvec, nsims, ndists = NULL, initcds = NULL, loopscompleted = 0,
              nloop = 1, nthin = 1, eprob = 0.05, lambda = 0,
              popcons = NULL, grouppopvec = NULL, ssdmat = NULL,
              betacompact = 0, betapop = 0,
              betaseg = 0, betasimilar = 0,
              temperbetacompact = 0, temperbetapop = 0,
              temperbetaseg = 0, temperbetasimilar = 0,
              betaseq = "powerlaw", betaseqlength = 10,
              betaweights = NULL,
              adjswaps = TRUE, rngseed = NULL, maxiterrsg = 5000,
              savename = NULL, verbose = TRUE)
}

\arguments{
\item{adjobj}{An adjacency matrix, list, or object of class
"SpatialPolygonsDataFrame."}
\item{popvec}{A vector containing the populations of each geographic
unit.}
\item{nsims}{The number of simulations run before a save point.}
\item{ndists}{The number of congressional districts. The default is \code{NULL}.}
\item{initcds}{A vector containing the congressional district labels
of each geographic unit. The default is \code{NULL}. If not provided, random and contiguous
congressional district assignments will be generated using \code{redist.rsg}.}
\item{loopscompleted}{Number of save points reached by the
algorithm. The default is \code{0}.}
\item{nloop}{The total number of save points for the algorithm. The
default is \code{1}. Note that the total number of simulations run
will be \code{nsims} * \code{nloop}.}
\item{nthin}{The amount by which to thin the Markov Chain. The default
is \code{1}.}
\item{eprob}{The probability of keeping an edge connected. The default
is \code{0.05}.}
\item{lambda}{The parameter determining the number of swaps to attempt
each iteration of the algorithm. The number of swaps each iteration is
equal to Pois(\code{lambda}) + 1. The default is \code{0}.}
\item{popcons}{The strength of the hard population
constraint. \code{popcons} = 0.05 means that any proposed swap that
brings a district more than 5\% away from population parity will be
rejected. The default is \code{NULL}.}
\item{grouppopvec}{A vector of populations for some sub-group of
interest. The default is \code{NULL}.}
\item{ssdmat}{A matrix of squared distances between geographic
units. The default is \code{NULL}.}
\item{betacompact}{Strength of the constraint to impose geographic
compactness. Values of \code{betacompact} less than zero impose more
compactness. The default is \code{0} (no constraint).}
\item{betapop}{Strength of the constraint to impose population
parity. Values of \code{betapop} less than zero bring the simulated
districts closer to population parity. The default is \code{0} (no constraint).}
\item{betaseg}{Strength of the constraint to "pack" a particular
subgroup into a single district. Values of \code{betaseg} greater than
zero will create districts with higher concentrations of the
subgroup. The default is \code{0} (no constraint).}
\item{betasimilar}{Strength of the constraint to simulate redistricting
similar to a specified plan. Values of \code{betasimilar} less than
zero will create districts more similar to the specified plan. The
default is \code{0} (no constraint).}
\item{temperbetacompact}{Flag to apply simulated tempering to the
\code{betacompact} constraint. The default is \code{0} (no tempering).}
\item{temperbetapop}{Flag to apply simulated tempering to the
\code{betacompact} constraint. The default is \code{0} (no
tempering).}
\item{temperbetaseg}{Flag to apply simulated tempering to the
\code{betaseg} constraint. The default is \code{0} (no tempering).}
\item{temperbetasimilar}{Flag to apply simulated tempering to the
\code{betasimilar} constraint. The default is \code{0} (no
tempering).}
\item{betaseq}{Sequence of beta values for tempering. The default is
\code{powerlaw} (see Fifield et. al (2015) for details).}
\item{betaseqlength}{Length of beta sequence desired for
tempering. The default is \code{10}.}
\item{betaweights}{Sequence of weights for different values of
beta. Allows the user to upweight certain values of beta over
others. The default is \code{NULL} (equal weighting).}
\item{adjswaps}{Flag to restrict swaps of beta so that only
values adjacent to current constraint are proposed. The default is
\code{TRUE}.}
\item{rngseed}{Allows the user to set the seed for the
simulations. Default is \code{NULL}.}
\item{maxiterrsg}{Maximum number of iterations for random seed-and-grow algorithm to generate starting values. Default is 5000.}
\item{savename}{Filename to save simulations. Default is \code{NULL}.}
\item{verbose}{Whether to print initialization statement. Default is \code{TRUE}.}
}

\details{
This function allows users to simulate redistricting plans
using Markov Chain Monte Carlo methods. Several constraints
correspoding to substantive requirements in the redistricting process
are implemented, including population parity and geographic
compactness. In addition, the function includes multiple-swap and
simulated tempering functionality to improve the mixing of the Markov
Chain.
} 


\value{
\code{redist.mcmc} returns an object of class "redist". The object
\code{redist} is a list that contains the folowing components (the
inclusion of some components is dependent on whether tempering
techniques are used):

\item{partitions}{Matrix of congressional district assignments generated by the algorithm. Each row corresponds to a geographic unit, and each column corresponds to a simulation.}
\item{distance_parity}{Vector containing the maximum distance from parity for a particular simulated redistricting plan.}
\item{mhdecisions}{A vector specifying whether a proposed redistricting plan was accepted (1) or rejected (0) in a given iteration.}
\item{mhprob}{A vector containing the Metropolis-Hastings acceptance probability for each iteration of the algorithm.}
\item{pparam}{A vector containing the draw of the \code{p} parameter for each simulation, which dictates the number of swaps attempted.}
\item{constraint_pop}{A vector containing the value of the population constraint for each accepted redistricting plan.}
\item{constraint_compact}{A vector containing the value of the compactness constraint for each accepted redistricting plan.}
\item{constraint_segregation}{A vector containing the value of the segregation constraint for each accepted redistricting plan.}
\item{constraint_similar}{A vector containing the value of the similarity constraint for each accepted redistricting plan.}
\item{beta_sequence}{A vector containing the value of beta for each iteration of the algorithm. Returned when tempering is being used.}
\item{mhdecisions_beta}{A vector specifying whether a proposed beta value was accepted (1) or rejected (0) in a given iteration of the algorithm. Returned when tempering is being used.}
\item{mhprob_beta}{A vector containing the Metropolis-Hastings acceptance probability for each iteration of the algorithm. Returned when tempering is being used.}

}

\references{
  Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander Tarr. (2015) "A New Automated Redistricting Simulator Using Markov Chain Monte Carlo." Working 
  Paper. Available at \url{http://imai.princeton.edu/research/files/redist.pdf}.
}

\examples{
\dontrun{
data(algdat.pfull)

## Code to run the simulations in Figure 4 in Fifield, Higgins, Imai and Tarr (2015)

## Get an initial partition
set.seed(1)
initcds <- algdat.pfull$cdmat[,sample(1:ncol(algdat.pfull$cdmat), 1)]

## Run the algorithm
alg_253 <- redist.mcmc(adjobj = algdat.pfull$adjlist,
                       popvec = algdat.pfull$precinct.data$pop,
                       initcds = initcds,
                       nsims = 10000)

}
}

\seealso{
\code{\link{redist.ipw}} for inverse probability weighting functionality.
}

\keyword{mcmc}
