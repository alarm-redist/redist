\name{redist-package}
\alias{redist-package}
\alias{redist}
\docType{package}
\title{R Package for the MCMC Redistricting Simulator}
\description{
  \code{redist} implements methods developed by Fifield, Higgins, Imai and
  Tarr (2016) to randomly sample congressional redistricting plans
  using Markov Chain Monte Carlo methods. The current version of this
  package implements the basic simulator and offers several
  improvements to improve the performance of the algorithm and to
  incorporate substantive constraints found in American congressional
  redistricting. First, it allows users to draw plans that are nearly
  equal in population. Second, users can apply constraints that
  increase the geographic compactness of redistricting plans. Third,
  it implements several tempering techniques to help efficiently
  explore the full distribution of redistricting plans. Finally, it
  allows users to generate standard diagnostics from the Markov Chain
  Monte Carlo literature in order to examine the performacne of the
  simulations.
  }

\details{
\tabular{ll}{
Package: \tab redist\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-08\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
  Benjamin Fifield, Department of Politics, Princeton University 
  \email{bfifield@princeton.edu}, \url{http://www.benfifield.com}

  Michael Higgins, Department of Statistics, Princeton University
  \email{mikehiggins@k-state.edu}, \url{http://www-personal.ksu.edu/~mikehiggins/}
  
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@princeton.edu}, \url{http://imai.princeton.edu}
  
  Alexander Tarr, Department of Electrical Engineering, Princeton University
  \email{atarr@princeton.edu}
  
  Maintainer: Ben Fifield <bfifield@princeton.edu>
}

\references{
Barbu, Adrian and Song-Chun Zhu. (2005) "Generalizing Swendsen-Wang to
Sampling Arbitrary Posterior Probabilities." IEEE Transactions on Pattern Analysis and
Machine Intelligence.

  Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
  Tarr. (2016) "A New Automated Redistricting Simulator Using Markov
  Chain Monte Carlo."   
  \emph{Working Paper.} Available at
  \url{http://imai.princeton.edu/research/files/redist.pdf}.

Swendsen, Robert and Jian-Sheng Wang. (1987) "Nonuniversal Critical
Dynamics in Monte Carlo Simulations." Physical Review Letters.
}
\keyword{package}

