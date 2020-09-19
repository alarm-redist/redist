#' Markov Chain Monte Carlo Methods for Redistricting Simulation
#'
#' Enables researchers to sample redistricting plans from a pre-specified
#' target distribution using a Markov Chain Monte Carlo algorithm. The package
#' allows for the implementation of various constraints in the redistricting
#' process such as geographic compactness and population parity requirements.
#' The algorithm also can be used in combination with efficient simulation
#' methods such as simulated and parallel tempering algorithms. Tools for
#' analysis such as inverse probability reweighting and plotting functionality
#' are included. The package implements methods described in Fifield, Higgins,
#' Imai and Tarr (2016) ``A New Automated Redistricting Simulator Using Markov
#' Chain Monte Carlo,'' working paper available at
#' <http://imai.fas.harvard.edu/research/files/redist.pdf>.
#'
#' \tabular{ll}{ Package: \tab redist\cr Type: \tab Package\cr
#' Version: \tab 2.0.1\cr
#' Date: \tab 2020-09-19\cr License: \tab GPL (>= 2)\cr }
#'
#' @name redist-package
#' @useDynLib redist, .registration = TRUE
#' @aliases redist-package redist
#' @docType package
#' @author Benjamin Fifield, Department of Politics, Princeton University
#' \email{benfifield@@gmail.com}, \url{https://www.benfifield.com}
#'
#' Michael Higgins, Department of Statistics, Princeton University
#' \email{mikehiggins@@k-state.edu}, \url{http://www-personal.k-state.edu/~mikehiggins/}
#'
#' Alexander Tarr, Department of Electrical Engineering, Princeton University
#' \email{atarr@@princeton.edu}
#'
#' Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@princeton.edu}, \url{http://imai.fas.harvard.edu}
#'
#' Maintainer: Ben Fifield <benfifield@gmail.com>
#'
#' @references
#' Barbu, Adrian and Song-Chun Zhu. (2005) "Generalizing Swendsen-Wang to
#' Sampling Arbitrary Posterior Probabilities." IEEE Transactions on
#' Pattern Analysis and Machine Intelligence.
#'
#' Fifield, Benjamin, Michael Higgins, Kosuke Imai and Alexander
#' Tarr. (2016) "A New Automated Redistricting Simulator Using Markov
#' Chain Monte Carlo." \emph{Working Paper.} Available at
#' \url{http://imai.princeton.edu/research/files/redist.pdf}.
#'
#' Swendsen, Robert and Jian-Sheng Wang. (1987) "Nonuniversal Critical
#' Dynamics in Monte Carlo Simulations." Physical Review Letters.
#'
#' @keywords package
#' @importFrom Rcpp evalCpp
#' @importFrom spdep poly2nb
#' @importFrom sp coordinates
#' @importFrom coda mcmc mcmc.list traceplot autocorr.plot densplot cumuplot gelman.plot
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom grDevices dev.off pdf
#' @importFrom stats median runif quantile
NULL
